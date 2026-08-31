/**
 * @file ccsds_turbo.c
 * @brief CCSDS 131.0-B-5 Turbo Code Implementation (NASA Standard)
 * 
 * 16-state RSC encoder, Dynamic K support (1784/3568/7136/8920), Rate 1/3
 * Polynomials: G0=10011 (feedback), G1=11011 (feedforward)
 * 
 * Interleaver: Exact CCSDS 131.0-B-5 algorithmic permutation
 */

#include "ccsds_turbo.h"
#include "config.h"
#include "csv_export.h"
#include "sim_runner.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// =================================================================
// --- 运行时 K 值和码率 (全局变量) ---
// =================================================================
int g_ccsds_k = CCSDS_K_1784;          // 默认 K=1784
CcsdsCodeRate g_ccsds_rate = CCSDS_RATE_1_3;  // 默认 R=1/3 (无打孔)

// 运行时计算的派生值
static int ccsds_message_length;   // K + 4
static int ccsds_codeword_length;  // 标准规定 n=(K+4)/r

// =================================================================
// --- CCSDS 全局数组 (使用最大尺寸静态分配) ---
// =================================================================
static int ccsds_message[CCSDS_K_MAX];                       // 原始消息
static int ccsds_message_padded[CCSDS_message_length_max];   // 补零后的消息
static int ccsds_interleaver[CCSDS_K_MAX];                   // 交织器: interleaver[i] = π⁻¹(i)
static int ccsds_deinterleaver[CCSDS_K_MAX];                 // 解交织器: deinterleaver[i] = π(i)
static int ccsds_codeword[CCSDS_codeword_length_max];        // 编码后的码字
static int ccsds_de_message[CCSDS_K_MAX];                    // 最终译码消息

static double ccsds_tx_symbol[CCSDS_codeword_length_max][2]; // 发送符号
static double ccsds_rx_symbol[CCSDS_codeword_length_max][2]; // 接收符号

// LLR 数组
static double ccsds_Lc_sys[CCSDS_message_length_max];
static double ccsds_Lc_par1[CCSDS_message_length_max];
static double ccsds_Lc_par2[CCSDS_message_length_max];
static double ccsds_La_1_to_2[CCSDS_message_length_max];
static double ccsds_La_2_to_1[CCSDS_message_length_max];
static double ccsds_Le_1[CCSDS_message_length_max];
static double ccsds_Le_2[CCSDS_message_length_max];
static double ccsds_L_APP[CCSDS_message_length_max];

// 信道参数
extern double N0;
extern double sgm;

// =================================================================
// --- CCSDS 交织器参数 (131.0-B-5 标准) ---
// =================================================================
// CCSDS 标准使用的 8 个质数 (1-indexed in algorithm, 0-indexed here for array access)
static const int CCSDS_PRIMES[9] = {0, 31, 37, 43, 47, 53, 59, 61, 67};

// K 值对应的 (k1, k2) 参数表
typedef struct {
    int k;
    int k1;
    int k2;
} CcsdsKParams;

static const CcsdsKParams CCSDS_K_TABLE[] = {
    {1784, 8, 223},   // 223 = 223 * 1
    {3568, 8, 446},   // 446 = 223 * 2
    {7136, 8, 892},   // 892 = 223 * 4
    {8920, 8, 1115},  // 1115 = 223 * 5
};
#define CCSDS_K_TABLE_SIZE (sizeof(CCSDS_K_TABLE) / sizeof(CCSDS_K_TABLE[0]))

double ccsds_get_nominal_rate(void) {
    switch (g_ccsds_rate) {
        case CCSDS_RATE_1_2: return 1.0 / 2.0;
        case CCSDS_RATE_1_4: return 1.0 / 4.0;
        case CCSDS_RATE_1_3:
        default:             return 1.0 / 3.0;
    }
}

int ccsds_get_codeword_length(void) {
    int symbol_count = g_ccsds_k + CCSDS_STATE_MEM;
    switch (g_ccsds_rate) {
        case CCSDS_RATE_1_2: return symbol_count * 2;
        case CCSDS_RATE_1_4: return symbol_count * 4;
        case CCSDS_RATE_1_3:
        default:             return symbol_count * 3;
    }
}

double ccsds_get_effective_rate(void) {
    return (double)g_ccsds_k / (double)ccsds_get_codeword_length();
}

/**
 * @brief 获取 K 对应的 k1, k2 参数
 */
static int ccsds_get_k_params(int k, int* k1, int* k2) {
    for (int i = 0; i < CCSDS_K_TABLE_SIZE; i++) {
        if (CCSDS_K_TABLE[i].k == k) {
            *k1 = CCSDS_K_TABLE[i].k1;
            *k2 = CCSDS_K_TABLE[i].k2;
            return 1;
        }
    }
    return 0;  // 未找到
}

/**
 * @brief 设置 CCSDS K 值
 * @param k 信息块长度 (1784, 3568, 7136, 8920)
 * @return 1 成功, 0 失败
 */
int ccsds_set_block_size(int k) {
    int k1, k2;
    if (!ccsds_get_k_params(k, &k1, &k2)) {
        fprintf(stderr, "[CCSDS] Error: Unsupported K=%d\n", k);
        return 0;
    }
    
    g_ccsds_k = k;
    ccsds_message_length = k + CCSDS_STATE_MEM;
    ccsds_codeword_length = ccsds_get_codeword_length();
    
    printf("[CCSDS] Block size set: K=%d (k1=%d, k2=%d)\n", k, k1, k2);
    return 1;
}

/**
 * @brief 生成 CCSDS 131.0-B-5 确定性交织器
 * 
 * 算法来自 CCSDS 131.0-B-5 Section 6.3g:
 * 对于 s = 1, 2, ..., K:
 *   m = (s - 1) mod 2
 *   i = floor((s - 1) / (2 * k2))
 *   j = floor((s - 1) / 2) - i * k2
 *   t = (19 * i + 1) mod (k1 / 2)
 *   q = (t mod 8) + 1
 *   c = (p[q] * j + 21 * m) mod k2
 *   pi(s) = 2 * (t + c * (k1 / 2) + 1) - m
 */
void ccsds_turbo_init(void) {
    int K = g_ccsds_k;
    int k1, k2;
    
    // 更新派生值
    ccsds_message_length = K + CCSDS_STATE_MEM;
    ccsds_codeword_length = ccsds_get_codeword_length();
    
    if (!ccsds_get_k_params(K, &k1, &k2)) {
        fprintf(stderr, "[CCSDS] Error: Cannot initialize interleaver for K=%d\n", K);
        return;
    }
    
    int k1_half = k1 / 2;  // k1/2 = 4
    
    // 生成交织表 (CCSDS 131.0-B-5 算法)
    for (int s = 1; s <= K; s++) {
        int m = (s - 1) % 2;
        int i = (s - 1) / (2 * k2);
        int j = ((s - 1) / 2) - (i * k2);
        
        int t = (19 * i + 1) % k1_half;
        int q = (t % 8) + 1;  // 1-indexed prime selector
        int c = (CCSDS_PRIMES[q] * j + 21 * m) % k2;
        
        int pi_s = 2 * (t + c * k1_half + 1) - m;
        
        // 存储 (0-based indexing for C arrays)
        // interleaver[s-1] = pi_s - 1 表示输出位置 s 对应的输入位置
        // 即 interleaved[s-1] = original[interleaver[s-1]]
        ccsds_interleaver[s - 1] = pi_s - 1;
    }
    
    // 构建 deinterleaver (反向映射)
    // deinterleaver[i] = j 表示原始位置 i 映射到交织位置 j
    // 即 interleaved[j] = original[i], 所以 j 满足 interleaver[j] = i
    for (int j = 0; j < K; j++) {
        ccsds_deinterleaver[ccsds_interleaver[j]] = j;
    }
    
    // 验证测试点
    #ifdef _DEBUG
    printf("[CCSDS] Interleaver verification:\n");
    printf("  pi(1) = %d (expected: 4, 0-based: 3)\n", ccsds_interleaver[0] + 1);
    printf("  pi(2) = %d (expected: 171, 0-based: 170)\n", ccsds_interleaver[1] + 1);
    if (K == 8920) {
        printf("  pi(8920) = %d (expected: 8749, 0-based: 8748)\n", ccsds_interleaver[8919] + 1);
    }
    #endif
    
    printf("[CCSDS] Interleaver initialized (K=%d, CCSDS 131.0-B-5 standard)\n", K);
}

/**
 * @brief LLR 交织
 */
static void ccsds_interleave_llr(double* in, double* out) {
    int K = g_ccsds_k;
    for (int i = 0; i < K; i++) {
        out[i] = in[ccsds_interleaver[i]];
    }
    for (int i = K; i < ccsds_message_length; i++) {
        out[i] = in[i];
    }
}

/**
 * @brief LLR 解交织
 */
static void ccsds_deinterleave_llr(double* in, double* out) {
    int K = g_ccsds_k;
    for (int i = 0; i < K; i++) {
        out[ccsds_interleaver[i]] = in[i];
    }
    for (int i = K; i < ccsds_message_length; i++) {
        out[i] = in[i];
    }
}

// =================================================================
// --- CCSDS 16-State RSC Encoder ---
// =================================================================
// CCSDS Polynomials:
// G0 = 10011 (binary) = 1 + D + D^4 (feedback)
// G1 = 11011 (binary) = 1 + D + D^3 + D^4 (feedforward)
//
// State: s0 is the newest bit and s3 is the oldest bit in the 4-stage register.
// With this state orientation, the CCSDS/AFF3CT connection vectors map to:
// Feedback: f = u ⊕ s2 ⊕ s3   (G0 = 10011, denominator)
// Parity:   c = f ⊕ s0 ⊕ s2 ⊕ s3  (G1 = 11011, numerator)
// State update: s3=s2, s2=s1, s1=s0, s0=f

/**
 * @brief CCSDS 16-state 分量 RSC 编码器
 */
/**
 * @brief CCSDS 16-state 分量 RSC 编码器 (带终止功能)
 * 
 * @param input_msg 输入位数组 (长度 K)
 * @param output_systematic 输出系统位数组 (长度 K+4)
 * @param output_parity 输出校验位数组 (长度 K+4)
 * @param length 输入长度 K
 */
static void ccsds_component_rsc_encoder(const int* input_msg,
                                        int* output_systematic,
                                        int* output_parity,
                                        int length) {
    int s0 = 0, s1 = 0, s2 = 0, s3 = 0;  // 4-bit state register
    
    // 1. 正常编码阶段
    for (int i = 0; i < length; i++) {
        output_systematic[i] = input_msg[i];
        // Feedback: f = u ⊕ s2 ⊕ s3 (CCSDS/AFF3CT register orientation)
        int f = input_msg[i] ^ s2 ^ s3;
        
        // Parity: c = f ⊕ s0 ⊕ s2 ⊕ s3
        output_parity[i] = f ^ s0 ^ s2 ^ s3;
        
        // State update
        s3 = s2;
        s2 = s1;
        s1 = s0;
        s0 = f;
    }

    // 2. 终止阶段 (Termination) - 4 steps
    // 目标: 使寄存器归零。
    // 对于反馈结构，令输入 u = s2 ⊕ s3，则反馈 f = 0。
    // 这样移入的就是 0，4次后状态全0。
    // 此时产生的 u 即为 "Switch Bits" (Tail Systematic Bits)
    for (int i = 0; i < CCSDS_STATE_MEM; i++) {
        // 计算 switch bit，使得 f = 0
        // f = u ^ s2 ^ s3  => 令 f=0 => u = s2 ^ s3
        int u_switch = s2 ^ s3;
        int f = 0; // By definition

        // 保存 Switch Bit (作为 systematic tail)
        output_systematic[length + i] = u_switch;

        // Parity: c = f ⊕ s0 ⊕ s2 ⊕ s3  (f=0)
        output_parity[length + i] = s0 ^ s2 ^ s3;

        // State update
        s3 = s2;
        s2 = s1;
        s1 = s0;
        s0 = f; // 0
    }
}

/** Build the CCSDS codeword by applying the standard multiplexer/puncturer. */
static int ccsds_multiplex_codeword(const int* systematic_a,
                                    const int* parity_a,
                                    const int* parity_b,
                                    int symbol_count,
                                    CcsdsCodeRate rate,
                                    int* output) {
    int out = 0;

    if (rate == CCSDS_RATE_1_2) {
        /* CCSDS 131.0-B-5 6.3j: (out0a,out1a,out0a,out1b), repeated. */
        for (int t = 0; t < symbol_count; t++) {
            output[out++] = systematic_a[t];
            output[out++] = (t % 2 == 0) ? parity_a[t] : parity_b[t];
        }
    } else {
        /* R=1/3: (out0a,out1a,out1b), repeated for all K+4 bit times. */
        for (int t = 0; t < symbol_count; t++) {
            output[out++] = systematic_a[t];
            output[out++] = parity_a[t];
            output[out++] = parity_b[t];
        }
    }

    return out;
}

/**
 * @brief CCSDS Turbo 编码器 (PCCC, 支持 R=1/3 和 R=1/2 打孔)
 */
/**
 * @brief CCSDS Turbo 编码器 (PCCC, 支持 R=1/3 和 R=1/2 打孔, 均含 Trellis Termination)
 */
void ccsds_turbo_encoder(void) {
    int K = g_ccsds_k;
    static int interleaved_msg[CCSDS_K_MAX]; // 仅存 K 位
    
    static int systematic1[CCSDS_message_length_max];
    static int systematic2[CCSDS_message_length_max];
    static int parity1[CCSDS_message_length_max];
    static int parity2[CCSDS_message_length_max];

    
    // 1. 交织 (仅对 K 个信息位进行)
    // 注意: ccsds_interleave_bits 原来是对 padded (K+4) 进行的，现在只需对 input (K) 进行
    // 但为了复用原交织器表(大小 K), 我们只取前 K 个即可。
    // ccsds_message_padded 现在只存 K 位有效信息, 后面不补零用于编码输入
    for (int i=0; i<K; i++) {
        interleaved_msg[i] = ccsds_message_padded[ccsds_interleaver[i]];
    }
    
    // 2. 编码器 1 
    ccsds_component_rsc_encoder(ccsds_message_padded, systematic1, parity1, K);
    
    // 3. 编码器 2 
    ccsds_component_rsc_encoder(interleaved_msg, systematic2, parity2, K);
    
    // 4. 码字复用 (Multiplexing)
    int k = ccsds_multiplex_codeword(systematic1, parity1, parity2,
                                     ccsds_message_length, g_ccsds_rate,
                                     ccsds_codeword);
    
    // 更新实际码字长度
    ccsds_codeword_length = k;
}

/**
 * @brief CCSDS Turbo BPSK 调制
 */
void ccsds_turbo_modulation(void) {
    for (int i = 0; i < ccsds_codeword_length; i++) {
        ccsds_tx_symbol[i][0] = -1.0 * (2.0 * ccsds_codeword[i] - 1.0); // 0->+1, 1->-1
        ccsds_tx_symbol[i][1] = 0.0;
    }
}

/**
 * @brief CCSDS Turbo AWGN 信道
 */
void ccsds_turbo_channel(void) {
    for (int i = 0; i < ccsds_codeword_length; i++) {
        for (int j = 0; j < 2; j++) {
            double u = (double)rand() / (double)RAND_MAX;
            if (u == 1.0) u = 0.999999;
            double r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));
            u = (double)rand() / (double)RAND_MAX;
            if (u == 1.0) u = 0.999999;
            double g = r * cos(2.0 * PI * u);
            ccsds_rx_symbol[i][j] = ccsds_tx_symbol[i][j] + g;
        }
    }
}

/**
 * @brief 随机生成 CCSDS 消息
 */
void ccsds_turbo_generate_message(void) {
    int K = g_ccsds_k;
    for (int i = 0; i < K; i++) {
        ccsds_message[i] = rand() % 2;
        ccsds_message_padded[i] = ccsds_message[i]; // 仅填充数据部分
    }
    // 注意: ccsds_message_padded[K..K+3] 将由编码器的 Termination 过程生成 (switch bits)
    // 仿真时不需要这里预设为0 (或者设为0也无妨，会被覆盖/忽略，但为了清晰，不管它)
}

/**
 * @brief 检查误比特
 */
long int ccsds_turbo_check_errors(void) {
    int K = g_ccsds_k;
    long int errors = 0;
    for (int i = 0; i < K; i++) {
        if (ccsds_message[i] != ccsds_de_message[i]) {
            errors++;
        }
    }
    return errors;
}

// =================================================================
// --- CCSDS 16-State Log-MAP Decoder ---
// =================================================================
// 16-state trellis: 32 edges (16 states × 2 inputs)

#define LOG_ZERO_CCSDS -1.0e10
#define LOG_ONE_CCSDS 0.0

// 16-state trellis structure
typedef struct {
    int input;      // 输入 (0 or 1)
    int output;     // 校验输出 (0 or 1)
    int from_state; // 起始状态 (0-15)
    int to_state;   // 终止状态 (0-15)
} CcsdsTrellisEdge;

// Generate trellis table for CCSDS 16-state RSC
// G0 = 10011, G1 = 11011
static CcsdsTrellisEdge ccsds_trellis[CCSDS_LINE_NUM];

static void ccsds_init_trellis(void) {
    int idx = 0;
    for (int state = 0; state < CCSDS_ST_NUM; state++) {
        for (int input = 0; input < 2; input++) {
            // State bits: s3 s2 s1 s0 (s0 is LSB, most recent)
            int s0 = (state >> 0) & 1;
            int s1 = (state >> 1) & 1;
            int s2 = (state >> 2) & 1;
            int s3 = (state >> 3) & 1;
            
            // Feedback: f = u ⊕ s2 ⊕ s3 (same orientation as the encoder)
            int f = input ^ s2 ^ s3;
            
            // G1=11011: c = f ⊕ s0 ⊕ s2 ⊕ s3
            int parity = f ^ s0 ^ s2 ^ s3;
            
            // Next state: shift in f
            int next_state = ((s2 << 3) | (s1 << 2) | (s0 << 1) | f) & 0xF;
            
            ccsds_trellis[idx].input = input;
            ccsds_trellis[idx].output = parity;
            ccsds_trellis[idx].from_state = state;
            ccsds_trellis[idx].to_state = next_state;
            idx++;
        }
    }
}

/**
 * @brief Jacobian Logarithm
 */
static double ccsds_jac_log(double a, double b) {
    if (a <= LOG_ZERO_CCSDS) return b;
    if (b <= LOG_ZERO_CCSDS) return a;
    
    double diff = a - b;
    if (diff > 10.0) return a;
    if (diff < -10.0) return b;
    
    return (a > b ? a : b) + log(1.0 + exp(-fabs(diff)));
}

/**
 * @brief CCSDS 16-state Log-MAP 译码器
 */
static void ccsds_log_map_decoder(double* Lc_sys, double* Lc_par, double* La_in, double* Le_out) {
    static double log_alpha[CCSDS_message_length_max + 1][CCSDS_ST_NUM];
    static double log_beta[CCSDS_message_length_max + 1][CCSDS_ST_NUM];
    static double log_gamma[CCSDS_message_length_max][CCSDS_LINE_NUM];
    
    int msg_len = ccsds_message_length;
    
    // 1. 计算分支度量 (Log-Gamma)
    for (int t = 0; t < msg_len; t++) {
        for (int idx = 0; idx < CCSDS_LINE_NUM; idx++) {
            int u = ccsds_trellis[idx].input;
            int c = ccsds_trellis[idx].output;
            
            double u_val = (u == 0) ? 1.0 : -1.0;
            double c_val = (c == 0) ? 1.0 : -1.0;
            
            log_gamma[t][idx] = 0.5 * ((Lc_sys[t] + La_in[t]) * u_val + Lc_par[t] * c_val);
        }
    }
    
    // 2. 计算前向概率 (Log-Alpha)
    for (int s = 0; s < CCSDS_ST_NUM; s++) {
        log_alpha[0][s] = LOG_ZERO_CCSDS;
    }
    log_alpha[0][0] = LOG_ONE_CCSDS;
    
    for (int t = 1; t < msg_len + 1; t++) {
        for (int s = 0; s < CCSDS_ST_NUM; s++) {
            log_alpha[t][s] = LOG_ZERO_CCSDS;
        }
        
        for (int idx = 0; idx < CCSDS_LINE_NUM; idx++) {
            int from_s = ccsds_trellis[idx].from_state;
            int to_s = ccsds_trellis[idx].to_state;
            double val = log_alpha[t-1][from_s] + log_gamma[t-1][idx];
            log_alpha[t][to_s] = ccsds_jac_log(log_alpha[t][to_s], val);
        }
        
        // Normalization
        double norm_val = log_alpha[t][0];
        for (int s = 0; s < CCSDS_ST_NUM; s++) {
            if (log_alpha[t][s] > LOG_ZERO_CCSDS)
                log_alpha[t][s] -= norm_val;
        }
    }
    
    // 3. 计算后向概率 (Log-Beta)
    // Trellis Termination: 使得最终状态必须为 0
    for (int s = 0; s < CCSDS_ST_NUM; s++) {
        if (s == 0)
            log_beta[msg_len][s] = LOG_ONE_CCSDS;  // State 0 prob = 1
        else
            log_beta[msg_len][s] = LOG_ZERO_CCSDS; // Others = 0
    }
    
    for (int t = msg_len - 1; t >= 0; t--) {
        for (int s = 0; s < CCSDS_ST_NUM; s++) {
            log_beta[t][s] = LOG_ZERO_CCSDS;
        }
        
        for (int idx = 0; idx < CCSDS_LINE_NUM; idx++) {
            int from_s = ccsds_trellis[idx].from_state;
            int to_s = ccsds_trellis[idx].to_state;
            double val = log_beta[t+1][to_s] + log_gamma[t][idx];
            log_beta[t][from_s] = ccsds_jac_log(log_beta[t][from_s], val);
        }
        
        // Normalization
        double norm_val = log_beta[t][0];
        for (int s = 0; s < CCSDS_ST_NUM; s++) {
            if (log_beta[t][s] > LOG_ZERO_CCSDS)
                log_beta[t][s] -= norm_val;
        }
    }
    
    // 4. 计算外在信息
    for (int t = 0; t < msg_len; t++) {
        double L_app_0 = LOG_ZERO_CCSDS;
        double L_app_1 = LOG_ZERO_CCSDS;
        
        for (int idx = 0; idx < CCSDS_LINE_NUM; idx++) {
            int u = ccsds_trellis[idx].input;
            int from_s = ccsds_trellis[idx].from_state;
            int to_s = ccsds_trellis[idx].to_state;
            
            double prob = log_alpha[t][from_s] + log_gamma[t][idx] + log_beta[t+1][to_s];
            
            if (u == 0) {
                L_app_0 = ccsds_jac_log(L_app_0, prob);
            } else {
                L_app_1 = ccsds_jac_log(L_app_1, prob);
            }
        }
        
        double L_all = L_app_0 - L_app_1;
        Le_out[t] = L_all - (Lc_sys[t] + La_in[t]);
        
        // Clipping
        if (Le_out[t] > 50.0) Le_out[t] = 50.0;
        if (Le_out[t] < -50.0) Le_out[t] = -50.0;
    }
}

/**
 * @brief CCSDS Turbo 迭代译码器 (支持 R=1/3 和 R=1/2 去打孔)
 */
/**
 * @brief CCSDS Turbo 迭代译码器 (支持 Trellis Termination)
 */
void ccsds_turbo_decoder(void) {
    int K = g_ccsds_k;
    int msg_len = ccsds_message_length;
    double Lc_factor = 4.0 / N0;
    
    // 1. 读取信道 LLR
    int k = 0;

    /* Start every absent (punctured) observation as an erasure. */
    for (int i = 0; i < msg_len; i++) {
        ccsds_Lc_sys[i] = 0.0;
        ccsds_Lc_par1[i] = 0.0;
        ccsds_Lc_par2[i] = 0.0;
    }

    if (g_ccsds_rate == CCSDS_RATE_1_2) {
        /* Same alternating puncturing phase is used through all four tail times. */
        for (int t = 0; t < msg_len; t++) {
            ccsds_Lc_sys[t] = Lc_factor * ccsds_rx_symbol[k++][0];
            double parity_llr = Lc_factor * ccsds_rx_symbol[k++][0];
            if (t % 2 == 0) ccsds_Lc_par1[t] = parity_llr;
            else            ccsds_Lc_par2[t] = parity_llr;
        }
    } else {
        for (int t = 0; t < msg_len; t++) {
            ccsds_Lc_sys[t]  = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_par1[t] = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_par2[t] = Lc_factor * ccsds_rx_symbol[k++][0];
        }
    }

    // 2. 初始化先验信息
    for (int i = 0; i < msg_len; i++) {
        ccsds_La_2_to_1[i] = 0.0;
        ccsds_La_1_to_2[i] = 0.0;
    }
    
    // Dec2 sees interleaved data systematics; its four switch bits are not transmitted.
    static double Lc_sys_for_dec2[CCSDS_message_length_max];
    
    // 4. 迭代译码
    for (int iter = 0; iter < CCSDS_ITERATIONS; iter++) {
        // --- DEC 1 ---
        // Input: Sys[0..K+3], Par1[0..K+3], La[0..K+3]
        // 注意: La 中的 tail 部分应该始终为 0 (因为 tail bit 是独立的，无交换信息)
        for(int i=0; i<4; i++) ccsds_La_2_to_1[K+i] = 0.0;
        
        ccsds_log_map_decoder(ccsds_Lc_sys, ccsds_Lc_par1, ccsds_La_2_to_1, ccsds_Le_1);
        
        // 交织外在信息: Le_1[0..K-1] -> La_1_to_2[0..K-1]
        // Tail 部分 Le_1[K..K+3] 丢弃
        ccsds_interleave_llr(ccsds_Le_1, ccsds_La_1_to_2);
        
        // --- DEC 2 ---
        // 准备 DEC2 的系统位输入: Interleaved(Sys[0..K-1]) + Tail2
        ccsds_interleave_llr(ccsds_Lc_sys, Lc_sys_for_dec2); 
        for(int i=0; i<4; i++) Lc_sys_for_dec2[K+i] = 0.0;
        
        for(int i=0; i<4; i++) ccsds_La_1_to_2[K+i] = 0.0;

        ccsds_log_map_decoder(Lc_sys_for_dec2, ccsds_Lc_par2, ccsds_La_1_to_2, ccsds_Le_2);
        
        // 解交织外在信息: Le_2[0..K-1] -> La_2_to_1[0..K-1]
        ccsds_deinterleave_llr(ccsds_Le_2, ccsds_La_2_to_1);
    }
    
    /* Final APP comes from the last decoder that ran, then returns to original order. */
    for (int i = 0; i < K; i++) {
        ccsds_L_APP[i] = Lc_sys_for_dec2[i] + ccsds_La_1_to_2[i] + ccsds_Le_2[i];
    }
    for (int i = 0; i < K; i++) {
        double app = ccsds_L_APP[ccsds_deinterleaver[i]];
        ccsds_de_message[i] = (app > 0.0) ? 0 : 1;
    }
}

int ccsds_turbo_self_test(void) {
    static const int input[] = {1,0,1,1,0,0,1,0,1,0,1,1};
    static const int expected_parity[] = {1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,1};
    static const int expected_systematic[] = {1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1};
    int systematic[16], parity[16];
    int ok = 1;

    ccsds_component_rsc_encoder(input, systematic, parity, 12);
    if (memcmp(systematic, expected_systematic, sizeof(systematic)) != 0 ||
        memcmp(parity, expected_parity, sizeof(parity)) != 0) {
        fprintf(stderr, "[CCSDS SELF-TEST] FAIL: fixed RSC parity/tail vector\n");
        ok = 0;
    } else {
        printf("[CCSDS SELF-TEST] PASS: fixed RSC parity/tail vector\n");
    }

    {
        static const int sys[] = {1,0,1,1,1,0,0,1};
        static const int pa[]  = {0,1,1,0,1,0,1,0};
        static const int pb[]  = {1,1,0,0,0,1,0,1};
        static const int expected[] = {1,0,0,1,1,1,1,0,1,1,0,1,0,1,1,1};
        int muxed[16];
        int length = ccsds_multiplex_codeword(sys, pa, pb, 8,
                                              CCSDS_RATE_1_2, muxed);
        if (length != 16 || memcmp(muxed, expected, sizeof(expected)) != 0) {
            fprintf(stderr, "[CCSDS SELF-TEST] FAIL: R=1/2 puncturing phase/tail\n");
            ok = 0;
        } else {
            printf("[CCSDS SELF-TEST] PASS: R=1/2 puncturing phase/tail\n");
        }
    }

    {
        int saved_k = g_ccsds_k;
        CcsdsCodeRate saved_rate = g_ccsds_rate;
        double saved_n0 = N0;
        double saved_sgm = sgm;

        g_ccsds_k = CCSDS_K_1784;
        g_ccsds_rate = CCSDS_RATE_1_2;
        ccsds_init_trellis();
        ccsds_turbo_init();
        for (int i = 0; i < g_ccsds_k; i++) {
            int bit = ((i * 13 + i / 7 + 3) % 17) < 8;
            ccsds_message[i] = bit;
            ccsds_message_padded[i] = bit;
        }
        ccsds_turbo_encoder();
        if (ccsds_codeword_length != 2 * (g_ccsds_k + CCSDS_STATE_MEM)) {
            fprintf(stderr, "[CCSDS SELF-TEST] FAIL: codeword length %d (expected 3576)\n",
                    ccsds_codeword_length);
            ok = 0;
        }
        ccsds_turbo_modulation();
        for (int i = 0; i < ccsds_codeword_length; i++) {
            ccsds_rx_symbol[i][0] = ccsds_tx_symbol[i][0];
            ccsds_rx_symbol[i][1] = 0.0;
        }
        N0 = 0.25;
        sgm = sqrt(N0 / 2.0);
        ccsds_turbo_decoder();
        if (ccsds_turbo_check_errors() != 0) {
            fprintf(stderr, "[CCSDS SELF-TEST] FAIL: noiseless R=1/2 round trip\n");
            ok = 0;
        } else {
            printf("[CCSDS SELF-TEST] PASS: n=3576 and noiseless R=1/2 round trip\n");
        }

        /* The polynomial/tail fix is shared by R=1/3, so guard that path too. */
        g_ccsds_rate = CCSDS_RATE_1_3;
        ccsds_turbo_init();
        ccsds_turbo_encoder();
        if (ccsds_codeword_length != 3 * (g_ccsds_k + CCSDS_STATE_MEM)) {
            fprintf(stderr, "[CCSDS SELF-TEST] FAIL: R=1/3 codeword length %d (expected 5364)\n",
                    ccsds_codeword_length);
            ok = 0;
        }
        ccsds_turbo_modulation();
        for (int i = 0; i < ccsds_codeword_length; i++) {
            ccsds_rx_symbol[i][0] = ccsds_tx_symbol[i][0];
            ccsds_rx_symbol[i][1] = 0.0;
        }
        ccsds_turbo_decoder();
        if (ccsds_turbo_check_errors() != 0) {
            fprintf(stderr, "[CCSDS SELF-TEST] FAIL: noiseless R=1/3 round trip\n");
            ok = 0;
        } else {
            printf("[CCSDS SELF-TEST] PASS: n=5364 and noiseless R=1/3 round trip\n");
        }

        g_ccsds_k = saved_k;
        g_ccsds_rate = saved_rate;
        N0 = saved_n0;
        sgm = saved_sgm;
        ccsds_message_length = g_ccsds_k + CCSDS_STATE_MEM;
        ccsds_codeword_length = ccsds_get_codeword_length();
    }

    return ok;
}

// =================================================================
// --- CCSDS 仿真主循环 ---
// =================================================================

void run_ccsds_turbo_simulation(SimConfig* cfg, FILE* csv_fp) {
    int K = g_ccsds_k;
    double nominal_rate = ccsds_get_nominal_rate();
    double code_rate = ccsds_get_effective_rate();
    
    long int bit_error, frame_error, seq;
    double BER, FER;
    
    const long MIN_FRAME_ERRORS = 100;
    
    // 初始化
    ccsds_init_trellis();
    ccsds_turbo_init();
    printf("[CCSDS] Nominal rate %.8f, effective rate %.8f, codeword %d bits\n",
           nominal_rate, code_rate, ccsds_get_codeword_length());
    
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    printf("|  Eb/N0   |  Bit Errors   |  Total Bits   |      BER       | Frame Errors  |      FER       |\n");
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    
    for (float SNR = cfg->start_snr; SNR <= cfg->end_snr + 0.001; SNR += cfg->snr_step) {
        N0 = (1.0 / code_rate) / pow(10.0, SNR / 10.0);
        sgm = sqrt(N0 / 2.0);
        
        bit_error = 0;
        frame_error = 0;
        long actual_frames = 0;
        
        for (seq = 1; seq <= cfg->num_frames; seq++) {
            ccsds_turbo_generate_message();
            ccsds_turbo_encoder();
            ccsds_turbo_modulation();
            ccsds_turbo_channel();
            ccsds_turbo_decoder();
            
            long errors_this_frame = ccsds_turbo_check_errors();
            bit_error += errors_this_frame;
            if (errors_this_frame > 0) {
                frame_error++;
            }
            actual_frames = seq;
            
            // Early termination
            if (frame_error >= MIN_FRAME_ERRORS && seq >= 1000) {
                break;
            }
            
            // Progress
            if (seq % 500 == 0) {
                printf("\r  [SNR=%.1fdB] Frame %ld/%ld (%.1f%%) - BitErr: %ld, FrameErr: %ld   ", 
                       SNR, seq, cfg->num_frames, 
                       100.0 * seq / cfg->num_frames, bit_error, frame_error);
                fflush(stdout);
            }
        }
        printf("\r                                                                        \r");
        
        long total_bits = (long)K * actual_frames;
        BER = (double)bit_error / (double)total_bits;
        FER = (double)frame_error / (double)actual_frames;
        
        printf("|  %5.1f   |  %11ld  |  %11ld  |  %12.4e  |  %11ld  |  %12.4e  |\n", 
               SNR, bit_error, total_bits, BER, frame_error, FER);
        
        if (csv_fp) {
            csv_append_row_with_fer(csv_fp, SNR, bit_error, total_bits, BER, 
                                    frame_error, actual_frames, FER);
        }
    }
    
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
}
