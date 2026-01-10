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
static int ccsds_codeword_length;  // 依赖码率: R=1/3 为 K*3+16, R=1/2 为 K*2+8

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
    ccsds_codeword_length = k * 3 + CCSDS_STATE_MEM * 4;
    
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
    ccsds_codeword_length = K * 3 + CCSDS_STATE_MEM * 4;
    
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
 * @brief 比特交织
 */
static void ccsds_interleave_bits(int* in, int* out) {
    int K = g_ccsds_k;
    for (int i = 0; i < K; i++) {
        out[i] = in[ccsds_interleaver[i]];
    }
    // 尾比特不交织
    for (int i = K; i < ccsds_message_length; i++) {
        out[i] = in[i];
    }
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
// G1 = 11011 (binary) = 1 + D + D^2 + D^4 (feedforward)
//
// State: s = (s0, s1, s2, s3) where s0 is D, s1 is D^2, s2 is D^3, s3 is D^4
// Feedback: f = u ⊕ s0 ⊕ s3   (from G0 = 1 + D + D^4)
// Parity:   c = f ⊕ s0 ⊕ s1 ⊕ s3  (from G1 = 1 + D + D^2 + D^4)
// State update: s3=s2, s2=s1, s1=s0, s0=f

/**
 * @brief CCSDS 16-state 分量 RSC 编码器
 */
static void ccsds_component_rsc_encoder(int* input_msg, int* output_parity, int length) {
    int s0 = 0, s1 = 0, s2 = 0, s3 = 0;  // 4-bit state register
    
    for (int i = 0; i < length; i++) {
        // Feedback: f = u ⊕ s0 ⊕ s3 (from G0 = 1 + D + D^4)
        int f = input_msg[i] ^ s0 ^ s3;
        
        // Parity: c = f ⊕ s0 ⊕ s1 ⊕ s3 (from G1 = 1 + D + D^2 + D^4)
        output_parity[i] = f ^ s0 ^ s1 ^ s3;
        
        // State update
        s3 = s2;
        s2 = s1;
        s1 = s0;
        s0 = f;
    }
}

/**
 * @brief CCSDS Turbo 编码器 (PCCC, 支持 R=1/3 和 R=1/2 打孔)
 */
void ccsds_turbo_encoder(void) {
    int K = g_ccsds_k;
    static int interleaved_msg[CCSDS_message_length_max];
    static int parity1[CCSDS_message_length_max];
    static int parity2[CCSDS_message_length_max];
    
    // 1. 交织
    ccsds_interleave_bits(ccsds_message_padded, interleaved_msg);
    
    // 2. 编码器 1 (非交织) - 产生 C_a
    ccsds_component_rsc_encoder(ccsds_message_padded, parity1, ccsds_message_length);
    
    // 3. 编码器 2 (交织) - 产生 C_b
    ccsds_component_rsc_encoder(interleaved_msg, parity2, ccsds_message_length);
    
    // 4. 码字复用 (根据码率选择)
    int k = 0;
    
    if (g_ccsds_rate == CCSDS_RATE_1_2) {
        // === R=1/2 打孔模式 ===
        // CCSDS 131.0-B-5: 对于 1-indexed 位置 s
        // s 奇数: 发送 C_a[s]  (s=1,3,5... 即 0-indexed i=0,2,4...)
        // s 偶数: 发送 C_b[s]  (s=2,4,6... 即 0-indexed i=1,3,5...)
        for (int i = 0; i < K; i++) {
            ccsds_codeword[k++] = ccsds_message_padded[i]; // u (systematic)
            if (i % 2 == 0) {
                // i=0,2,4... (s=1,3,5... 奇数) -> C_a
                ccsds_codeword[k++] = parity1[i];
            } else {
                // i=1,3,5... (s=2,4,6... 偶数) -> C_b
                ccsds_codeword[k++] = parity2[i];
            }
        }
        // 尾比特 (交替发送)
        for (int i = K; i < ccsds_message_length; i++) {
            if (i % 2 == 0) {
                ccsds_codeword[k++] = parity1[i];
            } else {
                ccsds_codeword[k++] = parity2[i];
            }
        }
    } else {
        // === R=1/3 无打孔模式 ===
        for (int i = 0; i < K; i++) {
            ccsds_codeword[k++] = ccsds_message_padded[i]; // u (systematic)
            ccsds_codeword[k++] = parity1[i];              // C_a
            ccsds_codeword[k++] = parity2[i];              // C_b
        }
        // 尾比特的校验位
        for (int i = K; i < ccsds_message_length; i++) {
            ccsds_codeword[k++] = parity1[i];
            ccsds_codeword[k++] = parity2[i];
        }
    }
    
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
        ccsds_message_padded[i] = ccsds_message[i];
    }
    // 补零 (trellis termination)
    for (int i = K; i < ccsds_message_length; i++) {
        ccsds_message_padded[i] = 0;
    }
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
            
            // Feedback: f = u ⊕ s0 ⊕ s3
            int f = input ^ s0 ^ s3;
            
            // Parity: c = f ⊕ s0 ⊕ s1 ⊕ s3
            int parity = f ^ s0 ^ s1 ^ s3;
            
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
    for (int s = 0; s < CCSDS_ST_NUM; s++) {
        log_beta[msg_len][s] = LOG_ONE_CCSDS;  // 等概率结束
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
void ccsds_turbo_decoder(void) {
    int K = g_ccsds_k;
    
    // 1. 计算信道 LLR (根据码率进行去打孔)
    double Lc_factor = 4.0 / N0;
    int k = 0;
    
    if (g_ccsds_rate == CCSDS_RATE_1_2) {
        // === R=1/2 去打孔模式 ===
        // 接收序列: [u, parity]
        // i=0,2,4... (s=1,3,5... 奇数) -> parity 是 C_a
        // i=1,3,5... (s=2,4,6... 偶数) -> parity 是 C_b
        for (int i = 0; i < K; i++) {
            ccsds_Lc_sys[i] = Lc_factor * ccsds_rx_symbol[k++][0];
            double parity_llr = Lc_factor * ccsds_rx_symbol[k++][0];
            
            if (i % 2 == 0) {
                // i=0,2,4... (s 奇数): 接收的是 C_a, C_b 被打孔
                ccsds_Lc_par1[i] = parity_llr;  // C_a present
                ccsds_Lc_par2[i] = 0.0;         // C_b punctured
            } else {
                // i=1,3,5... (s 偶数): 接收的是 C_b, C_a 被打孔
                ccsds_Lc_par1[i] = 0.0;         // C_a punctured
                ccsds_Lc_par2[i] = parity_llr;  // C_b present
            }
        }
        // 尾比特部分
        for (int i = K; i < ccsds_message_length; i++) {
            ccsds_Lc_sys[i] = 500.0;  // 强先验 (已知为0)
            double parity_llr = Lc_factor * ccsds_rx_symbol[k++][0];
            if (i % 2 == 0) {
                ccsds_Lc_par1[i] = parity_llr;
                ccsds_Lc_par2[i] = 0.0;
            } else {
                ccsds_Lc_par1[i] = 0.0;
                ccsds_Lc_par2[i] = parity_llr;
            }
        }
    } else {
        // === R=1/3 无打孔模式 ===
        for (int i = 0; i < K; i++) {
            ccsds_Lc_sys[i]  = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_par1[i] = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_par2[i] = Lc_factor * ccsds_rx_symbol[k++][0];
        }
        // 尾比特部分
        for (int i = K; i < ccsds_message_length; i++) {
            double rx_llr = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_sys[i] = rx_llr + 500.0;  // 强先验 (已知为0)
            ccsds_Lc_par1[i] = Lc_factor * ccsds_rx_symbol[k++][0];
            ccsds_Lc_par2[i] = 0.0;
        }
    }
    
    // 2. 初始化先验信息
    for (int i = 0; i < ccsds_message_length; i++) {
        ccsds_La_2_to_1[i] = 0.0;
        ccsds_La_1_to_2[i] = 0.0;
    }
    
    // 3. 迭代译码
    static double Lc_sys_interleaved[CCSDS_message_length_max];
    
    for (int iter = 0; iter < CCSDS_ITERATIONS; iter++) {
        // DEC1: 使用原始系统位 + C_a (parity1)
        ccsds_log_map_decoder(ccsds_Lc_sys, ccsds_Lc_par1, ccsds_La_2_to_1, ccsds_Le_1);
        
        // 交织外在信息
        ccsds_interleave_llr(ccsds_Le_1, ccsds_La_1_to_2);
        
        // 交织系统位
        ccsds_interleave_llr(ccsds_Lc_sys, Lc_sys_interleaved);
        
        // DEC2: 使用交织后系统位 + C_b (parity2)
        // 注意: parity2 是从交织后的消息编码得到的, 所以不需要额外交织
        ccsds_log_map_decoder(Lc_sys_interleaved, ccsds_Lc_par2, ccsds_La_1_to_2, ccsds_Le_2);
        
        // 解交织外在信息
        ccsds_deinterleave_llr(ccsds_Le_2, ccsds_La_2_to_1);
    }
    
    // 4. 最终判决
    for (int i = 0; i < K; i++) {
        ccsds_L_APP[i] = ccsds_Lc_sys[i] + ccsds_La_2_to_1[i] + ccsds_Le_1[i];
        ccsds_de_message[i] = (ccsds_L_APP[i] > 0.0) ? 0 : 1;
    }
}

// =================================================================
// --- CCSDS 仿真主循环 ---
// =================================================================

void run_ccsds_turbo_simulation(SimConfig* cfg, FILE* csv_fp) {
    int K = g_ccsds_k;
    double code_rate = (double)K / (double)(K * 3);
    
    long int bit_error, frame_error, seq;
    double BER, FER;
    
    const long MIN_FRAME_ERRORS = 100;
    
    // 初始化
    ccsds_init_trellis();
    ccsds_turbo_init();
    
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
