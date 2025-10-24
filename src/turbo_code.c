#include "turbo_code.h"
#include "log_map_decoder.h"
#include "trellis.h" // 我们将复用 CC 的 trellis 定义
#include "config.h"

// --- Turbo 码全局数组 ---
int turbo_message[TURBO_MESSAGE_BITS]; // 原始消息
int turbo_message_padded[TURBO_message_length]; // 补零后的消息
int turbo_interleaver[TURBO_MESSAGE_BITS]; // 交织器
int turbo_codeword[TURBO_codeword_length]; // 编码后的码字 [u1, c1, c2, u2, c2, c2, ...]
int turbo_de_message[TURBO_MESSAGE_BITS]; // 最终译码消息

double turbo_tx_symbol[TURBO_codeword_length][2]; // 发送符号
double turbo_rx_symbol[TURBO_codeword_length][2]; // 接收符号

// LLR 数组
double Lc_sys[TURBO_message_length];  // 信道 LLR (系统比特)
double Lc_par1[TURBO_message_length]; // 信道 LLR (校验比特1)
double Lc_par2[TURBO_message_length]; // 信道 LLR (校验比特2)
double La_1_to_2[TURBO_message_length]; // 解码器1 -> 2 的先验信息
double La_2_to_1[TURBO_message_length]; // 解码器2 -> 1 的先验信息
double Le_1[TURBO_message_length];    // 解码器1 的外在信息
double Le_2[TURBO_message_length];    // 解码器2 的外在信息
double L_APP[TURBO_message_length];   // 最终的后验 LLR

// 全局信道参数 (在 Turbo 循环中设置)
extern double N0;
extern double sgm;

/**
 * @brief 生成一个简单的随机交织器
 */
void turbo_generate_interleaver() {
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        turbo_interleaver[i] = i;
    }
    // Fisher-Yates shuffle 算法
    for (int i = TURBO_MESSAGE_BITS - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = turbo_interleaver[i];
        turbo_interleaver[i] = turbo_interleaver[j];
        turbo_interleaver[j] = temp;
    }
}

/**
 * @brief 比特交织
 */
void interleave_bits(int* in, int* out) {
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        out[i] = in[turbo_interleaver[i]];
    }
    // 尾比特不交织
    for(int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        out[i] = in[i];
    }
}

/**
 * @brief LLR 交织
 */
void interleave_llr(double* in, double* out) {
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        out[i] = in[turbo_interleaver[i]];
    }
     for(int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        out[i] = in[i];
    }
}

/**
 * @brief LLR 解交织
 */
void deinterleave_llr(double* in, double* out) {
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        out[turbo_interleaver[i]] = in[i];
    }
    for(int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        out[i] = in[i];
    }
}

/**
 * @brief 分量编码器 (RSC)
 * G = [1, (1+D^2)/(1+D+D^2)] (7,5 码的 RSC 形式)
 * 反馈: f = m[i] ^ s0 ^ s1
 * 校验: c = m[i] ^ s1 (即 G_par = 1+D^2)
 *
 * (为了匹配 trellis.c 中的 (7,5) 网格图，我们使用)
 * G = [1, (1+D^2)/(1+D+D^2)]
 * 反馈: f = m[i] ^ s0 ^ s1
 * 校验: c = m[i] ^ s1
 * s1_next = s0
 * s0_next = f
 */
void component_rsc_encoder(int* input_msg, int* output_parity) {
    int s0 = 0, s1 = 0;
    for (int i = 0; i < TURBO_message_length; i++) {
        int f = input_msg[i] ^ s0 ^ s1; // 反馈 (来自 1+D+D^2)
        output_parity[i] = input_msg[i] ^ s1; // 校验 (来自 1+D^2)
        
        s1 = s0;
        s0 = f;
    }
}

/**
 * @brief Turbo 码编码器 (PCCC, R=1/3)
 */
void turbo_encoder() {
    int interleaved_msg[TURBO_message_length];
    int parity1[TURBO_message_length];
    int parity2[TURBO_message_length];

    // 1. 交织
    interleave_bits(turbo_message_padded, interleaved_msg);

    // 2. 编码器 1 (非交织)
    component_rsc_encoder(turbo_message_padded, parity1);

    // 3. 编码器 2 (交织)
    component_rsc_encoder(interleaved_msg, parity2);

    // 4. 码字复用 [u, c1, c2]
    // (注意：尾比特处理)
    int k = 0;
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        turbo_codeword[k++] = turbo_message_padded[i]; // u
        turbo_codeword[k++] = parity1[i];             // c1
        turbo_codeword[k++] = parity2[i];             // c2
    }
    // 补上尾比特的校验位
    for (int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        turbo_codeword[k++] = parity1[i];
        turbo_codeword[k++] = parity2[i];
    }
}

/**
 * @brief Turbo 码 BPSK 调制 (R=1/3)
 */
void turbo_modulation() {
    for (int i = 0; i < TURBO_codeword_length; i++) {
        turbo_tx_symbol[i][0] = -1.0 * (2.0 * turbo_codeword[i] - 1.0); // 0->1, 1->-1
        turbo_tx_symbol[i][1] = 0.0;
    }
}

/**
 * @brief Turbo 码 AWGN 信道 (R=1/3)
 */
void turbo_channel() {
    for (int i = 0; i < TURBO_codeword_length; i++) {
        for (int j = 0; j < 2; j++) {
            double u = (double)rand() / (double)RAND_MAX;
            if (u == 1.0) u = 0.999999;
            double r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));
            u = (double)rand() / (double)RAND_MAX;
            if (u == 1.0) u = 0.999999;
            double g = (double)r * cos(2.0 * pi * u);
            turbo_rx_symbol[i][j] = turbo_tx_symbol[i][j] + g;
        }
    }
}

/**
 * @brief Turbo 码迭代译码器 (顶层封装)
 */
void turbo_decoder_wrapper() {
    
    // 1. 计算信道 LLR
    // Lc = 2*y / sgm^2 = 4*y / N0
    // (因为 y 是 +/- 1, 接收到的是 r, Lc = 4*r/N0)
    double Lc_factor = 4.0 / N0;
    int k = 0;
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        Lc_sys[i]  = Lc_factor * turbo_rx_symbol[k++][0];
        Lc_par1[i] = Lc_factor * turbo_rx_symbol[k++][0];
        Lc_par2[i] = Lc_factor * turbo_rx_symbol[k++][0];
    }
    // 尾比特的 LLR
    for (int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        Lc_sys[i] = 0; // 尾比特的系统信息为 0 (高可信度)
        Lc_par1[i] = Lc_factor * turbo_rx_symbol[k++][0];
        Lc_par2[i] = Lc_factor * turbo_rx_symbol[k++][0];
    }

    // 2. 初始化先验信息
    for (int i = 0; i < TURBO_message_length; i++) {
        La_2_to_1[i] = 0.0;
    }

    // 3. 迭代译码
    for (int iter = 0; iter < TURBO_ITERATIONS; iter++) {
        
        // --- 译码器 1 ---
        // 输入: Lc(u), Lc(c1), La(u) from Dec2
        log_map_decoder(Lc_sys, Lc_par1, La_2_to_1, Le_1);
        
        // --- 交织 ---
        // Le_1 -> La_1_to_2
        interleave_llr(Le_1, La_1_to_2);
        
        // --- 译码器 2 ---
        // 输入: Lc(u'), Lc(c2), La(u') from Dec1
        // (需要交织 Lc_sys)
        double Lc_sys_interleaved[TURBO_message_length];
        interleave_llr(Lc_sys, Lc_sys_interleaved);
        
        log_map_decoder(Lc_sys_interleaved, Lc_par2, La_1_to_2, Le_2);

        // --- 解交织 ---
        // Le_2 -> La_2_to_1
        deinterleave_llr(Le_2, La_2_to_1);
    }

    // 4. 最终判决
    // L_APP(u) = Lc(u) + La(u) + Le(u)
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        L_APP[i] = Lc_sys[i] + La_2_to_1[i] + Le_1[i];
        if (L_APP[i] > 0.0) {
            turbo_de_message[i] = 0; // LLR > 0, 判为 0
        } else {
            turbo_de_message[i] = 1; // LLR < 0, 判为 1
        }
    }
}

/**
 * @brief 随机生成 Turbo 码消息
 */
void turbo_generate_message() {
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        turbo_message[i] = rand() % 2;
        turbo_message_padded[i] = turbo_message[i];
    }
    // 补零
    for (int i = TURBO_MESSAGE_BITS; i < TURBO_message_length; i++) {
        turbo_message_padded[i] = 0;
    }
}

/**
 * @brief 检查 Turbo 码误比特
 */
long int turbo_check_errors() {
    long int errors = 0;
    for (int i = 0; i < TURBO_MESSAGE_BITS; i++) {
        if (turbo_message[i] != turbo_de_message[i]) {
            errors++;
        }
    }
    return errors;
}


/**
 * @brief Turbo 码仿真主循环
 */
void run_turbo_simulation(float start_snr, float finish_snr, long int seq_num)
{
    long int bit_error, seq;
	double BER;
	double progress;

    // Turbo 码率 (R=1/3)
    float code_rate = (float)TURBO_MESSAGE_BITS / (float)(TURBO_MESSAGE_BITS * 3);

    printf("=================================================================\n");
    printf("               当前译码器: Turbo 码 (Log-MAP, %d 次迭代)\n", TURBO_ITERATIONS);
    printf("               码型: (7,5)_8 RSC, R=1/3\n");
    printf("=================================================================\n");
    
    // 生成交织器 (一次即可)
    turbo_generate_interleaver();

    // 按 SNR 循环
	for (float SNR = start_snr; SNR <= finish_snr; SNR++)
	{
		// 计算信道噪声参数
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2.0); 
		
		bit_error = 0; // 误比特数清零

        // 按帧数循环
		for (seq = 1; seq <= seq_num; seq++)
		{
            // 1. 生成消息
            turbo_generate_message();
            
            // 2. 编码
            turbo_encoder();

            // 3. 调制
            turbo_modulation();

            // 4. 信道
            turbo_channel();

            // 5. 迭代译码
            turbo_decoder_wrapper();

            // 6. 计误
            bit_error += turbo_check_errors();

			progress = (double)(seq * 100) / (double)seq_num; 
			BER = (double)bit_error / (double)(TURBO_MESSAGE_BITS * seq);

			printf("进度=%5.1f%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\r", progress, SNR, bit_error, BER);
		}

		BER = (double)bit_error / (double)(TURBO_MESSAGE_BITS * seq_num);
		printf("进度=100.0%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\n", SNR, bit_error, BER);
	}
}