#include "sim_runner.h"
#include "csv_export.h"
#include "config.h"
#include "convolutional_code.h"
#include "trellis.h"
#include "viterbi.h"
#include "bcjr.h"
#include "turbo_code.h"
#include "ccsds_turbo.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

// 全局译码器类型变量定义
DecoderType g_decoder_type = DECODER_SOFT_VITERBI;

// 全局 Turbo 类型变量定义
TurboType g_turbo_type = TURBO_TYPE_75_OCT;  // 默认为原有 (7,5)_8 实现

// 外部变量声明
extern double N0;
extern double sgm;

// =================================================================
// --- 终端输出格式化函数 ---
// =================================================================

void print_sim_header(DecoderType decoder, SimConfig* cfg) {
    const char* decoder_name = DECODER_NAMES[decoder];
    
    int block_size;
    double code_rate;
    const char* turbo_variant = "";
    
    if (decoder == DECODER_TURBO) {
        if (g_turbo_type == TURBO_TYPE_CCSDS) {
            block_size = CCSDS_K;
            turbo_variant = " (CCSDS)";
        } else {
            block_size = TURBO_MESSAGE_BITS;
            turbo_variant = " ((7,5)_8)";
        }
        code_rate = 1.0 / 3.0;
    } else if (decoder == DECODER_UNCODED) {
        block_size = CC_MESSAGE_BITS;
        code_rate = 1.0;
    } else {
        block_size = CC_MESSAGE_BITS;
        code_rate = 0.5;
    }
    
    printf("\n");
    printf("+============================================================+\n");
    printf("|          Channel Coding Simulation System v2.0             |\n");
    printf("+============================================================+\n");
    printf("| Decoder:     %-44s |\n", decoder_name);
    printf("| Block Size:  %-4d bits                                     |\n", block_size);
    printf("| Code Rate:   %-6.4f                                       |\n", code_rate);
    printf("| SNR Range:   %.1f ~ %.1f dB (step=%.1f)                       |\n", 
           cfg->start_snr, cfg->end_snr, cfg->snr_step);
    printf("| Frames/SNR:  %-8ld                                      |\n", cfg->num_frames);
    printf("| Output CSV:  %-44s |\n", cfg->csv_filename);
    printf("+============================================================+\n\n");
}

void print_result_table_header(void) {
    printf("+----------+---------------+---------------+----------------+\n");
    printf("|  Eb/N0   |  Bit Errors   |  Total Bits   |      BER       |\n");
    printf("+----------+---------------+---------------+----------------+\n");
}

void print_result_row(double snr, long bit_errors, long total_bits, double ber) {
    printf("|  %5.1f   |  %11ld  |  %11ld  |  %12.4e  |\n", 
           snr, bit_errors, total_bits, ber);
}

void print_result_table_footer(void) {
    printf("+----------+---------------+---------------+----------------+\n");
}

// =================================================================
// --- Uncoded BPSK 仿真 ---
// =================================================================

void run_uncoded_simulation(SimConfig* cfg, FILE* csv_fp) {
    double code_rate = 1.0; // Uncoded, R=1
    int block_size = CC_MESSAGE_BITS;
    
    long int bit_error, frame_error, seq;
    double BER, FER;
    
    // Print header with FER column
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    printf("|  Eb/N0   |  Bit Errors   |  Total Bits   |      BER       | Frame Errors  |      FER       |\n");
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    
    for (float SNR = cfg->start_snr; SNR <= cfg->end_snr + 0.001; SNR += cfg->snr_step) {
        // 计算噪声参数
        N0 = 1.0 / pow(10.0, SNR / 10.0);
        sgm = sqrt(N0 / 2.0);
        
        bit_error = 0;
        frame_error = 0;
        
        for (seq = 1; seq <= cfg->num_frames; seq++) {
            long errors_this_frame = 0;
            
            // 生成随机比特并直接 BPSK 调制/解调
            for (int i = 0; i < block_size; i++) {
                int bit = rand() % 2;
                double tx = (bit == 0) ? 1.0 : -1.0; // BPSK: 0->+1, 1->-1
                
                // Box-Muller 产生高斯噪声
                double u1 = (double)rand() / RAND_MAX;
                if (u1 == 0.0) u1 = 0.0001;
                double u2 = (double)rand() / RAND_MAX;
                double noise = sgm * sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
                
                double rx = tx + noise;
                
                // 硬判决
                int decoded = (rx > 0) ? 0 : 1;
                
                if (bit != decoded) {
                    bit_error++;
                    errors_this_frame++;
                }
            }
            
            if (errors_this_frame > 0) frame_error++;
        }
        
        long total_bits = (long)block_size * cfg->num_frames;
        BER = (double)bit_error / (double)total_bits;
        FER = (double)frame_error / (double)cfg->num_frames;
        
        printf("|  %5.1f   |  %11ld  |  %11ld  |  %12.4e  |  %11ld  |  %12.4e  |\n", 
               SNR, bit_error, total_bits, BER, frame_error, FER);
        csv_append_row_with_fer(csv_fp, SNR, bit_error, total_bits, BER, frame_error, cfg->num_frames, FER);
    }
    
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
}

// =================================================================
// --- 卷积码仿真 (硬/软 Viterbi, BCJR) ---
// =================================================================

void run_cc_simulation_v2(DecoderType decoder, SimConfig* cfg, FILE* csv_fp) {
    double code_rate = (double)CC_MESSAGE_BITS / (double)CC_codeword_length;
    
    long int bit_error, frame_error, seq;
    double BER, FER;
    int i;
    
    // Print header with FER column
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    printf("|  Eb/N0   |  Bit Errors   |  Total Bits   |      BER       | Frame Errors  |      FER       |\n");
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    
    for (float SNR = cfg->start_snr; SNR <= cfg->end_snr + 0.001; SNR += cfg->snr_step) {
        // 计算噪声参数
        N0 = (1.0 / code_rate) / pow(10.0, SNR / 10.0);
        sgm = sqrt(N0 / 2.0);
        
        bit_error = 0;
        frame_error = 0;
        
        for (seq = 1; seq <= cfg->num_frames; seq++) {
            // 1. 生成随机消息
            for (i = 0; i < CC_MESSAGE_BITS; i++) {
                message[i] = rand() % 2;
            }
            // 2. 补零 (尾比特)
            for (i = CC_MESSAGE_BITS; i < CC_message_length; i++) {
                message[i] = 0;
            }
            
            // 3. 卷积编码
            encoder();
            
            // 4. BPSK 调制
            modulation();
            
            // 5. AWGN 信道
            channel();
            
            // 6. 解调/译码
            switch (decoder) {
                case DECODER_HARD_VITERBI:
                    demodulation();
                    hardDecoder(re_codeword, de_message, CC_message_length);
                    break;
                case DECODER_SOFT_VITERBI:
                    softDecode(rx_symbol, de_message, CC_message_length);
                    break;
                case DECODER_BCJR:
                    BCJR(rx_symbol, de_message, CC_message_length, CC_codeword_length);
                    break;
                default:
                    break;
            }
            
            // 7. 计误 (count errors per frame)
            long errors_this_frame = 0;
            for (i = 0; i < CC_MESSAGE_BITS; i++) {
                if (message[i] != de_message[i]) {
                    bit_error++;
                    errors_this_frame++;
                }
            }
            if (errors_this_frame > 0) frame_error++;
        }
        
        long total_bits = (long)CC_MESSAGE_BITS * cfg->num_frames;
        BER = (double)bit_error / (double)total_bits;
        FER = (double)frame_error / (double)cfg->num_frames;
        
        printf("|  %5.1f   |  %11ld  |  %11ld  |  %12.4e  |  %11ld  |  %12.4e  |\n", 
               SNR, bit_error, total_bits, BER, frame_error, FER);
        csv_append_row_with_fer(csv_fp, SNR, bit_error, total_bits, BER, frame_error, cfg->num_frames, FER);
    }
    
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
}

// =================================================================
// --- Turbo 码仿真 ---
// =================================================================

// 外部 Turbo 码函数声明
extern void turbo_generate_interleaver(void);
extern void turbo_generate_message(void);
extern void turbo_encoder(void);
extern void turbo_modulation(void);
extern void turbo_channel(void);
extern void turbo_decoder_wrapper(void);
extern long int turbo_check_errors(void);

void run_turbo_simulation_v2(SimConfig* cfg, FILE* csv_fp) {
    double code_rate = (double)TURBO_MESSAGE_BITS / (double)(TURBO_MESSAGE_BITS * 3);
    
    long int bit_error, frame_error, seq;
    double BER, FER;
    
    // Early termination threshold - stop when this many frame errors are collected
    const long MIN_FRAME_ERRORS = 100;
    
    // 生成交织器 (仅一次)
    turbo_generate_interleaver();
    
    // Print header with FER column
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    printf("|  Eb/N0   |  Bit Errors   |  Total Bits   |      BER       | Frame Errors  |      FER       |\n");
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
    
    for (float SNR = cfg->start_snr; SNR <= cfg->end_snr + 0.001; SNR += cfg->snr_step) {
        // 计算噪声参数
        N0 = (1.0 / code_rate) / pow(10.0, SNR / 10.0);
        sgm = sqrt(N0 / 2.0);
        
        bit_error = 0;
        frame_error = 0;
        long actual_frames = 0;
        
        for (seq = 1; seq <= cfg->num_frames; seq++) {
            turbo_generate_message();
            turbo_encoder();
            turbo_modulation();
            turbo_channel();
            turbo_decoder_wrapper();
            
            long errors_this_frame = turbo_check_errors();
            bit_error += errors_this_frame;
            if (errors_this_frame > 0) {
                frame_error++;  // Count frames with at least 1 error
            }
            actual_frames = seq;
            
            // Early termination: stop if we have enough frame errors
            if (frame_error >= MIN_FRAME_ERRORS && seq >= 1000) {
                break;
            }
            
            // Progress indicator every 1000 frames
            if (seq % 1000 == 0) {
                printf("\r  [SNR=%.1fdB] Frame %ld/%ld (%.1f%%) - BitErr: %ld, FrameErr: %ld   ", 
                       SNR, seq, cfg->num_frames, 
                       100.0 * seq / cfg->num_frames, bit_error, frame_error);
                fflush(stdout);
            }
        }
        printf("\r                                                                        \r");
        
        // Calculate BER and FER
        long total_bits = (long)TURBO_MESSAGE_BITS * actual_frames;
        BER = (double)bit_error / (double)total_bits;
        FER = (double)frame_error / (double)actual_frames;
        
        // Print result row with FER
        printf("|  %5.1f   |  %11ld  |  %11ld  |  %12.4e  |  %11ld  |  %12.4e  |\n", 
               SNR, bit_error, total_bits, BER, frame_error, FER);
        
        // Write to CSV with FER
        csv_append_row_with_fer(csv_fp, SNR, bit_error, total_bits, BER, frame_error, actual_frames, FER);
    }
    
    printf("+----------+---------------+---------------+----------------+---------------+----------------+\n");
}

// =================================================================
// --- 统一仿真入口 ---
// =================================================================

void run_simulation(DecoderType decoder, SimConfig* cfg) {
    // 确定码率和帧长
    double code_rate;
    int block_size;
    
    switch (decoder) {
        case DECODER_UNCODED:
            code_rate = 1.0;
            block_size = CC_MESSAGE_BITS;
            break;
        case DECODER_TURBO:
            code_rate = 1.0 / 3.0;
            block_size = (g_turbo_type == TURBO_TYPE_CCSDS) ? CCSDS_K : TURBO_MESSAGE_BITS;
            break;
        default:
            code_rate = 0.5;
            block_size = CC_MESSAGE_BITS;
            break;
    }
    
    // 生成 CSV 文件名 (如未指定)
    if (strlen(cfg->csv_filename) == 0) {
        csv_generate_filename(cfg->csv_filename, sizeof(cfg->csv_filename), decoder);
    }
    
    // 打印配置摘要
    print_sim_header(decoder, cfg);
    
    // 打开 CSV 文件
    FILE* csv_fp = csv_init(cfg->csv_filename, DECODER_NAMES[decoder], code_rate, block_size);
    if (!csv_fp) {
        fprintf(stderr, "[错误] 无法创建 CSV 文件，继续仿真但不保存数据...\n");
    }
    
    // 根据译码器类型运行对应仿真
    switch (decoder) {
        case DECODER_UNCODED:
            run_uncoded_simulation(cfg, csv_fp);
            break;
        case DECODER_HARD_VITERBI:
        case DECODER_SOFT_VITERBI:
        case DECODER_BCJR:
            run_cc_simulation_v2(decoder, cfg, csv_fp);
            break;
        case DECODER_TURBO:
            if (g_turbo_type == TURBO_TYPE_CCSDS) {
                run_ccsds_turbo_simulation(cfg, csv_fp);  // CCSDS 16-state
            } else {
                run_turbo_simulation_v2(cfg, csv_fp);     // Original (7,5)_8
            }
            break;
        default:
            fprintf(stderr, "[错误] 未知译码器类型: %d\n", decoder);
            break;
    }
    
    // 关闭 CSV 文件
    csv_close(csv_fp);
    
    printf("\n[完成] 结果已保存到: %s\n", cfg->csv_filename);
}
