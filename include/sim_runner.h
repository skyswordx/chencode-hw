#ifndef SIM_RUNNER_H_
#define SIM_RUNNER_H_

#include "config.h"
#include <stdio.h>

// =================================================================
// --- 仿真配置结构体 ---
// =================================================================
typedef struct {
    float start_snr;      // 起始 SNR (dB)
    float end_snr;        // 结束 SNR (dB)
    float snr_step;       // SNR 步长 (dB)
    long  num_frames;     // 每个 SNR 点的仿真帧数
    char  csv_filename[256]; // CSV 输出文件名
} SimConfig;

// =================================================================
// --- 仿真运行函数声明 ---
// =================================================================

/**
 * @brief 运行仿真的统一入口
 * @param decoder 译码器类型
 * @param cfg 仿真配置
 */
void run_simulation(DecoderType decoder, SimConfig* cfg);

/**
 * @brief Uncoded BPSK 仿真 (理论对照)
 */
void run_uncoded_simulation(SimConfig* cfg, FILE* csv_fp);

/**
 * @brief 卷积码仿真 (硬/软 Viterbi, BCJR)
 */
void run_cc_simulation_v2(DecoderType decoder, SimConfig* cfg, FILE* csv_fp);

/**
 * @brief Turbo 码仿真
 */
void run_turbo_simulation_v2(SimConfig* cfg, FILE* csv_fp);

/**
 * @brief 打印仿真配置摘要
 */
void print_sim_header(DecoderType decoder, SimConfig* cfg);

/**
 * @brief 打印仿真结果表格行
 */
void print_result_row(double snr, long bit_errors, long total_bits, double ber);

/**
 * @brief 打印仿真结果表头
 */
void print_result_table_header(void);

/**
 * @brief 打印仿真结果表尾
 */
void print_result_table_footer(void);

#endif // SIM_RUNNER_H_
