/**
 * @file ccsds_turbo.h
 * @brief CCSDS 131.0-B-5 Turbo Code Implementation (NASA Standard)
 * 
 * 16-state RSC encoder with K=1784, Rate 1/3
 * Polynomials: G0=10011 (feedback), G1=11011 (feedforward)
 */

#ifndef CCSDS_TURBO_H_
#define CCSDS_TURBO_H_

#include "config.h"
#include "sim_runner.h"
#include <stdio.h>

// =================================================================
// --- CCSDS Turbo Code API ---
// =================================================================

/**
 * @brief 初始化 CCSDS Turbo 码 (生成交织器等)
 */
void ccsds_turbo_init(void);

/**
 * @brief CCSDS Turbo 编码器 (16-state RSC, R=1/3)
 */
void ccsds_turbo_encoder(void);

/**
 * @brief CCSDS Turbo BPSK 调制
 */
void ccsds_turbo_modulation(void);

/**
 * @brief CCSDS Turbo AWGN 信道
 */
void ccsds_turbo_channel(void);

/**
 * @brief CCSDS Turbo 迭代译码器 (Log-MAP)
 */
void ccsds_turbo_decoder(void);

/**
 * @brief 随机生成 CCSDS Turbo 码消息
 */
void ccsds_turbo_generate_message(void);

/**
 * @brief 检查 CCSDS Turbo 码误比特
 * @return long int 当前帧的误比特数
 */
long int ccsds_turbo_check_errors(void);

/**
 * @brief CCSDS Turbo 码仿真主循环
 * @param cfg 仿真配置
 * @param csv_fp CSV 文件指针
 */
void run_ccsds_turbo_simulation(SimConfig* cfg, FILE* csv_fp);

#endif // CCSDS_TURBO_H_

