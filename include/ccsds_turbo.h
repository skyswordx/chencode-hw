/**
 * @file ccsds_turbo.h
 * @brief CCSDS 131.0-B-5 Turbo Code Implementation (NASA Standard)
 * 
 * 16-state RSC encoder with dynamic K and CCSDS nominal rates
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

/** Return the configured nominal code rate (for example, 1/2). */
double ccsds_get_nominal_rate(void);

/** Return the exact transmitted codeword length n=(K+4)/r. */
int ccsds_get_codeword_length(void);

/** Return the true finite-block rate K/n used for Eb/N0 normalization. */
double ccsds_get_effective_rate(void);

/** Run deterministic RSC, puncturing, and noiseless round-trip checks. */
int ccsds_turbo_self_test(void);

/**
 * @brief CCSDS Turbo 编码器 (16-state RSC, R=1/2 or R=1/3)
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

