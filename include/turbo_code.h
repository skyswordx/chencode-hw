#ifndef TURBO_CODE_H_
#define TURBO_CODE_H_

#include "config.h"

/**
 * @brief Turbo 码仿真主循环
 * @param start_snr 起始 SNR
 * @param finish_snr 结束 SNR
 * @param seq_num 仿真帧数
 */
void run_turbo_simulation(float start_snr, float finish_snr, long int seq_num);

/**
 * @brief 生成 Turbo 码的交织器
 */
void turbo_generate_interleaver();

/**
 * @brief Turbo 码编码器 (PCCC, R=1/3)
 */
void turbo_encoder();

/**
 * @brief Turbo 码 BPSK 调制 (R=1/3)
 */
void turbo_modulation();

/**
 * @brief Turbo 码 AWGN 信道 (R=1/3)
 */
void turbo_channel();

/**
 * @brief Turbo 码迭代译码器 (顶层封装)
 */
void turbo_decoder_wrapper();

/**
 * @brief 随机生成 Turbo 码消息
 */
void turbo_generate_message();

/**
 * @brief 检查 Turbo 码误比特
 * @return long int 当前帧的误比特数
 */
long int turbo_check_errors();


#endif // TURBO_CODE_H_