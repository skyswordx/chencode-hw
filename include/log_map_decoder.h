#ifndef LOG_MAP_DECODER_H_
#define LOG_MAP_DECODER_H_

#include "config.h"

/**
 * @brief Jacobian 对数 (max* 运算)
 * log(exp(a) + exp(b)) = max(a, b) + log(1 + exp(-abs(a-b)))
 * @param a 对数值 a
 * @param b 对数值 b
 * @return log(exp(a) + exp(b))
 */
double jac_log(double a, double b);

/**
 * @brief Log-MAP 核心译码器 (对数域 BCJR)
 * @param Lc_sys    [in]  系统比特的信道 LLR
 * @param Lc_par    [in]  校验比特的信道 LLR
 * @param La_in     [in]  先验信息 LLR (来自另一个译码器)
 * @param Le_out    [out] 外在信息 LLR (发送给另一个译码器)
 * @param trellis   [in]  所使用的网格图定义 (用于 RSC)
 */
void log_map_decoder(double* Lc_sys, double* Lc_par, double* La_in, double* Le_out);

#endif // LOG_MAP_DECODER_H_