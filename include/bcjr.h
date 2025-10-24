#ifndef BCJR_H_
#define BCJR_H_

#include "config.h"

/**
 * @brief BCJR (Bahl-Cocke-Jelinek-Raviv) 译码器
 * @param rx_sym      [in] 经过信道的接收符号 (I/Q)
 * @param de_message  [out] 译码后的消息序列
 * @param m_length    [in] 消息序列长度 (含尾比特)
 * @param c_length    [in] 码字序列长度
 */
void BCJR(double rx_sym[][softIn_st_num], int de_message[], int m_length, int c_length);

#endif // BCJR_H_