#ifndef VITERBI_H_
#define VITERBI_H_

#include "config.h"

/**
 * @brief Viterbi 硬判决译码器
 * @param re_codeword   [in] 硬判决解调后的码字序列
 * @param de_message    [out] 译码后的消息序列
 * @param ms_length     [in] 消息序列长度 (含尾比特)
 */
void hardDecoder(int re_codeword[], int de_message[], int ms_length);

/**
 * @brief Viterbi 软判决译码器
 * @param re_codewordSoft [in] 经过信道的接收符号 (I/Q)
 * @param de_message      [out] 译码后的消息序列
 * @param ms_length       [in] 消息序列长度 (含尾比特)
 */
void softDecode(double re_codewordSoft[][softIn_st_num], int de_message[], int ms_length);

#endif // VITERBI_H_