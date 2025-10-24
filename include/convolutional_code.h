#ifndef CONVOLUTIONAL_CODE_H_
#define CONVOLUTIONAL_CODE_H_

#include "config.h"

// --- 全局数组声明 ---
// 存放消息、码字和符号
extern int message[message_length];     // 消息序列 (含尾比特)
extern int codeword[codeword_length];   // 编码后的码字序列
extern int re_codeword[codeword_length];  // 硬判决解调后的码字序列
extern int de_message[message_length];  // 译码后的消息序列

extern double tx_symbol[codeword_length][2]; // BPSK 调制后的发送符号 (I/Q)
extern double rx_symbol[codeword_length][2]; // 经过 AWGN 信道后的接收符号 (I/Q)

// --- 函数原型声明 ---

/**
 * @brief (7,5)_8 卷积编码器
 */
void encoder();

/**
 * @brief BPSK 调制
 * 0 映射到 (1, 0), 1 映射到 (-1, 0)
 */
void modulation();

/**
 * @brief AWGN (加性高斯白噪声) 信道模型
 * 使用 Box-Muller 方法产生高斯噪声
 */
void channel();

/**
 * @brief BPSK 硬判决解调
 * 仅用于硬判决 Viterbi 译码器
 */
void demodulation();

#endif // CONVOLUTIONAL_CODE_H_