#ifndef CONFIG_H_
#define CONFIG_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// --- 译码器选择 ---
// 1 = 硬判决 Viterbi 译码器
// 2 = 软判决 Viterbi 译码器
// 3 = BCJR 译码器
#define DECODER_TYPE 2

// --- 码型参数 ---
#define MESSAGE_BITS 1024      // 核心信息比特数 (不含尾比特)
#define STATE_MEM 2            // 编码器寄存器个数 (记忆深度 m=2)
#define message_length (MESSAGE_BITS + STATE_MEM) // 补零后的消息总长度
#define codeword_length (message_length * 2)      // 码字长度 (码率为 1/2)
#define st_num 4               // 状态数 (2^STATE_MEM)
#define line_num 8             // 网格图中的转移边数 (2^(m+1))
#define softIn_st_num 2        // 软输入符号的维度 (I/Q 两路)

// --- 算法常量 ---
#define inf_int 1000000        // 用于硬判决 Viterbi 的“无穷大”路径度量
#define inf_double 1e18        // 用于软判决 Viterbi 的“无穷大”路径度量
#define pi 3.1415926           // 圆周率
#define p0 0.0                 // 概率初始值 (用于 BCJR)

#endif // CONFIG_H_