#ifndef CONFIG_H_
#define CONFIG_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// --- 译码器选择 ---
// 1 = 硬判决 Viterbi (卷积码)
// 2 = 软判决 Viterbi (卷积码)
// 3 = BCJR (卷积码)
// 4 = Turbo 码 (Log-MAP 迭代译码)
#define DECODER_TYPE 4

// --- 算法常量 ---
#define inf_int 1000000        // Viterbi 硬判决的“无穷大”
#define inf_double 1e18        // Viterbi 软判决的“无穷大”
#define pi 3.1415926           // 圆周率
#define p0 0.0                 // 概率初始值 (用于 BCJR)

// =================================================================
// --- 卷积码 (CC) 仿真参数 (用于 DECODER_TYPE 1, 2, 3) ---
// =================================================================
#define MESSAGE_BITS 1024       // 卷积码：核心信息比特数
#define STATE_MEM 2             // 卷积码：编码器寄存器个数 (m=2)
#define message_length (MESSAGE_BITS + STATE_MEM) // 卷积码：补零后的消息总长度
#define codeword_length (message_length * 2)      // 卷积码：码字长度 (R=1/2)
#define st_num 4                // 卷积码：状态数 (2^m)
#define line_num 8              // 卷积码：网格图转移边数
#define softIn_st_num 2         // 卷积码：软输入维度 (I/Q)


// =================================================================
// --- Turbo 码仿真参数 (用于 DECODER_TYPE 4) ---
// =================================================================
#define TURBO_MESSAGE_BITS 1024    // Turbo 码：核心信息比特数 (K)
#define TURBO_STATE_MEM 2          // Turbo 码：分量编码器寄存器数 (m=2)
#define TURBO_message_length (TURBO_MESSAGE_BITS + TURBO_STATE_MEM) // Turbo 码：补零后长度
#define TURBO_codeword_length (TURBO_MESSAGE_BITS * 3 + TURBO_STATE_MEM * 2) // Turbo 码：码字长度 R=1/3, [u, c1, c2], 加上尾比特
#define TURBO_st_num 4             // Turbo 码：分量编码器状态数
#define TURBO_line_num 8           // Turbo 码：分量编码器转移边数
#define TURBO_ITERATIONS 8         // Turbo 码：译码器迭代次数

#endif // CONFIG_H_