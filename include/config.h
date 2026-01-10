#ifndef CONFIG_H_
#define CONFIG_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// =================================================================
// --- 译码器类型枚举 (运行时选择) ---
// =================================================================
typedef enum {
    DECODER_UNCODED = 0,      // Uncoded BPSK (理论对照)
    DECODER_HARD_VITERBI = 1, // 硬判决 Viterbi
    DECODER_SOFT_VITERBI = 2, // 软判决 Viterbi
    DECODER_BCJR = 3,         // BCJR (MAP)
    DECODER_TURBO = 4         // Turbo 码 (Log-MAP)
} DecoderType;

// 全局译码器类型变量 (在 main.c 中定义)
extern DecoderType g_decoder_type;

// =================================================================
// --- Turbo 码变体选择 ---
// =================================================================
typedef enum {
    TURBO_TYPE_75_OCT = 0,    // (7,5)_8 4-state RSC (原有实现, K=1024)
    TURBO_TYPE_CCSDS = 1      // CCSDS 16-state RSC (NASA标准, K=1784)
} TurboType;

// 全局 Turbo 类型变量 (在 sim_runner.c 中定义)
extern TurboType g_turbo_type;

// =================================================================
// --- 算法常量 ---
// =================================================================
#define INF_INT 1000000        // Viterbi 硬判决的"无穷大"
#define INF_DOUBLE 1e18        // Viterbi 软判决的"无穷大"
#define PI 3.14159265358979    // 圆周率
#define P0 0.0                 // 概率初始值 (用于 BCJR)

// 兼容旧的小写宏名 (避免大量修改)
#define inf_int INF_INT
#define inf_double INF_DOUBLE
#define pi PI
#define p0 P0

// =================================================================
// --- 通用网格图参数 (CC 和 Turbo 共用) ---
// =================================================================
#define ST_NUM 4               // 状态数 (2^m)
#define LINE_NUM 8             // 网格图转移边数
#define SOFTIN_ST_NUM 2        // 软输入维度 (I/Q)
#define STATE_MEM 2            // 编码器寄存器个数 (记忆深度 m=2)

// 兼容旧的小写宏名
#define st_num ST_NUM
#define line_num LINE_NUM
#define softIn_st_num SOFTIN_ST_NUM

// =================================================================
// --- 卷积码 (CC) 仿真参数 ---
// =================================================================
#define CC_MESSAGE_BITS 1024       // 卷积码：核心信息比特数
#define CC_message_length (CC_MESSAGE_BITS + STATE_MEM) // 补零后的消息总长度
#define CC_codeword_length (CC_message_length * 2)      // 码字长度 (R=1/2)

// =================================================================
// --- Turbo 码仿真参数 (原有 (7,5)_8 实现) ---
// =================================================================
#define TURBO_MESSAGE_BITS 1024    // Turbo 码：核心信息比特数 (K)
#define TURBO_message_length (TURBO_MESSAGE_BITS + STATE_MEM) // 补零后长度
#define TURBO_codeword_length (TURBO_MESSAGE_BITS * 3 + STATE_MEM * 4) // R≈1/3
#define TURBO_ITERATIONS 8         // Turbo 码：译码器迭代次数

// =================================================================
// --- CCSDS Turbo 码参数 (NASA 标准 16-state) ---
// =================================================================
#define CCSDS_ST_NUM 16            // CCSDS: 16 状态 (2^4)
#define CCSDS_LINE_NUM 32          // CCSDS: 32 条转移边
#define CCSDS_STATE_MEM 4          // CCSDS: 记忆深度 m=4
#define CCSDS_ITERATIONS 10        // CCSDS: 译码器迭代次数 (NASA基准为10)

// 动态 K 值支持 (运行时选择)
#define CCSDS_K_MAX 8920           // 最大支持的 K 值 (用于静态数组分配)
#define CCSDS_K_1784  1784         // K=1784 (k1=8, k2=223)
#define CCSDS_K_3568  3568         // K=3568 (k1=8, k2=446)
#define CCSDS_K_7136  7136         // K=7136 (k1=8, k2=892)
#define CCSDS_K_8920  8920         // K=8920 (k1=8, k2=1115)

// 运行时 K 值 (在 ccsds_turbo.c 中定义)
extern int g_ccsds_k;

// 静态数组使用最大尺寸分配
#define CCSDS_message_length_max (CCSDS_K_MAX + CCSDS_STATE_MEM)
#define CCSDS_codeword_length_max (CCSDS_K_MAX * 3 + CCSDS_STATE_MEM * 4)

// =================================================================
// --- 译码器名称字符串 ---
// =================================================================
static const char* DECODER_NAMES[] = {
    "Uncoded BPSK",
    "Hard Viterbi (CC R=1/2)",
    "Soft Viterbi (CC R=1/2)",
    "BCJR / MAP (CC R=1/2)",
    "Turbo (PCCC R=1/3)"
};

#endif // CONFIG_H_