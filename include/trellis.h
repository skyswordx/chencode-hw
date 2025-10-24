#ifndef TRELLIS_H_
#define TRELLIS_H_

#include "config.h"

// --- 网格图数据结构定义 ---

/**
 * @brief 定义状态转移“边”的结构体
 */
typedef struct {
    int input;    // 引起该转移的输入比特 (0 或 1)
    int output;   // 该转移对应的输出码字索引 (c00, c01, c10, c11)
    int pt1;      // 转移的起始状态 (0, 1, 2, 3)
    int pt2;      // 转移的结束状态 (0, 1, 2, 3)
    int id;       // 边的唯一标识 ID (1 到 8)
} line;

/**
 * @brief 定义一个连接关系（点 + 边）
 */
typedef struct {
    int point;    // 连接的状态点
    int line;     // 连接的边
} connect;

/**
 * @brief 定义一个状态点的连接关系（Viterbi 用）
 * 描述到达某个状态有几条路径
 */
typedef struct {
    connect first;  // 第一条可能路径
    connect second; // 第二条可能路径
} graph_connect;


// --- 全局变量声明 (extern) ---
// extern 关键字表示“声明”一个变量，实体定义在 .c 文件中。

// (5,7)_8 卷积码的网格图常量表
extern const int state[4][2];          // 状态输出表
extern line stateTable[line_num];     // 状态转移表 (定义8条边)
extern graph_connect pathConn[st_num]; // 前向连接表 (用于 Viterbi ACS)
extern graph_connect pathConnB[st_num];// 后向连接表 (用于 BCJR Beta)

// 信道参数
extern double N0;  // 噪声功率谱密度
extern double sgm; // 噪声标准差 (sqrt(N0/2))

#endif // TRELLIS_H_