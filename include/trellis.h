#ifndef TRELLIS_H_
#define TRELLIS_H_

#include "config.h"

// --- 网格图数据结构定义 ---

typedef struct {
    int input;
    int output;
    int pt1;
    int pt2;
    int id;
} line;

typedef struct {
    int point;
    int line;
} connect;

typedef struct {
    connect first;
    connect second;
} graph_connect;


// --- 全局变量声明 (extern) ---
// (5,7)_8 卷积码的网格图常量表
extern const int state[4][2];          // 状态输出表
extern line stateTable[line_num];     // 状态转移表 (使用通用定义)
extern graph_connect pathConn[st_num]; // 前向连接表 (使用通用定义)
extern graph_connect pathConnB[st_num];// 后向连接表 (使用通用定义)

// 信道参数
extern double N0;  // 噪声功率谱密度
extern double sgm; // 噪声标准差 (sqrt(N0/2))

#endif // TRELLIS_H_