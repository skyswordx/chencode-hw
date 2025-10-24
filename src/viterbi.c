#include "viterbi.h"
#include "config.h"
#include "trellis.h" // Viterbi 算法需要网格图定义

// --- Viterbi 译码器所需的全局数组 ---
// 存放分支度量和路径度量
int branchTable[message_length][line_num];    // 分支度量表 (硬判决：汉明距离)
int pathTable[message_length + 1][st_num];  // 路径度量表 (累加度量)
int trellisTable[message_length][st_num];   // 幸存路径表 (存储 trellis 转移)
int minPath[message_length];                // 最终译码路径

double branchTableSoft[message_length][line_num]; // 分支度量表 (软判决：欧氏距离)
double pathTableSoft[message_length + 1][st_num];// 路径度量表 (软判决：累加度量)


/**
 * @brief Viterbi 硬判决译码器
 */
void hardDecoder(int re_codeword[], int de_message[], int ms_length)
{
    // 1. 计算分支度量 (汉明距离)
    int hammingDistance;
    for (int t = 0; t < ms_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            hammingDistance = 0;
            // 过滤掉不可能的转移 (用于处理 trellis 的开始和结束)
            // 这些判断是基于补零(Trellis 终止)的假设
            if (t == 0 && index != 0 && index != 1) { branchTable[t][index] = inf_int; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { branchTable[t][index] = inf_int; continue; }
            if (t == ms_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { branchTable[t][index] = inf_int; continue; }
            if (t == ms_length - 1 && index != 0 && index != 2) { branchTable[t][index] = inf_int; continue; }

            // 获取当前转移边对应的期望输出
            int nowCode1 = state[stateTable[index].output][0];
            int nowCode2 = state[stateTable[index].output][1];

            // 计算汉明距离
            if (nowCode1 != re_codeword[2 * t]) hammingDistance++;
            if (nowCode2 != re_codeword[2 * t + 1]) hammingDistance++;
            
            branchTable[t][index] = hammingDistance;
        }
    }

    // 2. 加-比-选 (Add-Compare-Select, ACS)
    int path1d, path2d, minpath;
    // 初始化 t=0 时刻的路径度量
    pathTable[0][0] = 0; // 起始状态 S0 度量为 0
    for (int fst = 1; fst < st_num; fst++) { pathTable[0][fst] = inf_int; } // 其他状态为无穷

    // 迭代计算 t = 1 到 t = ms_length
    for (int pt = 1; pt < ms_length + 1; pt++) // pt 是时刻
    {
        for (int st = 0; st < st_num; st++) // st 是当前状态
        {
            // 查找 pathConn[st] 得到达当前状态 st 的两条路径
            // 路径1：来自 pathConn[st].first.point 状态，经过 pathConn[st].first.line 边
            path1d = pathTable[pt - 1][pathConn[st].first.point] + branchTable[pt - 1][pathConn[st].first.line];
            // 路径2：
            path2d = pathTable[pt - 1][pathConn[st].second.point] + branchTable[pt - 1][pathConn[st].second.line];

            // 再次过滤不可能的转移
            if (pt == 1 && st != 0 && st != 2) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }
            if (pt == ms_length - 1 && st != 0 && st != 1) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }
            if (pt == ms_length && st != 0) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }

            // 比较 (Compare) 并选择 (Select)
            if (path1d < path2d) {
                minpath = path1d; // 存活路径的累加度量
                trellisTable[pt - 1][st] = stateTable[pathConn[st].first.line].id; // 记录存活的边 ID
            }
            else {
                minpath = path2d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].second.line].id;
            }
            pathTable[pt][st] = minpath; // 存储 (Add) 累加度量
        }
    }

    // 3. 回溯 (Traceback)
    // 补零终止后，最终状态必然是 S0 (st=0)
    int nowLine = 0, nowPoint = 0, prevPoint = 0;
    for (int tt = ms_length; tt > 0; tt--) // 从 t = ms_length 开始回溯
    {
        nowLine = trellisTable[tt - 1][nowPoint] - 1; // 获取存活的边 ID (数组索引-1)
        prevPoint = stateTable[nowLine].pt1;     // 找到这条边的起始状态
        nowPoint = prevPoint;                    // 更新回溯的状态点
        minPath[tt - 1] = nowLine + 1;           // 存储最优路径的边 ID
    }

    // 4. 根据最优路径输出译码消息
    for (int t = 0; t < ms_length; t++)
    {
        de_message[t] = stateTable[minPath[t] - 1].input; // 查找边 ID 对应的输入比特
    }
}


/**
 * @brief Viterbi 软判决译码器
 */
void softDecode(double re_codewordSoft[][softIn_st_num], int de_message[], int ms_length)
{
    // 1. 计算分支度量 (欧氏距离的平方)
    double euDistance;
    for (int t = 0; t < ms_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            euDistance = 0;
            // 过滤不可能的转移
            if (t == 0 && index != 0 && index != 1) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == ms_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == ms_length - 1 && index != 0 && index != 2) { branchTableSoft[t][index] = inf_double; continue; }

            int nowCode1 = state[stateTable[index].output][0];
            int nowCode2 = state[stateTable[index].output][1];

            // 0 映射到 +1, 1 映射到 -1
            double ideal_x1 = (nowCode1 == 0) ? 1.0 : -1.0;
            double ideal_x2 = (nowCode2 == 0) ? 1.0 : -1.0;

            // 计算欧氏距离的平方 (Squared Euclidean Distance)
            // (省略 sqrt 计算更快，且不影响大小比较)
            euDistance += pow(re_codewordSoft[2 * t][0] - ideal_x1, 2);
            // rx_symbol 的 [1] (Q路) 理论上是0，但也可能包含噪声
            euDistance += pow(re_codewordSoft[2 * t][1] - 0.0, 2); 
            euDistance += pow(re_codewordSoft[2 * t + 1][0] - ideal_x2, 2);
            euDistance += pow(re_codewordSoft[2 * t + 1][1] - 0.0, 2);

            branchTableSoft[t][index] = euDistance;
        }
    }

    // 2. 加-比-选 (ACS)
    // (与硬判决的逻辑完全相同，只是数据类型变为 double)
    double path1d, path2d, minpath;
    pathTableSoft[0][0] = 0;
    for (int fst = 1; fst < st_num; fst++) { pathTableSoft[0][fst] = inf_double; }

    for (int pt = 1; pt < ms_length + 1; pt++)
    {
        for (int st = 0; st < st_num; st++)
        {
            path1d = pathTableSoft[pt - 1][pathConn[st].first.point] + branchTableSoft[pt - 1][pathConn[st].first.line];
            path2d = pathTableSoft[pt - 1][pathConn[st].second.point] + branchTableSoft[pt - 1][pathConn[st].second.line];

            if (pt == 1 && st != 0 && st != 2) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            if (pt == ms_length - 1 && st != 0 && st != 1) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            if (pt == ms_length && st != 0) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            
            if (path1d < path2d) {
                minpath = path1d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].first.line].id;
            }
            else {
                minpath = path2d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].second.line].id;
            }
            pathTableSoft[pt][st] = minpath;
        }
    }

    // 3. 回溯 (Traceback)
    // (与硬判决完全相同)
    int nowLine = 0, nowPoint = 0, prevPoint = 0;
    for (int tt = ms_length; tt > 0; tt--)
    {
        nowLine = trellisTable[tt - 1][nowPoint] - 1; 
        prevPoint = stateTable[nowLine].pt1;     
        nowPoint = prevPoint;                         
        minPath[tt - 1] = nowLine + 1;                 
    }

    // 4. 输出译码消息
    // (与硬判决完全相同)
    for (int t = 0; t < ms_length; t++)
    {
        de_message[t] = stateTable[minPath[t] - 1].input;
    }
}