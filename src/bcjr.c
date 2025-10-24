#include "bcjr.h"
#include "config.h"
#include "trellis.h" // BCJR 算法需要网格图定义

// --- BCJR 译码器所需的全局数组 ---
// 存放前向、后向和分支概率
double pCh[codeword_length][2];     // 信道观察概率 P(y|x) [x=0 或 x=1]
double pLine[message_length][line_num]; // 分支转移概率 Gamma
double pA[message_length + 1][st_num];  // 前向概率 Alpha
double pB[message_length + 1][st_num];  // 后向概率 Beta


/**
 * @brief 计算欧氏距离的平方 (静态辅助函数)
 * @param rx_x 接收符号的 I 分量
 * @param rx_y 接收符号的 Q 分量
 * @param sym  期望发送的符号 (0 或 1)
 * @return 欧氏距离的平方
 */
static double euDist(double rx_x, double rx_y, int sym)
{
    // 0 映射到 +1, 1 映射到 -1
    double ideal_x = (sym == 0) ? 1.0 : -1.0;
    // 返回欧氏距离的平方
    return pow(rx_x - ideal_x, 2) + pow(rx_y - 0.0, 2);
}

/**
 * @brief 计算信道观察概率 P(y|x) (静态辅助函数)
 * P(y | x) = (1 / sqrt(pi*N0)) * exp(-d^2(y,x) / N0)
 * @param eu_sq 欧氏距离的平方 d^2(y,x)
 * @param n0    噪声功率谱密度 N0
 * @return P(y|x)
 */
static double chObs(double eu_sq, double n0)
{
    // 注意: N0 = 2 * sgm^2
    return exp(-eu_sq / n0) / (sqrt(pi * n0));
}


/**
 * @brief BCJR 译码器主函数
 */
void BCJR(double rx_sym[][softIn_st_num], int de_message[], int m_length, int c_length)
{
	// 1. 计算信道观察概率 P(y|x)
    // pCh[t][0] = P(y_t | x=0), pCh[t][1] = P(y_t | x=1)
    double sum;
    for (int t = 0; t < c_length; t++)
    {
        sum = 0;
        // 计算 P(y_t | 0)
        pCh[t][0] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 0), N0);
        // 计算 P(y_t | 1)
        pCh[t][1] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 1), N0);
        
        // 归一化处理，防止概率连乘导致浮点数下溢
        sum = pCh[t][0] + pCh[t][1];
        if (sum > 1e-300) { // 避免除零
            pCh[t][0] = pCh[t][0] / sum;
            pCh[t][1] = pCh[t][1] / sum;
        } else {
            pCh[t][0] = 0.5;
            pCh[t][1] = 0.5;
        }
    }
    
    // 2. 计算分支转移概率 (Gamma)
    // Gamma_t(s', s) = P(u_t) * P(y_t | s', s)
    // P(u_t) = 0.5 (假设先验概率 P(0)=P(1)=0.5)
    // P(y_t | s', s) = P(y_t1 | c1) * P(y_t2 | c2)
    for (int t = 0; t < m_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            // 过滤不可能的转移
            if (t == 0 && index != 0 && index != 1) { pLine[t][index] = p0; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { pLine[t][index] = p0; continue; }
            if (t == m_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { pLine[t][index] = p0; continue; }
            if (t == m_length - 1 && index != 0 && index != 2) { pLine[t][index] = p0; continue; }
            
            int out_idx = stateTable[index].output; // 获取输出索引
            int c1 = state[out_idx][0]; // 获取输出的第1个码字
            int c2 = state[out_idx][1]; // 获取输出的第2个码字

            // P(u_t) * P(y_t1 | c1) * P(y_t2 | c2)
            pLine[t][index] = 0.5 * pCh[2 * t][c1] * pCh[2 * t + 1][c2];
        }   
    }
    
	// 3. 计算前向概率 (Alpha)
    // Alpha_t(s) = sum_{s'} [ Alpha_{t-1}(s') * Gamma_t(s', s) ]
    pA[0][0] = 1.0; pA[0][1] = 0.0; pA[0][2] = 0.0; pA[0][3] = 0.0; // t=0 时刻，S0 概率为1
    double sumA;
    for (int t = 1; t < m_length + 1; t++)
    {
        sumA = 0;
        for (int st = 0; st < st_num; st++) // 遍历所有结束状态 st
        {
            // 累加所有能到达 st 的路径概率
            pA[t][st] = pA[t - 1][pathConn[st].first.point] * pLine[t - 1][pathConn[st].first.line] +
                        pA[t - 1][pathConn[st].second.point] * pLine[t - 1][pathConn[st].second.line];
            sumA += pA[t][st];
        }
        // 归一化，防止下溢
        if (sumA > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pA[t][st] = pA[t][st] / sumA;
            }
        }
    }
    
    // 4. 计算后向概率 (Beta)
    // Beta_t(s) = sum_{s''} [ Beta_{t+1}(s'') * Gamma_{t+1}(s, s'') ]
    pB[m_length][0] = 1.0; pB[m_length][1] = 0.0; pB[m_length][2] = 0.0; pB[m_length][3] = 0.0; // 终止时刻
    double sumB;
    for (int t = m_length - 1; t > -1; t--)
    {
        sumB = 0;
        for (int st = 0; st < st_num; st++) // 遍历所有起始状态 st
        {
            // 累加所有从 st 出发的路径概率
            pB[t][st] = pB[t + 1][pathConnB[st].first.point] * pLine[t][pathConnB[st].first.line] +
                        pB[t + 1][pathConnB[st].second.point] * pLine[t][pathConnB[st].second.line];
            sumB += pB[t][st];
        }
        // 归一化
        if (sumB > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pB[t][st] = pB[t][st] / sumB;
            }
        }
    }
    
    // 5. 计算后验概率 (APP) 并判决
    // LLR(u_t) = log [ P(u_t=0 | y) / P(u_t=1 | y) ]
    // P(u_t=i | y) = sum_{s', s s.t. input=i} [ Alpha_{t-1}(s') * Gamma_t(s',s) * Beta_t(s) ]
    
    double Po, Pz; // Po = P(u_t = 1 | y), Pz = P(u_t = 0 | y)
    for (int t = 0; t < m_length; t++)
	{
        // 累加所有 input = 0 的转移边的概率
		Pz = pA[t][0] * pLine[t][0] * pB[t+1][0]   // S0->S0 (边 1, input=0)
           + pA[t][1] * pLine[t][2] * pB[t+1][0]   // S1->S0 (边 3, input=0)
           + pA[t][2] * pLine[t][4] * pB[t+1][1]   // S2->S1 (边 5, input=0)
           + pA[t][3] * pLine[t][6] * pB[t+1][1];  // S3->S1 (边 7, input=0)

        // 累加所有 input = 1 的转移边的概率
		Po = pA[t][0] * pLine[t][1] * pB[t+1][2]   // S0->S2 (边 2, input=1)
           + pA[t][1] * pLine[t][3] * pB[t+1][2]   // S1->S2 (边 4, input=1)
           + pA[t][2] * pLine[t][5] * pB[t+1][3]   // S2->S3 (边 6, input=1)
           + pA[t][3] * pLine[t][7] * pB[t+1][3];  // S3->S3 (边 8, input=1)

        // MAP 判决：选择概率最大的
		if(Pz < Po) {
			de_message[t] = 1;
	    } else {
			de_message[t] = 0;
		}		
	}
}