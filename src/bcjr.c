#include "bcjr.h"
#include "config.h"
#include "trellis.h" // BCJR 算法需要网格图定义

// --- BCJR 译码器所需的全局数组 ---
// 存放前向、后向和分支概率
double pCh[CC_codeword_length][2];     // 信道观察概率 P(y|x) [x=0 或 x=1]
double pLine[CC_message_length][line_num]; // 分支转移概率 Gamma
double pA[CC_message_length + 1][st_num];  // 前向概率 Alpha
double pB[CC_message_length + 1][st_num];  // 后向概率 Beta


/**
 * @brief 计算欧氏距离的平方 (静态辅助函数)
 */
static double euDist(double rx_x, double rx_y, int sym)
{
    double ideal_x = (sym == 0) ? 1.0 : -1.0;
    return pow(rx_x - ideal_x, 2) + pow(rx_y - 0.0, 2);
}

/**
 * @brief 计算信道观察概率 P(y|x) (静态辅助函数)
 */
static double chObs(double eu_sq, double n0)
{
    return exp(-eu_sq / n0) / (sqrt(pi * n0));
}


/**
 * @brief BCJR 译码器主函数
 */
void BCJR(double rx_sym[][softIn_st_num], int de_message[], int m_length, int c_length)
{
	// 1. 计算信道观察概率 P(y|x)
    double sum;
    for (int t = 0; t < c_length; t++)
    {
        sum = 0;
        pCh[t][0] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 0), N0);
        pCh[t][1] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 1), N0);
        
        sum = pCh[t][0] + pCh[t][1];
        if (sum > 1e-300) { 
            pCh[t][0] = pCh[t][0] / sum;
            pCh[t][1] = pCh[t][1] / sum;
        } else {
            pCh[t][0] = 0.5;
            pCh[t][1] = 0.5;
        }
    }
    
    // 2. 计算分支转移概率 (Gamma)
    for (int t = 0; t < m_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            // 过滤不可能的转移
            if (t == 0 && index != 0 && index != 1) { pLine[t][index] = p0; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { pLine[t][index] = p0; continue; }
            if (t == m_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { pLine[t][index] = p0; continue; }
            if (t == m_length - 1 && index != 0 && index != 2) { pLine[t][index] = p0; continue; }
            
            int out_idx = stateTable[index].output; 
            int c1 = state[out_idx][0]; 
            int c2 = state[out_idx][1]; 

            pLine[t][index] = 0.5 * pCh[2 * t][c1] * pCh[2 * t + 1][c2];
        }   
    }
    
	// 3. 计算前向概率 (Alpha)
    pA[0][0] = 1.0; pA[0][1] = 0.0; pA[0][2] = 0.0; pA[0][3] = 0.0; 
    double sumA;
    for (int t = 1; t < m_length + 1; t++)
    {
        sumA = 0;
        for (int st = 0; st < st_num; st++) 
        {
            pA[t][st] = pA[t - 1][pathConn[st].first.point] * pLine[t - 1][pathConn[st].first.line] +
                        pA[t - 1][pathConn[st].second.point] * pLine[t - 1][pathConn[st].second.line];
            sumA += pA[t][st];
        }
        if (sumA > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pA[t][st] = pA[t][st] / sumA;
            }
        }
    }
    
    // 4. 计算后向概率 (Beta)
    pB[m_length][0] = 1.0; pB[m_length][1] = 0.0; pB[m_length][2] = 0.0; pB[m_length][3] = 0.0; 
    double sumB;
    for (int t = m_length - 1; t > -1; t--)
    {
        sumB = 0;
        for (int st = 0; st < st_num; st++) 
        {
            pB[t][st] = pB[t + 1][pathConnB[st].first.point] * pLine[t][pathConnB[st].first.line] +
                        pB[t + 1][pathConnB[st].second.point] * pLine[t][pathConnB[st].second.line];
            sumB += pB[t][st];
        }
        if (sumB > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pB[t][st] = pB[t][st] / sumB;
            }
        }
    }
    
    // 5. 计算后验概率 (APP) 并判决
    double Po, Pz; 
    for (int t = 0; t < m_length; t++)
	{
		Pz = pA[t][0] * pLine[t][0] * pB[t+1][0]   
           + pA[t][1] * pLine[t][2] * pB[t+1][0]   
           + pA[t][2] * pLine[t][4] * pB[t+1][1]   
           + pA[t][3] * pLine[t][6] * pB[t+1][1];  

		Po = pA[t][0] * pLine[t][1] * pB[t+1][2]   
           + pA[t][1] * pLine[t][3] * pB[t+1][2]   
           + pA[t][2] * pLine[t][5] * pB[t+1][3]   
           + pA[t][3] * pLine[t][7] * pB[t+1][3];  

		if(Pz < Po) {
			de_message[t] = 1;
	    } else {
			de_message[t] = 0;
		}		
	}
}