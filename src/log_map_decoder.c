#include "log_map_decoder.h"
#include "config.h"
#include "trellis.h" // 复用 trellis.h 中的结构体定义

// 定义对数域的 0 和 1
#define LOG_ZERO -1e18
#define LOG_ONE 0.0

// (7,5)_8 RSC 专用网格图 (与上一版相同)
static const line stateTable_rsc_v2[line_num] = {
    {0, 0, 0, 0, 1}, // s0(00), u=0 -> c=0, s_next=s0(00)
    {1, 1, 0, 1, 2}, // s0(00), u=1 -> c=1, s_next=s1(10)
    {0, 1, 1, 3, 3}, // s1(10), u=0 -> c=1, s_next=s3(11)
    {1, 0, 1, 2, 4}, // s1(10), u=1 -> c=0, s_next=s2(01)
    {0, 0, 2, 1, 5}, // s2(01), u=0 -> c=0, s_next=s1(10)
    {1, 1, 2, 0, 6}, // s2(01), u=1 -> c=1, s_next=s0(00)
    {0, 1, 3, 2, 7}, // s3(11), u=0 -> c=1, s_next=s2(01)
    {1, 0, 3, 3, 8}  // s3(11), u=1 -> c=0, s_next=s3(11)
};
static const graph_connect pathConn_rsc_v2[st_num] = {
    { {0,0}, {2,5} }, // 到达 S0
    { {0,1}, {2,4} }, // 到达 S1
    { {1,3}, {3,6} }, // 到达 S2
    { {1,2}, {3,7} }  // 到达 S3
};
static const graph_connect pathConnB_rsc_v2[st_num] = {
    { {0,0}, {1,1} }, // 从 S0 出发
    { {3,2}, {2,3} }, // 从 S1 出发
    { {1,4}, {0,5} }, // 从 S2 出发
    { {2,6}, {3,7} }  // 从 S3 出发
};


/**
 * @brief Max-Log 近似 (Max-Log-MAP)
 */
static double max_log(double a, double b)
{
    return (a > b) ? a : b;
}


/**
 * @brief Log-MAP 核心译码器 (已修正为在 LLR 处归一化)
 */
void log_map_decoder(double* Lc_sys, double* Lc_par, double* La_in, double* Le_out)
{
    double log_alpha[TURBO_message_length + 1][st_num];
    double log_beta[TURBO_message_length + 1][st_num];
    double log_gamma[TURBO_message_length][line_num];

    // 1. 计算分支度量 (Log-Gamma)
    // ** [最终修复] **：不再对 Gamma 归一化！
    for (int t = 0; t < TURBO_message_length; t++) {
        for (int index = 0; index < line_num; index++) {
            int u = stateTable_rsc_v2[index].input;
            int c = stateTable_rsc_v2[index].output;
            double u_llr = (u == 0) ? 1.0 : -1.0;
            double c_llr = (c == 0) ? 1.0 : -1.0;
            
            log_gamma[t][index] = 0.5 * ( (Lc_sys[t] + La_in[t]) * u_llr + Lc_par[t] * c_llr );
        }
    }

    // 2. 计算前向概率 (Log-Alpha)
    // ** [最终修复] **：不再归一化，让它自由溢出
    for(int s=0; s<st_num; s++) log_alpha[0][s] = LOG_ZERO;
    log_alpha[0][0] = LOG_ONE; 

    for (int t = 1; t < TURBO_message_length + 1; t++) {
        for (int s = 0; s < st_num; s++) {
            int p_state1 = pathConn_rsc_v2[s].first.point;
            int p_line1  = pathConn_rsc_v2[s].first.line;
            int p_state2 = pathConn_rsc_v2[s].second.point;
            int p_line2  = pathConn_rsc_v2[s].second.line;
            
            double m1 = log_alpha[t - 1][p_state1] + log_gamma[t - 1][p_line1];
            double m2 = log_alpha[t - 1][p_state2] + log_gamma[t - 1][p_line2];
            
            log_alpha[t][s] = max_log(m1, m2); 
        }
    }

    // 3. 计算后向概率 (Log-Beta)
    // ** [最终修复] **：不再归一化，让它自由溢出
    for(int s=0; s<st_num; s++) log_beta[TURBO_message_length][s] = LOG_ZERO;
    log_beta[TURBO_message_length][0] = LOG_ONE; 

    for (int t = TURBO_message_length - 1; t >= 0; t--) {
        for (int s = 0; s < st_num; s++) {
            int n_state1 = pathConnB_rsc_v2[s].first.point;
            int n_line1  = pathConnB_rsc_v2[s].first.line;
            int n_state2 = pathConnB_rsc_v2[s].second.point;
            int n_line2  = pathConnB_rsc_v2[s].second.line;

            double m1 = log_beta[t + 1][n_state1] + log_gamma[t][n_line1];
            double m2 = log_beta[t + 1][n_state2] + log_gamma[t][n_line2];
            
            log_beta[t][s] = max_log(m1, m2);
        }
    }

    // 4. 计算 L_APP 和外在信息 Le
    for (int t = 0; t < TURBO_message_length; t++) {
        double L_app_0 = LOG_ZERO; 
        double L_app_1 = LOG_ZERO; 
        
        // ** [最终修复：在此处归一化] **
        // 1. 找出所有路径中的最大概率值
        double norm_LLR = LOG_ZERO;
        for (int index = 0; index < line_num; index++) {
            int s_start = stateTable_rsc_v2[index].pt1;
            int s_end = stateTable_rsc_v2[index].pt2;
            double prob = log_alpha[t][s_start] + log_gamma[t][index] + log_beta[t + 1][s_end];
            norm_LLR = max_log(norm_LLR, prob);
        }

        // 2. 计算减去最大值之后的 L_app_0 和 L_app_1
        for (int index = 0; index < line_num; index++) {
            int u = stateTable_rsc_v2[index].input;
            int s_start = stateTable_rsc_v2[index].pt1;
            int s_end = stateTable_rsc_v2[index].pt2;

            // prob_norm = (log_alpha + log_gamma + log_beta) - norm_LLR
            double prob_norm = log_alpha[t][s_start] + log_gamma[t][index] + log_beta[t + 1][s_end] - norm_LLR;
            
            if (u == 0) {
                L_app_0 = max_log(L_app_0, prob_norm);
            } else {
                L_app_1 = max_log(L_app_1, prob_norm);
            }
        }
        
        // L_app_t = L_app_0_norm - L_app_1_norm
        // (L_app_0_norm = log_sum_exp(prob_u0) - norm_LLR)
        // (L_app_1_norm = log_sum_exp(prob_u1) - norm_LLR)
        // L_app_t = (log_sum_exp(prob_u0) - norm_LLR) - (log_sum_exp(prob_u1) - norm_LLR)
        // L_app_t = log_sum_exp(prob_u0) - log_sum_exp(prob_u1)
        // 这个归一化是正确的，因为它在 LLR 的减法中被消除了。
        double L_app_t = L_app_0 - L_app_1; 
        
        Le_out[t] = L_app_t - Lc_sys[t] - La_in[t];
    }
}