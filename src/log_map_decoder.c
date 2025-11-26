#include "log_map_decoder.h"
#include "config.h"
#include "trellis.h" 
#include <math.h>

// 定义对数域的“零概率” (即 log(0) -> -inf)
// 使用一个足够小的负数，但不要让它在加法中导致下溢失效
#define LOG_ZERO -1.0e10 
#define LOG_ONE 0.0

// (7,5)_8 RSC 专用网格图
// 状态映射: Index 0(00), 1(10), 2(01), 3(11)  <-- 注意这里的二进制是 (s1, s0)
// 您的编码器逻辑: f = u^s0^s1; out = f^s1; s1=s0; s0=f;
// 经过验证，您的 trellis 表逻辑是正确的，可以直接复用。
static const line stateTable_rsc_v2[line_num] = {
    {0, 0, 0, 0, 1}, 
    {1, 1, 0, 1, 2}, 
    {0, 1, 1, 3, 3}, 
    {1, 0, 1, 2, 4}, 
    {0, 0, 2, 1, 5}, 
    {1, 1, 2, 0, 6}, 
    {0, 1, 3, 2, 7}, 
    {1, 0, 3, 3, 8}  
};
static const graph_connect pathConn_rsc_v2[st_num] = {
    { {0,0}, {2,5} }, 
    { {0,1}, {2,4} }, 
    { {1,3}, {3,6} }, 
    { {1,2}, {3,7} }  
};
static const graph_connect pathConnB_rsc_v2[st_num] = {
    { {0,0}, {1,1} }, 
    { {3,2}, {2,3} }, 
    { {1,4}, {0,5} }, 
    { {2,6}, {3,7} }  
};

/**
 * @brief Jacobian Logarithm 实现
 * max*(a, b) = max(a, b) + log(1 + exp(-|a-b|))
 */
double jac_log(double a, double b)
{
    if (a <= LOG_ZERO) return b;
    if (b <= LOG_ZERO) return a;
    
    double diff = a - b;
    // 优化: 当差值很大时，修正项趋近于0，可忽略以提高速度
    if (diff > 10.0) return a;
    if (diff < -10.0) return b;
    
    return (a > b ? a : b) + log(1.0 + exp(-fabs(diff)));
}

/**
 * @brief Log-MAP 核心译码器 
 */
void log_map_decoder(double* Lc_sys, double* Lc_par, double* La_in, double* Le_out)
{
    // 栈内存优化：如果在嵌入式或栈受限环境，这些大数组应改为静态或动态分配
    // 但在PC仿真上通常没问题
    static double log_alpha[TURBO_message_length + 1][st_num];
    static double log_beta[TURBO_message_length + 1][st_num];
    static double log_gamma[TURBO_message_length][line_num];

    // 1. 计算分支度量 (Log-Gamma)
    for (int t = 0; t < TURBO_message_length; t++) {
        for (int index = 0; index < line_num; index++) {
            int u = stateTable_rsc_v2[index].input;
            int c = stateTable_rsc_v2[index].output;
            
            // 调制映射: 0->+1, 1->-1 
            // LLR 定义: ln(P(x=+1)/P(x=-1))
            // Gamma = 0.5 * (u_llr * u_soft + c_llr * c_soft)
            double u_val = (u == 0) ? 1.0 : -1.0;
            double c_val = (c == 0) ? 1.0 : -1.0;
            
            // 注意：这里不需要再乘0.5，因为传入的 Lc 已经是 LLR (4*Es/No * y)。
            // 公式推导中，gamma = 0.5 * x * LLR
            log_gamma[t][index] = 0.5 * ( (Lc_sys[t] + La_in[t]) * u_val + Lc_par[t] * c_val );
        }
    }

    // 2. 计算前向概率 (Log-Alpha)
    // 初始化：起始状态必须是0
    for(int s=0; s<st_num; s++) log_alpha[0][s] = LOG_ZERO;
    log_alpha[0][0] = LOG_ONE; 

    for (int t = 1; t < TURBO_message_length + 1; t++) {
        for (int s = 0; s < st_num; s++) {
            int prev_s1 = pathConn_rsc_v2[s].first.point;
            int line1   = pathConn_rsc_v2[s].first.line;
            int prev_s2 = pathConn_rsc_v2[s].second.point;
            int line2   = pathConn_rsc_v2[s].second.line;
            
            double m1 = log_alpha[t - 1][prev_s1] + log_gamma[t - 1][line1];
            double m2 = log_alpha[t - 1][prev_s2] + log_gamma[t - 1][line2];
            
            // 使用 jac_log 代替 max_log
            log_alpha[t][s] = jac_log(m1, m2); 
        }
        
        // Alpha 归一化 (防止数值溢出)
        // 每一时刻减去状态0的值 (或其他任意状态的值)
        double norm_val = log_alpha[t][0];
        for(int s=0; s<st_num; s++) {
             // 避免减去 LOG_ZERO (虽然不太可能发生，除非SNR极高且路径断裂)
             if(log_alpha[t][s] > LOG_ZERO)
                log_alpha[t][s] -= norm_val;
        }
    }

    // 3. 计算后向概率 (Log-Beta)
    // 【关键修复】: 由于编码器没有做 trellis termination (归零)，
    // 尾比特仅仅是补0数据，这不能保证结束状态是0。
    // 因此，Beta 初始化必须设为“所有状态等概率”。
    for(int s=0; s<st_num; s++) {
        log_beta[TURBO_message_length][s] = LOG_ONE; // log(1) = 0, 表示等概率
    }

    for (int t = TURBO_message_length - 1; t >= 0; t--) {
        for (int s = 0; s < st_num; s++) {
            int next_s1 = pathConnB_rsc_v2[s].first.point;
            int line1   = pathConnB_rsc_v2[s].first.line;
            int next_s2 = pathConnB_rsc_v2[s].second.point;
            int line2   = pathConnB_rsc_v2[s].second.line;

            double m1 = log_beta[t + 1][next_s1] + log_gamma[t][line1];
            double m2 = log_beta[t + 1][next_s2] + log_gamma[t][line2];
            
            // 使用 jac_log
            log_beta[t][s] = jac_log(m1, m2);
        }
        
        // Beta 归一化
        double norm_val = log_beta[t][0];
        for(int s=0; s<st_num; s++) {
             if(log_beta[t][s] > LOG_ZERO)
                log_beta[t][s] -= norm_val;
        }
    }

    // 4. 计算外在信息 Le (L_ext)
    for (int t = 0; t < TURBO_message_length; t++) {
        double L_app_0 = LOG_ZERO; 
        double L_app_1 = LOG_ZERO; 
        
        // 遍历所有边，累加概率到 u=0 或 u=1 的集合中
        for (int index = 0; index < line_num; index++) {
            int u = stateTable_rsc_v2[index].input;
            int s_start = stateTable_rsc_v2[index].pt1;
            int s_end = stateTable_rsc_v2[index].pt2;
            
            // Log-MAP: prob = alpha + gamma + beta
            // 注意：这里 gamma 包含先验信息(La)和系统位信道信息(Lc_sys)。
            // 计算 Extrinsic 时，后面会减去它们。
            double prob = log_alpha[t][s_start] + log_gamma[t][index] + log_beta[t + 1][s_end];
            
            if (u == 0) {
                L_app_0 = jac_log(L_app_0, prob);
            } else {
                L_app_1 = jac_log(L_app_1, prob);
            }
        }
        
        // L_APP = log(P(u=0)/P(u=1)) = L_app_0 - L_app_1
        // 注意：我们的 gamma 计算方式使用了 input 0 -> +1, 1 -> -1
        // 所以如果 L_app_0 大，说明 0 的概率大。
        // 标准 LLR 定义为 log(P(0)/P(1))。
        double L_all = L_app_0 - L_app_1;
        
        // Le = L_all - (Lc_sys + La_in)
        // 这里的 Lc_sys 和 La_in 已经在 gamma 中加过了，所以要减去
        Le_out[t] = L_all - (Lc_sys[t] + La_in[t]);
        
        // 限制幅度，防止数值发散
        if (Le_out[t] > 50.0) Le_out[t] = 50.0;
        if (Le_out[t] < -50.0) Le_out[t] = -50.0;
    }
}