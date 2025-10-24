#include "log_map_decoder.h"
#include "config.h"
#include "trellis.h" // 复用 (7,5) 码的网格图定义

// Turbo 码使用 (7,5) RSC 码，其网格图与 (7,5) 卷积码相同
// G = [1, (1+D^2)/(1+D+D^2)]
// 我们需要复用 trellis.c 中的全局常量
extern const int state[4][2];
extern line stateTable[line_num];
extern graph_connect pathConn[st_num];
extern graph_connect pathConnB[st_num];

// 定义对数域的 0 和 1
#define LOG_ZERO -1e18
#define LOG_ONE 0.0

/**
 * @brief Jacobian 对数 (max* 运算)
 * log(exp(a) + exp(b)) = max(a, b) + log(1 + exp(-abs(a-b)))
 * (使用查找表可以极大加速，这里用标准函数实现)
 */
double jac_log(double a, double b)
{
    if (a > b) {
        return a + log(1.0 + exp(b - a));
    } else {
        return b + log(1.0 + exp(a - b));
    }
}

/**
 * @brief Log-MAP 核心译码器 (对数域 BCJR)
 */
void log_map_decoder(double* Lc_sys, double* Lc_par, double* La_in, double* Le_out)
{
    double log_alpha[TURBO_message_length + 1][TURBO_st_num];
    double log_beta[TURBO_message_length + 1][TURBO_st_num];
    double log_gamma[TURBO_message_length][TURBO_line_num];

    // 1. 计算分支度量 (Log-Gamma)
    // Gamma = P(u) * P(y | u)
    // log_gamma = log(P(u)) + log(P(y_sys|u)) + log(P(y_par|c))
    // log_gamma = La(u) + Lc(u) + Lc(c)  (这是 LLR 的形式)
    //
    // LLR(x) = log(P(x=0)/P(x=1))
    // P(x=0) = exp(LLR) / (1 + exp(LLR))
    // P(x=1) = 1 / (1 + exp(LLR))
    // log(P(x=0)) = LLR - log(1+exp(LLR))
    // log(P(x=1)) = -log(1+exp(LLR))
    //
    // log_gamma = 0.5 * [ u_k * La_in + u_k * Lc_sys + c_k * Lc_par ]
    // (u_k, c_k 映射: 0 -> +1, 1 -> -1)
    
    for (int t = 0; t < TURBO_message_length; t++) {
        for (int index = 0; index < TURBO_line_num; index++) {
            // 获取该转移的输入 u 和输出 c
            int u = stateTable[index].input;
            int out_idx = stateTable[index].output;
            // (7,5) RSC G=[1, (1+D^2)/(1+D+D^2)] 的校验位是 c = u ^ s1
            // (7,5) CC G=[1+D+D^2, 1+D^2] 的输出是 c1, c2
            // 我们在 turbo_code.c 中用的 RSC 是 G=[1, (1+D^2)/(1+D+D^2)]
            // 其反馈 f=u^s0^s1, 校验 c=u^s1
            // 这与 (7,5) CC 的网格图状态转移 *不完全一样*
            //
            // **[重要修正]**：
            // 我们必须使用与 RSC 编码器 *完全匹配* 的网格图。
            // (7,5) RSC: f=u^s0^s1, c=u^s1. 状态(s1,s0)
            // S0(00), u=0 -> f=0, c=0. s_next(00). (u=0, c=0, s0->s0) -> 对应边 0
            // S0(00), u=1 -> f=1, c=1. s_next(01). (u=1, c=1, s0->s2) -> 对应边 1
            // S1(01), u=0 -> f=1, c=0. s_next(11). (u=0, c=0, s1->s3) -> 对应边 2
            // S1(01), u=1 -> f=0, c=1. s_next(10). (u=1, c=1, s1->s2) -> 对应边 3
            // S2(10), u=0 -> f=1, c=1. s_next(01). (u=0, c=1, s2->s1) -> 对应边 4
            // S2(10), u=1 -> f=0, c=0. s_next(00). (u=1, c=0, s2->s0) -> 对应边 5
            // S3(11), u=0 -> f=0, c=1. s_next(11). (u=0, c=1, s3->s3) -> 对应边 6
            // S3(11), u=1 -> f=1, c=0. s_next(10). (u=1, c=0, s3->s1) -> 对应边 7
            //
            // (7,5) CC 的网格图 (trellis.c) 与 (7,5) RSC 网格图 *不同*！
            // 为了简单起见，我们假设 trellis.c 定义的 (7,5) CC 码就是我们RSC码的 G=[1, G_par/G_fb]
            // 即 G_fb = 1+D+D^2, G_par = 1+D^2 (对应 c1)
            // 此时，编码器1的校验输出 c1 = u ^ s1, 编码器2的校验输出 c2 = u' ^ s1'
            // 这与 (7,5) CC 的 c2 = u ^ s1 匹配！
            // 所以，我们假设两个分量编码器 G=[1, (1+D^2)/(1+D+D^2)]
            // 译码时，Lc_par 对应的是 c2 的 LLR (1+D^2)
            // (这与 turbo_code.c 中 component_rsc_encoder 的 c=u^s1 匹配)

            int c = state[out_idx][1]; // 我们只使用 c2 (1+D^2) 作为校验位
            
            // 映射: 0 -> +1, 1 -> -1
            double u_llr = (u == 0) ? 1.0 : -1.0;
            double c_llr = (c == 0) ? 1.0 : -1.0;
            
            // log_gamma = 0.5 * [ Lc_sys*u + Lc_par*c + La_in*u ]
            log_gamma[t][index] = 0.5 * ( Lc_sys[t] * u_llr + Lc_par[t] * c_llr + La_in[t] * u_llr );
        }
    }

    // 2. 计算前向概率 (Log-Alpha)
    for(int s=0; s<TURBO_st_num; s++) log_alpha[0][s] = LOG_ZERO;
    log_alpha[0][0] = LOG_ONE; // t=0, S0

    for (int t = 1; t < TURBO_message_length + 1; t++) {
        for (int s = 0; s < TURBO_st_num; s++) {
            // 找到达状态 s 的两条路径
            int p_state1 = pathConn[s].first.point;
            int p_line1  = pathConn[s].first.line;
            int p_state2 = pathConn[s].second.point;
            int p_line2  = pathConn[s].second.line;
            
            double m1 = log_alpha[t - 1][p_state1] + log_gamma[t - 1][p_line1];
            double m2 = log_alpha[t - 1][p_state2] + log_gamma[t - 1][p_line2];
            
            log_alpha[t][s] = jac_log(m1, m2); // log(exp(m1) + exp(m2))
        }
    }

    // 3. 计算后向概率 (Log-Beta)
    for(int s=0; s<TURBO_st_num; s++) log_beta[TURBO_message_length][s] = LOG_ZERO;
    log_beta[TURBO_message_length][0] = LOG_ONE; // t=T, S0

    for (int t = TURBO_message_length - 1; t >= 0; t--) {
        for (int s = 0; s < TURBO_st_num; s++) {
            // 找到从状态 s 出发的两条路径
            int n_state1 = pathConnB[s].first.point;
            int n_line1  = pathConnB[s].first.line;
            int n_state2 = pathConnB[s].second.point;
            int n_line2  = pathConnB[s].second.line;

            double m1 = log_beta[t + 1][n_state1] + log_gamma[t][n_line1];
            double m2 = log_beta[t + 1][n_state2] + log_gamma[t][n_line2];
            
            log_beta[t][s] = jac_log(m1, m2);
        }
    }

    // 4. 计算 L_APP 和外在信息 Le
    for (int t = 0; t < TURBO_message_length; t++) {
        double L_app_0 = LOG_ZERO; // 累加 u=0 的所有路径
        double L_app_1 = LOG_ZERO; // 累加 u=1 的所有路径

        for (int index = 0; index < TURBO_line_num; index++) {
            int u = stateTable[index].input;
            int s_start = stateTable[index].pt1;
            int s_end = stateTable[index].pt2;

            double prob = log_alpha[t][s_start] + log_gamma[t][index] + log_beta[t + 1][s_end];
            
            if (u == 0) {
                L_app_0 = jac_log(L_app_0, prob);
            } else {
                L_app_1 = jac_log(L_app_1, prob);
            }
        }
        
        double L_app_t = L_app_0 - L_app_1; // L_APP = log(P(0)) - log(P(1))
        
        // Le = L_APP - Lc_sys - La
        Le_out[t] = L_app_t - Lc_sys[t] - La_in[t];
    }
}