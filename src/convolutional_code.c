#include "convolutional_code.h"
#include "config.h"
#include "trellis.h" // 需要 trellis.h 里的 N0/sgm

// --- 全局数组定义 ---
int message[message_length];
int codeword[codeword_length];
int re_codeword[codeword_length];
int de_message[message_length];

double tx_symbol[codeword_length][softIn_st_num];
double rx_symbol[codeword_length][softIn_st_num];


/**
 * @brief (7,5)_8 卷积编码器
 * 对应 G(D) = [1+D+D^2, 1+D^2] -> [111, 101] -> (7, 5) 八进制
 * 状态寄存器 s = (s1, s0)
 * s0_next = message[i]
 * s1_next = s0
 * 输出 c1 = m[i] ^ s0 ^ s1  (对应 111)
 * 输出 c2 = m[i] ^ s1       (对应 101)
 * 这与 trellis.c 中定义的 (7,5) 网格图一致。
 */
void encoder()
{
    int c1, c2;
    int s0 = 0; // 寄存器 1 (D)
    int s1 = 0; // 寄存器 2 (D^2)
    for (int i = 0; i < message_length; i++)
    {
        c1 = message[i] ^ s0 ^ s1; // G_c1 = 1 + D + D^2
        c2 = message[i] ^ s1;      // G_c2 = 1 + D^2
        
        s1 = s0;         // 移位
        s0 = message[i]; // 移位
        codeword[2 * i] = c1;
        codeword[2 * i + 1] = c2;
    }
}

/**
 * @brief BPSK 调制
 */
void modulation()
{
	int i;
	// 0 映射到 (1, 0)
	// 1 映射到 (-1, 0)
	for (i = 0; i<codeword_length; i++)
	{
		tx_symbol[i][0] = -1.0 * (2.0 * codeword[i] - 1.0); // 0->1, 1->-1
		tx_symbol[i][1] = 0.0;
	}
}

/**
 * @brief AWGN 信道模型
 */
void channel()
{
	int i, j;
	double u, r, g;

	for (i = 0; i<codeword_length; i++)
	{
		for (j = 0; j<softIn_st_num; j++) // 对 I, Q 两路加噪声 (虽然Q路是0)
		{
            // Box-Muller 变换生成高斯分布
			u = (double)rand() / (double)RAND_MAX;
			if (u == 1.0) u = 0.999999;
			r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

			u = (double)rand() / (double)RAND_MAX;
			if (u == 1.0) u = 0.999999;
			g = (double)r * cos(2.0 * pi * u);

			rx_symbol[i][j] = tx_symbol[i][j] + g; // 接收符号 = 发送符号 + 噪声
		}
	}
}

/**
 * @brief BPSK 硬判决解调
 */
void demodulation()
{
	int i;
	double d1, d2;
	for (i = 0; i<codeword_length; i++)
	{
        // 0 对应 +1, 1 对应 -1
		// d1 = 接收符号到 +1 (即 0) 的欧氏距离的平方
		d1 = pow(rx_symbol[i][0] - 1.0, 2) + pow(rx_symbol[i][1] - 0.0, 2);
		// d2 = 接收符号到 -1 (即 1) 的欧氏距离的平方
		d2 = pow(rx_symbol[i][0] + 1.0, 2) + pow(rx_symbol[i][1] - 0.0, 2);
		
        if (d1 < d2)
			re_codeword[i] = 0; // 离 +1 更近，判为 0
		else
			re_codeword[i] = 1; // 离 -1 更近，判为 1
	}
}