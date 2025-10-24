/***************************************************
信道编码课程作业：卷积码
此程序模板提供了消息生成器、BPSK调制、AWGN信道模型和BPSK解调。
此文件是整合了 Viterbi (硬/软) 和 BCJR 译码器的完整仿真框架。
码型: (7,5)_8 卷积码, 码率 R=1/2, 4状态 (m=2)
***************************************************/

#define  _CRT_SECURE_NO_WARNINGS // 兼容 MSVC

// 包含所有自定义模块的头文件
#include "config.h"
#include "convolutional_code.h"
#include "trellis.h"
#include "viterbi.h"
#include "bcjr.h"

// 码率计算 (仅计算信息比特)
float code_rate = (float)MESSAGE_BITS / (float)codeword_length;

// 编码器状态数 (在 main 函数中设置)
int state_num;

// --- 函数原型声明 ---
void statetable(); // 状态表生成 (在此项目中未使用，使用硬编码)
void decoder();    // 译码器总入口


/**
 * @brief 主函数
 */
int main()
{
	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;

    // 设置编码器状态数 (用于消息补零)
    state_num = STATE_MEM;

	// 生成状态表 (未使用，因为我们硬编码了 trellis.c)
	statetable();

	// 设置随机数种子
	srand((int)time(0));

	// 输入 SNR 范围和仿真帧数
	printf("\n请输入起始 SNR (dB): ");
	scanf("%f", &start);
	printf("\n请输入结束 SNR (dB): ");
	scanf("%f", &finish);
	printf("\n请输入要仿真的消息帧数: ");
	scanf("%ld", &seq_num);

    printf("=================================================================\n");
    #if DECODER_TYPE == 1
        printf("               当前译码器: 硬判决 Viterbi\n");
    #elif DECODER_TYPE == 2
        printf("               当前译码器: 软判决 Viterbi\n");
    #elif DECODER_TYPE == 3
        printf("               当前译码器: BCJR (MAP)\n");
    #endif
    printf("=================================================================\n");


	// 按 SNR 循环
	for (SNR = start; SNR <= finish; SNR++)
	{
		// 计算信道噪声参数
        // N0 = Eb/N0 * R = 1 / (R * 10^(SNR/10))
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2.0); // 高斯噪声标准差
		
		bit_error = 0; // 误比特数清零

        // 按帧数循环
		for (seq = 1; seq <= seq_num; seq++)
		{
			// 1. 随机生成二进制消息
			for (i = 0; i < MESSAGE_BITS; i++)
			{
				message[i] = rand() % 2;
			}
            // 2. 补零 (尾比特)，用于 trellis 终止
			for (i = MESSAGE_BITS; i < message_length; i++)
			{
				message[i] = 0;
			}

			// 3. 卷积编码
			encoder();

			// 4. BPSK 调制
			modulation();

			// 5. AWGN 信道
			channel();

            // 6. BPSK 解调 (仅在硬判决时需要)
            #if DECODER_TYPE == 1
			    demodulation();
            #endif

			// 7. 卷积译码
			decoder();

			// 8. 计算误比特数
            // 只比较 MESSAGE_BITS 个信息比特，不比较尾比特
			for (i = 0; i < MESSAGE_BITS; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num; // 计算进度
            // 计算当前 BER
			BER = (double)bit_error / (double)(MESSAGE_BITS * seq);

            // 打印实时仿真结果 (使用 \r 回车符覆盖上一行)
			printf("进度=%5.1f%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\r", progress, SNR, bit_error, BER);
		}

		// 计算当前 SNR 点的最终 BER
		BER = (double)bit_error / (double)(MESSAGE_BITS * seq_num);
        // 打印最终结果 (使用 \n 换行)
		printf("进度=100.0%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\n", SNR, bit_error, BER);
	}
	system("pause"); // 保持窗口打开
    return 0;
}

/**
 * @brief 状态表生成函数 (空)
 * 因为我们使用了 trellis.c 中的硬编码常量表，所以此函数无需执行任何操作。
 */
void statetable()
{
    // 空
}

/**
 * @brief 译码器总入口
 * 根据 config.h 中定义的 DECODER_TYPE 来选择调用哪个译码函数。
 */
void decoder()
{ 
    #if DECODER_TYPE == 1
        // 调用硬判决 Viterbi
        hardDecoder(re_codeword, de_message, message_length);
    
    #elif DECODER_TYPE == 2
        // 调用软判决 Viterbi
        softDecode(rx_symbol, de_message, message_length);
    
    #elif DECODER_TYPE == 3
        // 调用 BCJR
        BCJR(rx_symbol, de_message, message_length, codeword_length);
    
    #else
        // 错误处理
        printf("错误: 未知的 DECODER_TYPE! 请检查 config.h。\n");
        exit(1); // 退出程序
    #endif
}