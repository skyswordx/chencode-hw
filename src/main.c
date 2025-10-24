/***************************************************
信道编码课程作业：卷积码 & Turbo 码
此程序模板提供了消息生成器、BPSK调制、AWGN信道模型和BPSK解调。
此文件是整合了 Viterbi (硬/软)、BCJR 和 Turbo 译码器的完整仿真框架。
***************************************************/

#define  _CRT_SECURE_NO_WARNINGS // 兼容 MSVC

// 包含所有自定义模块的头文件
#include "config.h"
#include "convolutional_code.h"
#include "trellis.h"
#include "viterbi.h"
#include "bcjr.h"
#include "turbo_code.h" // <-- 新增

// 码率 (根据译码器类型在 main 中设置)
float code_rate;

// 编码器状态数 (在 main 中设置)
int state_num;

// --- 函数原型声明 ---
void statetable(); // 状态表生成 (在此项目中未使用，使用硬编码)
void decoder();    // 卷积码译码器总入口

/**
 * @brief 卷积码 (CC) 仿真流程
 */
void run_cc_simulation(float start_snr, float finish_snr, long int seq_num)
{
    int i;
	long int bit_error, seq;
	double BER;
	double progress;

    // CC 码率 (仅计算信息比特)
    code_rate = (float)MESSAGE_BITS / (float)codeword_length;

    printf("=================================================================\n");
    #if DECODER_TYPE == 1
        printf("               当前译码器: 硬判决 Viterbi\n");
    #elif DECODER_TYPE == 2
        printf("               当前译码器: 软判决 Viterbi\n");
    #elif DECODER_TYPE == 3
        printf("               当前译码器: BCJR (MAP)\n");
    #endif
    printf("               码型: (7,5)_8 卷积码, R=1/2\n");
    printf("=================================================================\n");

    // 按 SNR 循环
	for (float SNR = start_snr; SNR <= finish_snr; SNR++)
	{
		// 计算信道噪声参数
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2.0); 
		
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

			progress = (double)(seq * 100) / (double)seq_num; 
			BER = (double)bit_error / (double)(MESSAGE_BITS * seq);

			printf("进度=%5.1f%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\r", progress, SNR, bit_error, BER);
		}

		BER = (double)bit_error / (double)(MESSAGE_BITS * seq_num);
		printf("进度=100.0%%, SNR=%4.1f dB, 误比特=%8ld, BER=%E\n", SNR, bit_error, BER);
	}
}


/**
 * @brief 主函数
 */
int main()
{
	float start, finish;
	long int seq_num;

    // 设置编码器状态数 (用于 CC 消息补零)
    state_num = STATE_MEM;

	// 生成状态表 (未使用)
	statetable();

	// 设置随机数种子
	srand((int)time(0));

	// 输入 SNR 范围和仿真帧数
	printf("\n请输入起始 SNR (dB): ");
	scanf("%f", &start);
	printf("\n请输入结束 SNR (dB): ");
	scanf("%f", &finish);
	printf("\n请输入要仿真的消息帧数: ");
	scanf("%ld", &seq_num); // <-- 修复 Bug：使用 %ld 读取 long int

    // 根据选择的译码器类型，执行不同的仿真流程
    #if DECODER_TYPE <= 3
        // 执行卷积码仿真
        run_cc_simulation(start, finish, seq_num);
    
    #elif DECODER_TYPE == 4
        // 执行 Turbo 码仿真
        run_turbo_simulation(start, finish, seq_num);
    
    #endif

	// 修复 Bug：使用 getchar() 替换 system("pause")，提高兼容性
    printf("仿真完成，请按 Enter 键退出...\n");
    while (getchar() != '\n'); // 清空输入缓冲区
    getchar(); // 等待用户按键
    
    return 0;
}

/**
 * @brief 状态表生成函数 (空)
 */
void statetable()
{
    // 空
}

/**
 * @brief 卷积码译码器总入口
 */
void decoder()
{ 
    #if DECODER_TYPE == 1
        hardDecoder(re_codeword, de_message, message_length);
    
    #elif DECODER_TYPE == 2
        softDecode(rx_symbol, de_message, message_length);
    
    #elif DECODER_TYPE == 3
        BCJR(rx_symbol, de_message, message_length, codeword_length);
    
    #else
       // Turbo 码不使用此函数
    #endif
}