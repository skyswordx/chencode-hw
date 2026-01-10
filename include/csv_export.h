#ifndef CSV_EXPORT_H_
#define CSV_EXPORT_H_

#include <stdio.h>

/**
 * @brief 初始化 CSV 文件，写入元数据和表头
 * @param filename 文件路径
 * @param decoder_name 译码器名称
 * @param code_rate 码率
 * @param block_size 帧长
 * @return FILE* 文件指针，失败返回 NULL
 */
FILE* csv_init(const char* filename, const char* decoder_name, 
               double code_rate, int block_size);

/**
 * @brief 追加一行 BER 数据到 CSV
 * @param fp 文件指针
 * @param snr 信噪比 (dB)
 * @param bit_errors 误比特数
 * @param total_bits 总比特数
 * @param ber 误码率
 */
void csv_append_row(FILE* fp, double snr, long bit_errors, 
                    long total_bits, double ber);

/**
 * @brief 关闭 CSV 文件
 * @param fp 文件指针
 */
void csv_close(FILE* fp);

/**
 * @brief 生成带时间戳的 CSV 文件名
 * @param buffer 输出缓冲区
 * @param buf_size 缓冲区大小
 * @param decoder_type 译码器类型编号
 */
void csv_generate_filename(char* buffer, int buf_size, int decoder_type);

#endif // CSV_EXPORT_H_
