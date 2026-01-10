#include "csv_export.h"
#include "config.h"
#include <time.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(dir) _mkdir(dir)
#else
#include <sys/stat.h>
#define MKDIR(dir) mkdir(dir, 0755)
#endif

/**
 * @brief 确保输出目录存在
 */
static void ensure_output_dir() {
    MKDIR("output");
}

FILE* csv_init(const char* filename, const char* decoder_name,
               double code_rate, int block_size) {
    ensure_output_dir();
    
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "[CSV] 无法创建文件: %s\n", filename);
        return NULL;
    }
    
    // 写入元数据注释 (MATLAB readtable 支持 # 注释)
    time_t now = time(NULL);
    struct tm* t = localtime(&now);
    
    fprintf(fp, "# ============================================\n");
    fprintf(fp, "# Channel Coding Simulation Results\n");
    fprintf(fp, "# ============================================\n");
    fprintf(fp, "# Decoder: %s\n", decoder_name);
    fprintf(fp, "# Code Rate: %.4f\n", code_rate);
    fprintf(fp, "# Block Size: %d bits\n", block_size);
    fprintf(fp, "# Generated: %04d-%02d-%02d %02d:%02d:%02d\n",
            t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
            t->tm_hour, t->tm_min, t->tm_sec);
    fprintf(fp, "# ============================================\n");
    
    // CSV 表头
    fprintf(fp, "SNR_dB,Bit_Errors,Total_Bits,BER\n");
    
    return fp;
}

void csv_append_row(FILE* fp, double snr, long bit_errors,
                    long total_bits, double ber) {
    if (!fp) return;
    fprintf(fp, "%.1f,%ld,%ld,%.6e\n", snr, bit_errors, total_bits, ber);
    fflush(fp); // 实时刷新，防止崩溃丢失数据
}

void csv_close(FILE* fp) {
    if (fp) {
        fflush(fp);
        fclose(fp);
    }
}

void csv_generate_filename(char* buffer, int buf_size, int decoder_type) {
    time_t now = time(NULL);
    struct tm* t = localtime(&now);
    
    const char* type_names[] = {
        "uncoded", "hard_viterbi", "soft_viterbi", "bcjr", "turbo"
    };
    
    const char* type_name = (decoder_type >= 0 && decoder_type <= 4) 
                            ? type_names[decoder_type] 
                            : "unknown";
    
    snprintf(buffer, buf_size, "output/ber_%s_%04d%02d%02d_%02d%02d%02d.csv",
             type_name,
             t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
             t->tm_hour, t->tm_min, t->tm_sec);
}
