#include "time.h"

// 定义全局变量来存储开始时间
clock_t tic_time;

// 函数tic用于记录开始时间
void tic() {
    tic_time = clock();
}

// 函数toc用于计算经过的时间并打印结果
double toc() {
    clock_t toc_time = clock();
    double elapsed_time = (double)(toc_time - tic_time) / CLOCKS_PER_SEC;
    
    return elapsed_time;
    //printf("Elapsed time: %.6f seconds\n", elapsed_time);
}