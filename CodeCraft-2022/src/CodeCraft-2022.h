#pragma once 

#include "stdc++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include <atomic>
#include <thread>   

using namespace std;

#define __ONLINE__ 0

#define MAX_CLIENTNODE_NUM 35 // 最大客户节点数量
#define MAX_EDGENODE_NUM 135  // 最大边缘节点数量
#define MAX_T_LEN 8950        // 最大时间序列长度
#define MAX_STREAM_TYPES 100  // 最大流数量
#define CPU_CORES 1           // 设置的CPU核心数
#define MAX_CPU_CORES 4       // 支持的最大CPU核心数 (注意内存)
#define TRY_SHUT_DOWN 1       // 是否尝试关机以搜索更优解
#define MAX_FUN_PTR_NUM 100   // 函数指针数组最大长度
#define GLOBAL_MAX_TIME_MS (275 * 1000) // 运行时间限制(ms)

#define MAX_CLIENT_DEMAND 550000
#define MAX_EDGE_BANDWIDTH 1000000

#define CACHE_P1 0.05         // 决赛新增: 常规边缘节点的缓存保留率
#define CACHE_P2 0.01         // 决赛新增: 被选择的20台节点的缓存保留率

typedef unsigned char uint8;
typedef char int8;
typedef int int32;
typedef long long int64;

#define LOG(x) cout << x << endl 
#define ERROR(x) {perror(x); abort();}

#define NOT_ALLOCATION -1 // 不能改为其他值


// min \Sum allocMatrix
// beam search? 