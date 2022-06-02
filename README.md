# Huawei-CodeCraft-2022

2022华为软件精英挑战赛 - 杭厦赛区 - 土豪法称霸杭厦 - 决赛季军



### 运行

测试环境: ubuntu20.04; gcc 7.5; cmake 3.16.

进入目录执行

```
./build_and_run.sh
```

判题器集成在算法中（验证解的合法性和算分）

官方数据成本在9.7w左右



### 框架

- 通过arrange_tasks编排每个线程的搜索任务

- 多线程启动任务gen_solution_mthreads搜索解

- 通过reduce聚合所有结果并选择最优解

- 调用rescue矫正缓存，修正不合法的解

  

### 配置

主要配置都在CodeCraft-2022.h中：

```c++
#define MAX_CLIENTNODE_NUM 35 // 最大客户节点数量
#define MAX_EDGENODE_NUM 135  // 最大边缘节点数量
#define MAX_T_LEN 8950        // 最大时间序列长度
#define MAX_STREAM_TYPES 100  // 最大流数量
#define CPU_CORES 1           // 设置的CPU核心数
#define MAX_CPU_CORES 4       // 支持的最大CPU核心数 (注意内存)
#define TRY_SHUT_DOWN 1       // 是否尝试关机以搜索更优解
#define MAX_FUN_PTR_NUM 100   // 函数指针数组最大长度
#define GLOBAL_MAX_TIME_MS (275 * 1000) // 运行时间限制(ms)
```



### 其他

- 一开始大多数变量直接预先开了数组为了贪图一点效率，实际上快不了多少并且使得代码很丑见谅orz。
- 有问题欢迎issue讨论~

