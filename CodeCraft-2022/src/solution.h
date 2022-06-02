#pragma once

#include "CodeCraft-2022.h"
#include <unordered_map>
#include <unordered_set>
#include "heap.hpp"


// 时间片相关的操作
typedef struct Timeslice {
	int32 tid; // 时间片id
	union {
		int32 bandwidth;
		int32 max_load;
		int32 sum_demand;
		int32 score; // universe
	};
	Timeslice(int32 _tid, int32 bd) {
		tid = _tid;
		bandwidth = bd;
	}
	Timeslice() {
		tid = -1;
		bandwidth = 0;
	}
	int32 getkey() {
		return tid;
	}

	bool operator <(const Timeslice& b) const {
		if (sum_demand != b.sum_demand) return sum_demand < b.sum_demand; // 大顶堆
		else return tid < b.tid;
	}

	static bool sort_by_bandwidth_ascend (const Timeslice& a, const Timeslice& b)  {
		if (a.bandwidth != b.bandwidth) return a.bandwidth < b.bandwidth; // 升序. 带宽小的时间片放前面 (评价算分用)
		else return a.tid < b.tid;
	}
	static bool sort_by_bandwidth_descend(const Timeslice& a, const Timeslice& b) {
		if (a.bandwidth != b.bandwidth) return a.bandwidth > b.bandwidth; 
		else return a.tid < b.tid;
	}

	static bool sort_by_max_load_descend(const Timeslice& a, const Timeslice& b) {
		if (a.max_load != b.max_load) return a.max_load > b.max_load; // 最大负载大的时间片放前面.
		else return a.tid < b.tid;
	}
	static bool sort_by_sum_demand_descend(const Timeslice& a, const Timeslice& b) {
		if (a.sum_demand != b.sum_demand) return a.sum_demand > b.sum_demand; // 最大需求大的时间片放前面
		else return a.tid < b.tid;
	}
	static bool sort_by_sum_demand_ascend(const Timeslice& a, const Timeslice& b) {
		if (a.sum_demand != b.sum_demand) return a.sum_demand < b.sum_demand; // 最大需求大的时间片放前面
		else return a.tid < b.tid;
	}
}Timeslice;




struct cmp_min_heap {
	bool operator()(Timeslice& a, Timeslice& b) {
		if (a.bandwidth != b.bandwidth)  return a.bandwidth < b.bandwidth; // 小的在堆顶
		//else return a.tid > b.tid;
		else return false; // 如果值一样直接不变
	}
};
struct cmp_max_heap {
	bool operator()(Timeslice& a, Timeslice& b) {
		if (a.bandwidth != b.bandwidth)  return a.bandwidth > b.bandwidth; // 大的在堆顶
		//else return a.tid > b.tid; 
		else return false; // 如果值一样直接不变
	}
};


// 边缘节点
typedef struct EdgeNode{

	uint8 id;
	uint8 degree;
	int32 max_bandwidth;
	int32 v95;
	int32 v95_tid;
	double score;
	bool operator <(const EdgeNode& b) const {
		if (max_bandwidth != b.max_bandwidth) return max_bandwidth > b.max_bandwidth; // 带宽大优先
		else if (degree != b.degree) return degree > b.degree; // 度大优先
		else return id < b.id;
	}
	static bool sort_by_max_bandwidth_degree(const EdgeNode& a, const EdgeNode& b) {
		if (a.max_bandwidth != b.max_bandwidth) return a.max_bandwidth > b.max_bandwidth;
		else if (a.degree != b.degree) return a.degree < b.degree;
		else return a.id < b.id;
	}
	static bool sort_by_max_bandwidth_ascend(const EdgeNode& a, const EdgeNode& b) {
		if (a.max_bandwidth != b.max_bandwidth) return a.max_bandwidth < b.max_bandwidth;
		else return a.id < b.id;
	}
	static bool sort_by_max_bandwidth_descend(const EdgeNode& a, const EdgeNode& b) {
		if (a.max_bandwidth != b.max_bandwidth) return a.max_bandwidth > b.max_bandwidth;
		else return a.id < b.id;
	}
	static bool sort_by_degree_ascend(const EdgeNode& a, const EdgeNode& b) {
		if (a.degree != b.degree) return a.degree < b.degree;
		else return a.id < b.id;
	}
	static bool sort_by_degree_max_bandwidth_ascend(const EdgeNode& a, const EdgeNode& b) {
		if (a.degree != b.degree) return a.degree < b.degree;
		else if (a.max_bandwidth != b.max_bandwidth) return a.max_bandwidth < b.max_bandwidth;
		else return a.id < b.id;
	}
	static bool sort_by_degree_descend(const EdgeNode& a, const EdgeNode& b) {
		if (a.degree != b.degree) return a.degree > b.degree;
		else return a.id < b.id;
	}
	static bool sort_by_score_ascend(const EdgeNode& a, const EdgeNode& b) {
		if (a.score != b.score) return a.score < b.score; 
		else return a.id < b.id;
	}
	static bool sort_by_score_descend(const EdgeNode& a, const EdgeNode& b) {
		if (a.score != b.score) return a.score > b.score;
		else return a.id < b.id;
	}
}EdgeNode;

// 客户节点
typedef struct ClientNode{
	uint8 id;
	uint8 degree;
	double score;
	bool operator <(const ClientNode& b) const {
		if (degree != b.degree) return degree < b.degree;     // 度小优先
		//else if (demand != b.demand) return demand < b.demand; // 需求量小优先
		else return id < b.id; 
	}
	static bool sort_by_degree_ascend(const ClientNode& a, const ClientNode& b) {
		if (a.degree != b.degree) return a.degree < b.degree;
		else return a.id < b.id;
	}
	static bool sort_by_degree_descend(const ClientNode& a, const ClientNode& b) {
		if (a.degree != b.degree) return a.degree > b.degree;
		else return a.id < b.id;
	}
	static bool sort_by_score_ascend(const ClientNode& a, const ClientNode& b) {
		if (a.score != b.score) return a.score < b.score;
		else return a.id < b.id;
	}
	static bool sort_by_score_descend(const ClientNode& a, const ClientNode& b) {
		if (a.score != b.score) return a.score > b.score;
		else return a.id < b.id;
	}	
}ClientNode;

// 流量需求
typedef struct StreamDemand {
	uint8 stream_id;
	uint8 client_id;
	uint8 client_degree;
	int32 sum_stream;

	int32 demand;

	bool operator <(const StreamDemand& b) const {
		if (demand != b.demand) return demand > b.demand;     // 需求大的优先
		else return stream_id < b.stream_id;
	}

	static bool sort_by_demand_descend(const StreamDemand& a, const StreamDemand& b) {
		int deg_demand_a = a.demand / (a.client_degree);
		int deg_demand_b = b.demand / (b.client_degree);

		if (a.sum_stream != b.sum_stream) return a.sum_stream > b.sum_stream;
		if (a.stream_id != b.stream_id) return a.stream_id < b.stream_id;
		//if (a.stream_id != b.stream_id) return a.client_id < b.client_id;
		//else return a.demand > b.demand;
		else if (a.demand != b.demand) return a.demand > b.demand;
		//else if (deg_demand_a != deg_demand_b) return deg_demand_a > deg_demand_b;
		
		else return a.client_id < b.client_id;
	}

	static bool sort_by_demand_descend2(const StreamDemand& a, const StreamDemand& b) {

		int deg_demand_a = a.demand / a.client_degree;
		int deg_demand_b = b.demand / b.client_degree;

		if (a.demand != b.demand) return a.demand > b.demand;
		//if (deg_demand_a != deg_demand_b) return deg_demand_a > deg_demand_b;
		else if (a.client_id != b.client_id) return a.client_id < b.client_id;
		else return a.stream_id < b.stream_id;
	}


	static bool sort_by_demand_ascend(const StreamDemand& a, const StreamDemand& b) {
		if (a.demand != b.demand) return a.demand < b.demand;
		else if (a.client_id != b.client_id) return a.client_id < b.client_id;
		else return a.stream_id < b.stream_id;
	}
}StreamDemand;

typedef struct SearchState {
	Timeslice edge_usage[MAX_EDGENODE_NUM][MAX_T_LEN]; // 边缘节点每个时间片的带宽使用
	int32 edge_cache[MAX_EDGENODE_NUM][MAX_T_LEN];
	short alloc_mat[MAX_T_LEN][MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES]; // 反过来存
	EdgeNode edge_nodes[MAX_EDGENODE_NUM]; // 边缘节点状态集合
	uint8 demands_alloc_flag[MAX_T_LEN][MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES];
	int32 edge_max_used_bandwidth[MAX_EDGENODE_NUM];
	uint8 use_flag[MAX_EDGENODE_NUM][MAX_T_LEN]; // 标记某个时间片, 某个边缘节点是否已被使用过
	uint8 edge_use_flag[MAX_EDGENODE_NUM]; // 标记某个边缘节点是否被使用过
	int32 select_20_edge[MAX_T_LEN][MAX_EDGENODE_NUM]; // 标记留用率为0.01的机器
}SearchState;


typedef struct Adjancency {
	uint8 nbr_id;
	uint8 nbr_degree;
	static bool sort_by_degree_ascend(const Adjancency& a, const Adjancency& b) {
		if (a.nbr_degree != b.nbr_degree) return a.nbr_degree < b.nbr_degree;
		else return a.nbr_id < b.nbr_id;
	}
	static bool sort_by_degree_decend(const Adjancency& a, const Adjancency& b) {
		if (a.nbr_degree != b.nbr_degree) return a.nbr_degree > b.nbr_degree;
		else return a.nbr_id < b.nbr_id;
	}
}Adjancency;


// 打包比较
typedef struct Package{
	int level;
	double score;
	double score2;
	int select_node;
}Package;

// 每个线程私有的搜索状态变量
typedef struct ThreadLocal {
	EdgeNode edge_nodes[MAX_EDGENODE_NUM];
	EdgeNode sorted_edge_nodes[MAX_EDGENODE_NUM];       // 按优先级排序的边缘节点
	ClientNode sorted_client_nodes[MAX_CLIENTNODE_NUM]; // 按优先级排序的客户节点

	Timeslice edge_usage[MAX_EDGENODE_NUM][MAX_T_LEN];     // 边缘节点带宽使用序列
	Timeslice edge_usage_bak[MAX_EDGENODE_NUM][MAX_T_LEN]; // 备份边缘节点带宽使用序列
	int32 edge_cache[MAX_EDGENODE_NUM][MAX_T_LEN];         // 每个节点的缓存

	uint8 demands_alloc_flag[MAX_T_LEN][MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES];  // 标记客户某个流量类型是否已分配

	bool center_use[MAX_T_LEN]; // 标记是否是中心节点的后5%
	short alloc_mat[MAX_T_LEN][MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES]; // 分配矩阵, 初始化为-1

	int32 temp_each_stream_maxval[MAX_EDGENODE_NUM][MAX_STREAM_TYPES]; // 临时变量: 用于统计一个时刻内每个边缘每个流的最大使用量
	Timeslice temp_center_usage[MAX_T_LEN];                            // 临时变量: 中心节点每个时刻的流量

	int32 each_stream_usage[MAX_T_LEN][MAX_EDGENODE_NUM][MAX_STREAM_TYPES];  // 用于统计所有时刻内每个边缘每个流的总使用量
	uint8 each_stream_count[MAX_T_LEN][MAX_EDGENODE_NUM][MAX_STREAM_TYPES];  // 用于统计所有时刻内每个边缘每个流的总数量
	int32 each_stream_maxval[MAX_T_LEN][MAX_EDGENODE_NUM][MAX_STREAM_TYPES]; // 动态维护所有时刻每个边缘每个流的最大流量 (迁移后无效)
	Timeslice center_usage[MAX_T_LEN];  // 动态维护中心节点每个时刻的流量 (迁移后无效)

	uint8 use_flag[MAX_EDGENODE_NUM][MAX_T_LEN]; // 标记某个时间片, 某个边缘节点是否已被使用过 (或者说后5%中除缓存外的白嫖标记)
	uint8 edge_use_flag[MAX_EDGENODE_NUM];       // 标记某个边缘节点是否被使用过

	Timeslice edge_max_load[MAX_EDGENODE_NUM][MAX_T_LEN];   // 预分配处排序: 每个时刻 每个边缘节点能吸收最大带宽/度
	Timeslice edge_max_load_2[MAX_EDGENODE_NUM][MAX_T_LEN]; // 预分配处排序: 每个时刻 每个边缘节点能吸收最大带宽

	EdgeNode sorted_edges_importance[MAX_EDGENODE_NUM]; // 关机顺序, 会影响关机效果
	uint8 ban_flag[MAX_EDGENODE_NUM];                   // 节点禁用标记
	uint8 ban_flag_temp[MAX_EDGENODE_NUM];

	// 动态维护第k大
	Heap<Timeslice, cmp_min_heap> min_heaps[MAX_EDGENODE_NUM]; // 每个边缘节点维护v95的堆 (迁移后无效)
	Heap<Timeslice, cmp_min_heap> center_min_heap; // 中心节点维护v95的堆 (迁移后无效)
	Heap<Timeslice, cmp_min_heap> migra_min_heap;  // 迁移用的小顶堆
	Heap<Timeslice, cmp_max_heap> migra_max_heap;  // 迁移用的大顶堆 (注意初始化表给定初始化参)

	int32 edge_max_used_bandwidth[MAX_EDGENODE_NUM];  // 所有时刻中的最大使用流量

	int32 center_v95; // 中心节点的v95
	int32 center_v95_t; // 中心节点v95的时刻t

	SearchState state;  // 搜索状态备份
	SearchState state2; // 搜索状态备份 (用于迁移过程2)

	int32 cost;    // 总分
	int32 ban_cnt; // 关机数量

	
	uint8 tag[MAX_T_LEN]; // 预处理中标记是否已使用的临时数组
	int32 select_20_edge[MAX_T_LEN][MAX_EDGENODE_NUM]; // 每个时刻, 标记节点是否被选中 (选中后缓存保留率为0.01)

	ThreadLocal() : migra_max_heap(MAX_T_LEN) {

		memset(tag, 0, sizeof(tag));
		memset(alloc_mat, NOT_ALLOCATION, sizeof(alloc_mat));
		memset(demands_alloc_flag, 0, sizeof(demands_alloc_flag));
		memset(edge_usage, 0, sizeof(edge_usage));
		memset(edge_usage_bak, 0, sizeof(edge_usage_bak));
		memset(edge_cache, 0, sizeof(edge_cache));
		memset(edge_max_used_bandwidth, 0, sizeof(edge_max_used_bandwidth));
		memset(edge_max_load, 0, sizeof(edge_max_load));
		memset(edge_max_load_2, 0, sizeof(edge_max_load_2));
		memset(use_flag, 0, sizeof(use_flag));
		memset(edge_use_flag, 0, sizeof(edge_use_flag));
		memset(ban_flag, 0, sizeof(ban_flag));
		memset(temp_each_stream_maxval, 0, sizeof(temp_each_stream_maxval));
		memset(each_stream_usage, 0, sizeof(each_stream_usage));
		memset(each_stream_count, 0, sizeof(each_stream_count));
		memset(temp_center_usage, 0, sizeof(temp_center_usage));
		memset(each_stream_maxval, 0, sizeof(each_stream_maxval));
		memset(center_usage, 0, sizeof(center_usage));
		memset(center_use, 0, sizeof(center_use));
		memset(select_20_edge, 0, sizeof(select_20_edge));
		center_v95 = 0, center_v95_t = -1;
	}
}ThreadLocal;


class Solution {

public:
	ThreadLocal *locals[CPU_CORES];
	int32 best_local_id;

	// 全局参数
	int32 qos_limit; // 全局流量限制
	int32 base_cost; // 基本费用
	int32 T;         // 时间片数量
	int32 max_free_usage_count; // 最大免费使用次数
	int32 max_free_usage_90_count;
	int32 g95_index; // 95分位索引值
	int32 g90_index;
	int32 edge_node_count;   // 边缘节点数量
	int32 client_node_count; // 客户节点数量
	double center_cost;      // 中心节点系数
	int32 conn_edge_count;   // 可联通的边缘节点数量
	int32 max_max_bandwidth; // 全局最大带宽
	int32 max_degree;        // 全局最大度

	Adjancency client_adj[MAX_CLIENTNODE_NUM][MAX_EDGENODE_NUM]; // 客户节点的邻居 
	Adjancency edge_adj[MAX_EDGENODE_NUM][MAX_CLIENTNODE_NUM];   // 边缘节点的邻居 

	// 用于输出映射
	unordered_map<int32, string> edge_id_str;
	unordered_map<int32, string> client_id_str;
	unordered_map<string, int32> edge_str_id;
	unordered_map<string, int32> client_str_id;

	int32 qos_mat[MAX_CLIENTNODE_NUM][MAX_EDGENODE_NUM]; // 客户与边缘节点的连接质量 

	// 边缘节点
	EdgeNode edge_nodes[MAX_EDGENODE_NUM];
	// 客户节点
	ClientNode client_nodes[MAX_CLIENTNODE_NUM];
	// 流量相关
	int32 stream_types_num[MAX_T_LEN];                  // 每个时刻的类型数
	int32 demands_num[MAX_T_LEN];                       // 每个时刻需求数
	vector<vector<string>> stream_types_name;           // 每个时刻类型名

	StreamDemand(*demands_all_types)[MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES];           // 未排序的所有需求
	StreamDemand(*sorted_demands_type_desc)[MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES];   // 需求先按流类型排(不同类型按流总量大小). 然后同类型从大到小排
	StreamDemand(*sorted_demands_desc)[MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES];        // 需求从大到小排

	Timeslice sum_demands[MAX_T_LEN]; // 每个时间片的总带宽
	Timeslice client_sum_demand[MAX_CLIENTNODE_NUM][MAX_T_LEN]; // 客户总流量

private:

	// 线程计数器
	atomic<int> thread_counter;

	// 分配与迁移的函数指针数组 在init中初始化.
	bool (Solution::*ptr_alloc_func[MAX_CPU_CORES])(int strategy_id); // 函数指针数组
	int (Solution::*ptr_migra_func[MAX_CPU_CORES])();
	int strategy_num[MAX_CPU_CORES];

	void (Solution::* ptr_pre_alloc_func[MAX_FUN_PTR_NUM])();      // 预分配策略列表
	bool (Solution::* ptr_remain_alloc_func[MAX_FUN_PTR_NUM])();   // 剩余需求分配策略列表
	int pre_alloc_func_num, remain_alloc_func_num;

	// 分配线程私有id
	inline int get_local_id() {
		static thread_local int id = -1;
		if (id == -1) id = assign_local_id();
		if (best_local_id != -1) return best_local_id; // 如果已经确定最优id. 直接返回最优id..
		return id;
	}
	inline int assign_local_id() {
		return (thread_counter++) % CPU_CORES;
	}

	bool alloc_func_list_1(int strategy_id);
	bool alloc_func_list_2(int strategy_id);
	bool alloc_func_list_3(int strategy_id);
	bool alloc_func_list_4(int strategy_id);

	void mount_func(); // 挂载分配策略到函数指针上
	void arrange_tasks(); // 编排不同线程的任务



	// 在 init函数的末尾:
	// 1. 修改 alloc_func_list_1 或 migration_list
	// 2. 设置好每个list中的策略数 strategy_num


	int alloc(int e_node, int c_node, int bandwidth, int time, bool update_95, int remain); // 将客户节点的x流量放置到边缘节点
	void move(int from_e, int to_e, int c_node, int bandwidth, int time, bool update_to_95); // 将客户节点的x流量从A边缘移动到B边缘

	void backup(SearchState& state);
	void recovery(SearchState& state);

	void reset();

	int canuse_bandwidth(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN],int flag = 0);
	int freeuse_bandwidth(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN], bool infer_all = false);
	void update_future_cache(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN], bool update_heap, bool update_all = false);

	Timeslice get_v95(int eid);

	void pre_alloc_strategy_001();
	void pre_alloc_strategy_002();
	bool alloc_strategy_002();
	bool alloc_strategy_001();
	void migration_edge();
	void migration_center();
	int migration_process();
	int migration_process_2();

	void select_today_20(int t);
	void recal_cache(int t);

	int resort_cal_cost();
	double cal_cost(int edge_node);

	void cal_edge_importance(int best_strategy);
	void gen_solution_mthreads();
	void rescue();
	void verify(); // 验证解是否合法
	void reduce();


public:
	void init();
	void solve();
	void evaluate();

	Solution() {
		thread_counter = 0;

		max_free_usage_count = 0, g95_index = 0;
		edge_node_count = client_node_count = 0;
		qos_limit = T = 0;

		demands_all_types = new StreamDemand[MAX_T_LEN][MAX_CLIENTNODE_NUM][MAX_STREAM_TYPES];
		sorted_demands_type_desc = new StreamDemand[MAX_T_LEN][MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES];
		sorted_demands_desc = new StreamDemand[MAX_T_LEN][MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES];

		for (int i = 0; i < CPU_CORES; ++i) locals[i] = new ThreadLocal();

		memset(demands_all_types, 0, sizeof(StreamDemand) * MAX_T_LEN * MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES);
		memset(sorted_demands_type_desc, 0, sizeof(StreamDemand) * MAX_T_LEN * MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES);
		memset(sorted_demands_desc, 0, sizeof(StreamDemand) * MAX_T_LEN * MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES);
		memset(edge_nodes, 0, sizeof(edge_nodes));
		memset(qos_mat, 0, sizeof(qos_mat));
		memset(demands_num, 0, sizeof(demands_num));
		memset(sum_demands, 0, sizeof(sum_demands));
		memset(client_sum_demand, 0, sizeof(client_sum_demand));
		memset(client_adj, 0, sizeof(client_adj));
		memset(edge_adj, 0, sizeof(edge_adj));
	}

	~Solution() {

		delete[] demands_all_types;
		delete[] sorted_demands_type_desc;
		delete[] sorted_demands_desc;

		for (int i = 0; i < CPU_CORES; ++i) delete locals[i];
	}

	
};