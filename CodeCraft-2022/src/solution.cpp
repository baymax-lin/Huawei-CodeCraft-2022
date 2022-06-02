#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include "solution.h"
#include "heap.hpp"
#include "reader.h"
#include "CodeCraft-2022.h"

using namespace std;


/*
  决赛新增: 选择20台的缓存保留率为0.01
*/

void Solution::select_today_20(int t) {
	if (t == 0) return;
	ThreadLocal* local = locals[get_local_id()];

	vector<pair<int, int>> temp;
	for (int i = 0; i < edge_node_count; ++i) {
		int usage = local->edge_usage[i][t - 1].bandwidth + local->edge_cache[i][t - 1];
		temp.push_back(make_pair(i, usage));
	}

	sort(temp.begin(), temp.end(), [](const pair<int, int>& p1, const pair<int, int>& p2) {
		if (p1.second != p2.second) return p1.second > p2.second;
		else return p1.first < p2.first;
		});

	int canuse = 20;

	memset(local->select_20_edge[t], 0, sizeof(local->select_20_edge[0]));
	for (int i = 0; i < edge_node_count; ++i) {
		if (local->edge_nodes[i].degree == 0) continue;
		if (local->use_flag[i][t - 1] == 1) {
			canuse--;
			local->select_20_edge[t][i] = 1;
		}
		if (canuse < 0) {
			ERROR("too much selected nodes in the same time");
		}
	}


	int num = min(canuse, edge_node_count);
	for (int i = 0; i < temp.size(); ++i) {
		int edge = temp[i].first;
		if (local->select_20_edge[t][edge] == 1) continue;
		local->select_20_edge[t][edge] = 1;
		canuse--;
		if (canuse == 0) break;
	}

	// 选择的20个节点必须cover use_flag中白嫖占用的节点

}

/*
* 选择出保留率为0.01的边缘节点后
* 需要更新未来几个时刻的缓存值
*/
void Solution::recal_cache(int t) {
	ThreadLocal* local = locals[get_local_id()];
	for (int i = 0; i < edge_node_count; ++i) {
		update_future_cache(i, t, local->edge_usage, local->edge_cache, false);
	}
}


/*
	将不同预分配/分配挂载到函数指针数组上
*/

void Solution::mount_func() {

	ptr_pre_alloc_func[0] = &Solution::pre_alloc_strategy_001;
	//ptr_pre_alloc_func[1] = &Solution::pre_alloc_strategy_002;

	pre_alloc_func_num = 2;

	ptr_remain_alloc_func[0] = &Solution::alloc_strategy_001;
	ptr_remain_alloc_func[1] = &Solution::alloc_strategy_002;

	remain_alloc_func_num = 2;
}

/*
  本文件放置计算任务编排
  多核任务分配
  生成解的过程会调用一次 ptr_alloc_func 生成合法解
  而后调用一次 ptr_migra_func 进行迁移降低成本

  ============== tasks ==============
*/

void Solution::arrange_tasks() {

	// 设置多线程搜索任务 
	ptr_alloc_func[0] = &Solution::alloc_func_list_1;
	ptr_alloc_func[1] = &Solution::alloc_func_list_2;
	ptr_alloc_func[2] = &Solution::alloc_func_list_3;
	ptr_alloc_func[3] = &Solution::alloc_func_list_4;

	ptr_migra_func[0] = &Solution::migration_process;
	ptr_migra_func[1] = &Solution::migration_process;
	ptr_migra_func[2] = &Solution::migration_process;
	ptr_migra_func[3] = &Solution::migration_process;

	strategy_num[0] = 2;
	strategy_num[1] = 2;
	strategy_num[2] = 2;
	strategy_num[3] = 2;

}

/*
* alloc_func_list1-4 对应cpu1-4的调用的分配过程
* strategy_id表示选择该分配过程下的第几种预分配/分配策略
* 
*/
bool Solution::alloc_func_list_1(int strategy_id) {

	ThreadLocal* local = locals[get_local_id()];

	bool res = true;

	for (int i = 0; i < edge_node_count; ++i) {
		local->sorted_edge_nodes[i].score = local->sorted_edge_nodes[i].max_bandwidth * 1.0 * (local->sorted_edge_nodes[i].degree + 1e-2);
	}
	sort(local->sorted_edge_nodes, local->sorted_edge_nodes + edge_node_count, EdgeNode::sort_by_score_descend);

	int pre_func_id = strategy_id / 2;
	int alloc_fun_id = strategy_id % 2;
	std::invoke(ptr_pre_alloc_func[pre_func_id], this);
	res = std::invoke(ptr_remain_alloc_func[alloc_fun_id], this); 

	return res;
}

bool Solution::alloc_func_list_2(int strategy_id) {

	ThreadLocal* local = locals[get_local_id()];

	bool res = true;

	for (int i = 0; i < edge_node_count; ++i) {
		local->sorted_edge_nodes[i].score = local->sorted_edge_nodes[i].max_bandwidth * 1.0 * sqrt(sqrt(local->sorted_edge_nodes[i].degree + 1e-2));
	}
	sort(local->sorted_edge_nodes, local->sorted_edge_nodes + edge_node_count, EdgeNode::sort_by_score_descend);
	
	int pre_func_id = strategy_id / 2;
	int alloc_fun_id = strategy_id % 2;

	std::invoke(ptr_pre_alloc_func[pre_func_id], this);
	res = std::invoke(ptr_remain_alloc_func[alloc_fun_id], this); 

	return res;
}

bool Solution::alloc_func_list_3(int strategy_id) {

	ThreadLocal* local = locals[get_local_id()];

	bool res = true;

	for (int i = 0; i < edge_node_count; ++i) {
		local->sorted_edge_nodes[i].score = local->sorted_edge_nodes[i].degree;
	}
	sort(local->sorted_edge_nodes, local->sorted_edge_nodes + edge_node_count, EdgeNode::sort_by_score_descend);

	int pre_func_id = strategy_id / 2;
	int alloc_fun_id = strategy_id % 2;

	std::invoke(ptr_pre_alloc_func[pre_func_id], this);
	res = std::invoke(ptr_remain_alloc_func[alloc_fun_id], this); 

	return res;
}


bool Solution::alloc_func_list_4(int strategy_id) {

	ThreadLocal* local = locals[get_local_id()];

	bool res = true;

	for (int i = 0; i < edge_node_count; ++i) {
		local->sorted_edge_nodes[i].score = sqrt(sqrt(local->sorted_edge_nodes[i].max_bandwidth)) * 1.0 * local->sorted_edge_nodes[i].degree;
	}
	sort(local->sorted_edge_nodes, local->sorted_edge_nodes + edge_node_count, EdgeNode::sort_by_score_descend);

	int pre_func_id = strategy_id / 2;
	int alloc_fun_id = strategy_id % 2;

	std::invoke(ptr_pre_alloc_func[pre_func_id], this);
	res = std::invoke(ptr_remain_alloc_func[alloc_fun_id], this);

	return res;
}


/*
  评价成本
*/
void Solution::evaluate() {

	int64 score_int = 0;
	double score = 0;

	int degree_0_count = 0;

	if (best_local_id == -1) {
		ERROR("best_local_id == -1");
	}

	ThreadLocal* local = locals[best_local_id];

	for (int i = 0; i < edge_node_count; ++i) {
		for (int t = 0; t < T; ++t) {
			local->edge_usage_bak[i][t].tid = t;
			local->edge_usage_bak[i][t].bandwidth = local->edge_usage[i][t].bandwidth + local->edge_cache[i][t];
		}
	}

	int v95sum = 0;
	for (int32 i = 0; i < edge_node_count; ++i) {
		sort(local->edge_usage_bak[i], local->edge_usage_bak[i] + T, Timeslice::sort_by_bandwidth_ascend);

		int select_t_count = max_free_usage_count;

		if (local->edge_nodes[i].degree == 0) degree_0_count++;

		int index = g95_index;

		local->edge_nodes[i].v95_tid = local->edge_usage_bak[i][index].tid;

		if (local->edge_usage_bak[i][T - 1].bandwidth > 0) {
			int32 W = local->edge_usage_bak[i][index].bandwidth;
			v95sum += W;
			if (W >= 0 && W <= base_cost) {
				score += base_cost;
			}
			else {
				score += (W - base_cost) * 1.0 * (W - base_cost) / local->edge_nodes[i].max_bandwidth + W;
			}
		}
	}
	int edge_score = int64(score + 0.5);

	// 计算中心节点代价
	// 首先从分配矩阵中统计每天中心节点使用量
	for (int t = 0; t < T; ++t) {
		memset(local->temp_each_stream_maxval, 0, sizeof(local->temp_each_stream_maxval));

		for (int j = 0; j < client_node_count; ++j) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				int i = local->alloc_mat[t][j][s];
				if (i != NOT_ALLOCATION) {
					int batch_demand = demands_all_types[t][j][s].demand;
					local->temp_each_stream_maxval[i][s] = max(local->temp_each_stream_maxval[i][s], batch_demand);
				}
			}
		}
		int32 sum_bd = 0;

		for (int i = 0; i < edge_node_count; ++i) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				sum_bd += local->temp_each_stream_maxval[i][s];
			}
		}
		local->temp_center_usage[t].tid = t;
		local->temp_center_usage[t].bandwidth = sum_bd;
	}

	sort(local->temp_center_usage, local->temp_center_usage + T, Timeslice::sort_by_bandwidth_ascend);
	int index = g95_index;

	int32 W = local->temp_center_usage[index].bandwidth;
	score += W * center_cost;

	score_int = int64(score + 0.5);

	LOG("edge v95 sum: " << v95sum);
	LOG("edge score = " << edge_score);
	LOG("center score = " << score_int - edge_score);
	LOG("total score = " << score_int);

	return;

}

/*
  重置搜索状态
*/
void Solution::reset() {

	ThreadLocal* local = locals[get_local_id()];

	memset(local->edge_max_used_bandwidth, 0, sizeof(local->edge_max_used_bandwidth));
	memset(local->alloc_mat, NOT_ALLOCATION, sizeof(local->alloc_mat));
	memset(local->demands_alloc_flag, 0, sizeof(local->demands_alloc_flag));
	memset(local->edge_usage, 0, sizeof(local->edge_usage));

	for (int cc = 0; cc < edge_node_count; cc++) {
		for (int tt = 0; tt < T; tt++) {
			local->edge_usage[cc][tt].tid = tt;
		}
	}

	memset(local->use_flag, 0, sizeof(local->use_flag));
	memset(local->edge_use_flag, 0, sizeof(local->edge_use_flag));
	memset(local->edge_cache, 0, sizeof(local->edge_cache));
	memset(local->each_stream_maxval, 0, sizeof(local->each_stream_maxval));
	memset(local->each_stream_usage, 0, sizeof(local->each_stream_usage));
	memset(local->each_stream_count, 0, sizeof(local->each_stream_count));
	memset(local->temp_each_stream_maxval, 0, sizeof(local->temp_each_stream_maxval));
	memset(local->center_usage, 0, sizeof(local->center_usage));
	memset(local->select_20_edge, 0, sizeof(local->select_20_edge));

	for (int i = 0; i < edge_node_count; ++i) {
		local->edge_nodes[i].v95 = 0;
		local->edge_nodes[i].v95_tid = -1;
		local->min_heaps[i].reset();
	}

	local->center_min_heap.reset();

}


void Solution::init() {

	// 初始化最优local为-1
	best_local_id = -1;

	// 计算客户节点到边缘节点的度
	for (int j = 0; j < client_node_count; ++j) {
		int conn = 0;
		for (int i = 0; i < edge_node_count; ++i) {
			if (qos_mat[j][i] < qos_limit) {
				conn++;
			}
		}
		client_nodes[j].id = j;
		client_nodes[j].degree = conn;
	}

	max_degree = 0;
	max_max_bandwidth = 0;
	// 1. 计算边缘节点的度
	int disconn = 0;
	conn_edge_count = 0;
	for (int i = 0; i < edge_node_count; ++i) {
		int conn = 0;
		for (int j = 0; j < client_node_count; ++j) {
			if (qos_mat[j][i] < qos_limit) {
				conn++;
			}
		}
		if (conn == 0) {
			disconn++;
		}
		else conn_edge_count++;

		int score = edge_nodes[i].max_bandwidth * sqrt(sqrt(conn)); 

		edge_nodes[i].id = i;
		edge_nodes[i].degree = conn;
		edge_nodes[i].score = score;
		max_degree = max(conn, max_degree);
		max_max_bandwidth = max(edge_nodes[i].max_bandwidth, max_max_bandwidth);

	}
	LOG("disconn nodes: " << disconn);

	// 构建邻接表 (客户)
	for (int j = 0; j < client_node_count; ++j) {
		int count = 0;
		for (int i = 0; i < edge_node_count; ++i) {
			if (qos_mat[j][i] < qos_limit) {
				client_adj[j][count].nbr_id = i;
				client_adj[j][count].nbr_degree = edge_nodes[i].degree;
				count++;
			}
		}
		sort(client_adj[j], client_adj[j] + count, Adjancency::sort_by_degree_ascend);
	}
	// 构建邻接表 (边缘)
	for (int i = 0; i < edge_node_count; ++i) {
		int count = 0;
		for (int j = 0; j < client_node_count; ++j) {
			if (qos_mat[j][i] < qos_limit) {
				edge_adj[i][count].nbr_id = j;
				edge_adj[i][count].nbr_degree = client_nodes[j].degree;
				count++;
			}
		}
		sort(edge_adj[i], edge_adj[i] + count, Adjancency::sort_by_degree_ascend);
	}

	g95_index = int32(ceil(T * 0.95)) - 1;
	max_free_usage_count = T - g95_index - 1; // 非核心节点的最大使用次数 T = 100 , 就只能用5次

	g90_index = int32(ceil(T * 0.90 + 0.5)) - 1;
	max_free_usage_90_count = T - g90_index - 1;


	// 抽象需求:
	for (int t = 0; t < T; ++t) {
		int num = 0;
		int sum_demand = 0;

		for (int s = 0; s < stream_types_num[t]; ++s) {
			int sum_demand_stream = 0;
			for (int j = 0; j < client_node_count; ++j) {
				sum_demand_stream += demands_all_types[t][j][s].demand;
			}
			for (int j = 0; j < client_node_count; ++j) {
				sorted_demands_type_desc[t][num].client_id = j;
				sorted_demands_type_desc[t][num].stream_id = s;
				sorted_demands_type_desc[t][num].client_degree = client_nodes[j].degree;
				sorted_demands_type_desc[t][num].demand = demands_all_types[t][j][s].demand;
				sorted_demands_type_desc[t][num].sum_stream = sum_demand_stream;
				++num;
				sum_demand += demands_all_types[t][j][s].demand;

				// 更新客户某时刻的总容量
				client_sum_demand[j][t].tid = t;
				client_sum_demand[j][t].sum_demand += demands_all_types[t][j][s].demand;
			}
		}
		sum_demands[t].sum_demand = sum_demand; // 记录总需求量
		sum_demands[t].tid = t;

		demands_num[t] = num;
		sort(sorted_demands_type_desc[t], sorted_demands_type_desc[t] + num, StreamDemand::sort_by_demand_descend);
	}

	memcpy(sorted_demands_desc, sorted_demands_type_desc, sizeof(StreamDemand) * MAX_T_LEN * MAX_CLIENTNODE_NUM * MAX_STREAM_TYPES);
	//memcpy(time_demands, sum_demands, sizeof(sum_demands));
	for (int t = 0; t < T; ++t) {
		sort(sorted_demands_desc[t], sorted_demands_desc[t] + demands_num[t], StreamDemand::sort_by_demand_descend2);
	}

	sort(sum_demands, sum_demands + T, Timeslice::sort_by_sum_demand_ascend); // 每个时刻总需求升序(缓慢扩容, 尽量填满)

	// 线程局部 locals 信息初始化
	for (int c = 0; c < CPU_CORES; ++c) {
		ThreadLocal* local = locals[c];

		//1. 初始化各个分数为-1
		local->cost = -1;
		local->ban_cnt = 0;

		for (int i = 0; i < edge_node_count; ++i) {
			local->min_heaps[i].set_K(max_free_usage_count + 1); // 小顶堆维护top k个最大 (k >= 2)
		}
		local->center_min_heap.set_K(max_free_usage_count + 1); // 中心节点的小顶堆

		// edge_nodes的信息也得拷贝到local中 (v95...等等一系列信息的更新)
		memcpy(local->edge_nodes, edge_nodes, sizeof(edge_nodes));

		memcpy(local->sorted_edge_nodes, edge_nodes, sizeof(edge_nodes));
		memcpy(local->sorted_client_nodes, client_nodes, sizeof(client_nodes));

		// 防止无流量时 排序错误
		for (int cc = 0; cc < edge_node_count; cc++) {
			for (int tt = 0; tt < T; tt++) {
				local->edge_usage[cc][tt].tid = tt;
			}
		}
	}

	mount_func();
	arrange_tasks();



}

/*
* 某个时刻, 将需求分配到一个边缘节点
* 更新对应的中心节点代价
* 更新对应的v95(堆)
* 更新缓存
*/
int Solution::alloc(int e_node, int c_node, int stream_id, int time, bool update_95, int remain = -1) {

	ThreadLocal* local = locals[get_local_id()];

	int free = 0;
	if (remain == -1) {
		free = canuse_bandwidth(e_node, time, local->edge_usage, local->edge_cache);
	}
	else free = remain;

	int prev_free = free;
	int bd = demands_all_types[time][c_node][stream_id].demand;
	local->edge_usage[e_node][time].tid = time;


	if (bd < 0) {
		return 0;
		ERROR("bd < 0");
	}

	if (free >= bd) {
		local->edge_usage[e_node][time].bandwidth += bd;
		local->alloc_mat[time][c_node][stream_id] = e_node;
		local->demands_alloc_flag[time][c_node][stream_id] = 1;
		free -= bd;
	}
	else {
		ERROR("free < bd !");
	}

	local->edge_max_used_bandwidth[e_node] = max(local->edge_max_used_bandwidth[e_node], 
		local->edge_usage[e_node][time].bandwidth + local->edge_cache[e_node][time]);

	Timeslice use = { time, local->edge_usage[e_node][time].bandwidth + local->edge_cache[e_node][time] };
	local->min_heaps[e_node].push_topK(use);

	update_future_cache(e_node, time, local->edge_usage, local->edge_cache, true);

	// 动态更新目标节点对该类流量的最大带宽需求量
	if (local->each_stream_maxval[time][e_node][stream_id] < bd) {
		int diff = bd - local->each_stream_maxval[time][e_node][stream_id];
		local->center_usage[time].tid = time;
		local->center_usage[time].bandwidth += diff;
		local->each_stream_maxval[time][e_node][stream_id] = bd;
		local->center_min_heap.push_topK(local->center_usage[time]);
	}

	return prev_free - free;
}

/*
* 将某个需求从边缘节点A迁移到边缘节点B
* 更新边缘节点B中维护v95的堆 
* 边缘节点A的堆不再维护, 通过另外的方法算出
* 更新缓存
*/
void Solution::move(int from_e, int to_e, int c_node, int stream_id, int time, bool update_to_95) {

	ThreadLocal* local = locals[get_local_id()];

	int to_free = canuse_bandwidth(to_e, time, local->edge_usage, local->edge_cache);
	int to_prev_free = to_free;
	int from_free = canuse_bandwidth(from_e, time, local->edge_usage, local->edge_cache);
	int from_prev_free = from_free;
	int bd = demands_all_types[time][c_node][stream_id].demand;


	if (bd < 0) ERROR("bd < 0");

	if (to_free >= bd) {
		local->alloc_mat[time][c_node][stream_id] = to_e;

		local->edge_usage[from_e][time].bandwidth -= bd;
		local->edge_usage[to_e][time].bandwidth += bd;

		from_free += bd;

	}
	else {
		return;
		ERROR("to_free < bd!");
	}


	local->edge_max_used_bandwidth[to_e] = max(local->edge_max_used_bandwidth[to_e],
		local->edge_usage[to_e][time].bandwidth + local->edge_cache[to_e][time]);

	update_future_cache(to_e, time, local->edge_usage, local->edge_cache, false);
	update_future_cache(from_e, time, local->edge_usage, local->edge_cache, false);

	// 更新目标节点对该类流量的最大带宽需求量
	if (local->each_stream_maxval[to_e][time][stream_id] < bd) {
		int diff = bd - local->each_stream_maxval[to_e][time][stream_id];
		local->center_usage[time].tid = time;
		local->center_usage[time].bandwidth += diff;
		local->each_stream_maxval[to_e][time][stream_id] = bd;
	}


}


/*
* 备份搜索状态
*/
void Solution::backup(SearchState& state) {
	ThreadLocal* local = locals[get_local_id()];

	memcpy(state.alloc_mat, local->alloc_mat, sizeof(local->alloc_mat));
	memcpy(state.edge_usage, local->edge_usage, sizeof(local->edge_usage));
	memcpy(state.edge_nodes, local->edge_nodes, sizeof(local->edge_nodes));
	memcpy(state.demands_alloc_flag, local->demands_alloc_flag, sizeof(local->demands_alloc_flag));
	memcpy(state.edge_max_used_bandwidth, local->edge_max_used_bandwidth, sizeof(local->edge_max_used_bandwidth));
	memcpy(state.use_flag, local->use_flag, sizeof(local->use_flag));
	memcpy(state.edge_use_flag, local->edge_use_flag, sizeof(local->edge_use_flag));
	memcpy(state.edge_cache, local->edge_cache, sizeof(local->edge_cache));
	memcpy(state.select_20_edge, local->select_20_edge, sizeof(local->select_20_edge));
}

/*
* 恢复备份区的搜索状态到主内存
*/
void Solution::recovery(SearchState & state) {
	ThreadLocal* local = locals[get_local_id()];

	memcpy(local->alloc_mat, state.alloc_mat, sizeof(local->alloc_mat));
	memcpy(local->edge_usage, state.edge_usage, sizeof(local->edge_usage));
	memcpy(local->edge_nodes, state.edge_nodes, sizeof(local->edge_nodes));
	memcpy(local->demands_alloc_flag, state.demands_alloc_flag, sizeof(local->demands_alloc_flag));
	memcpy(local->edge_max_used_bandwidth, state.edge_max_used_bandwidth, sizeof(local->edge_max_used_bandwidth));
	memcpy(local->use_flag, state.use_flag, sizeof(local->use_flag));
	memcpy(local->edge_use_flag, state.edge_use_flag, sizeof(local->edge_use_flag));
	memcpy(local->edge_cache, state.edge_cache, sizeof(local->edge_cache));
	memcpy(local->select_20_edge, state.select_20_edge, sizeof(local->select_20_edge));

}

/*
* 重新计算代价
*/
int Solution::resort_cal_cost() {

	ThreadLocal* local = locals[get_local_id()];

	// 合并cache和使用量
	for (int i = 0; i < edge_node_count; ++i) {
		for (int t = 0; t < T; ++t) {
			local->edge_usage_bak[i][t].tid = t;
			local->edge_usage_bak[i][t].bandwidth = local->edge_usage[i][t].bandwidth + local->edge_cache[i][t];
		}
	}

	double score = 0;
	int score_int = 0;
	for (int32 i = 0; i < edge_node_count; ++i) {
		sort(local->edge_usage_bak[i], local->edge_usage_bak[i] + T, Timeslice::sort_by_bandwidth_ascend);

		int index = g95_index;

		if (local->edge_usage_bak[i][T - 1].bandwidth > 0) {
			int32 W = local->edge_usage_bak[i][index].bandwidth;
			if (W >= 0 && W <= base_cost) {
				score += base_cost;
			}
			else {
				score += (W - base_cost) * 1.0 * (W - base_cost) / local->edge_nodes[i].max_bandwidth + W;
			}
		}
	}

	// 计算中心节点代价
	// 首先从分配矩阵中统计每天中心节点使用量
	for (int t = 0; t < T; ++t) {
		memset(local->temp_each_stream_maxval, 0, sizeof(local->temp_each_stream_maxval));

		for (int j = 0; j < client_node_count; ++j) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				int i = local->alloc_mat[t][j][s];
				if (i != NOT_ALLOCATION) {
					int batch_demand = demands_all_types[t][j][s].demand;
					local->temp_each_stream_maxval[i][s] = max(local->temp_each_stream_maxval[i][s], batch_demand);
				}
			}
		}
		int32 sum_bd = 0;
		for (int i = 0; i < edge_node_count; ++i) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				sum_bd += local->temp_each_stream_maxval[i][s];
			}
		}
		local->temp_center_usage[t].tid = t;
		local->temp_center_usage[t].bandwidth = sum_bd;

	}

	sort(local->temp_center_usage, local->temp_center_usage + T, Timeslice::sort_by_bandwidth_ascend);
	int index = g95_index;
	int32 W = local->temp_center_usage[index].bandwidth;
	score += W * center_cost;



	score_int = int64(score + 0.5);
	return score_int;
}

/*
* 计算某个边缘节点的成本
*/
double Solution::cal_cost(int edge_node) {

	ThreadLocal* local = locals[get_local_id()];

	double cost = 0;
	int max_use = local->edge_max_used_bandwidth[edge_node];
	int i = edge_node;
	if (max_use > 0) {
		int32 W = get_v95(edge_node).bandwidth;
		if (W >= 0 && W <= base_cost) {
			cost += base_cost;
		}
		else {
			cost += (W - base_cost) * 1.0 * (W - base_cost) / local->edge_nodes[i].max_bandwidth + W;
		}
	}
	return cost;
}

/*
* 仅在迁移时使用(分配阶段的v95应从堆取)
* 仅考虑未来4个时刻的约束以提高效率
*/
int Solution::freeuse_bandwidth(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN], bool infer_all) {

	ThreadLocal* local = locals[get_local_id()];

	int max_bd = local->edge_nodes[eid].max_bandwidth;
	int x0 = 0, x1 = 0;
	int ret2 = 0;

	// 对于未来. 要多预测几天 100w最多影响 (~4.6 = 4-5天) 从第四天开始回溯条件
	// (Ct + Ut + Xt) <= v95 or maxbd
	// A1 = (Ct + Ut + Xt) * 0.05 + Ut+1 <= v95 or maxbd
	// A2 = A1 * 0.05 + Ut+2 <= v95 or maxbd
	// A3 = A2 * 0.05 + Ut+3 <= v95 or maxbd
	// A4 = A3 * 0.05 + Ut+4 <= v95 or maxbd
	// 从A4开始回溯

	double cache_p = CACHE_P1;
	int v95 = max(base_cost, local->edge_nodes[eid].v95);
	v95 = min(v95, max_bd);

	if (time == T - 1 ) {
		x0 = local->use_flag[eid][time] ? max_bd : v95;
		x0 = x0 - cache[eid][time] - usage[eid][time].bandwidth;
		ret2 = x0;
	}
	else {

		int RT = min(time + 4, T - 1);
		if (infer_all) RT = T - 1;

		int r_constrait = local->use_flag[eid][RT] == 0 ? v95 : max_bd;
		for (int t = RT; t > time; --t) {
			if (local->select_20_edge[t][eid] == 1)
				cache_p = CACHE_P2;
			else
				cache_p = CACHE_P1;
			int r = (r_constrait + 1 - usage[eid][t].bandwidth) * 1.0 / max(cache_p, 1e-3) - 1;

			if (local->use_flag[eid][t - 1] == 1) {  // 白嫖时刻的可用量
				r_constrait = min((int)max_bd, r);
			}
			else r_constrait = min((int)v95, r);

			// 因缓存原因可能会算到负, 为了避免负数继续 除以 cache_p 导致溢出到正, 这里限制一下
			// 官方数据的影响实际上并不大..
			r_constrait = max(r_constrait, 0);
		}
		x1 = r_constrait - cache[eid][time] - usage[eid][time].bandwidth;
		x0 = max_bd - cache[eid][time] - usage[eid][time].bandwidth;
		ret2 = min(x0, x1);
	}
	return max(ret2, 0);
}

/*
* 计算最大可用带宽
* 仅考虑未来4个时刻的约束以提高效率
*/
int Solution::canuse_bandwidth(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN], int flag) {

	ThreadLocal* local = locals[get_local_id()];

	int max_bd = local->edge_nodes[eid].max_bandwidth;

	// 预分配标记
	// 如果是预分配 每个节点都留空一点.. (防止预分配块前面一两个时刻完全无法打入)
	if (flag) {
		max_bd = local->edge_nodes[eid].max_bandwidth * (1 - CACHE_P2) ;
	}

	double cache_p = CACHE_P1;

	int x0 = 0, x1 = 0;
	int ret1 = max_bd - usage[eid][time].bandwidth;
	int ret2 = 0;

	// 对于未来. 要多预测几天 100w最多影响 (~4.6 = 4-5天) 从第四天开始回溯条件
	// A1 = (Ct + Ut + Xt) * 0.05 + Ut+1 <= max
	// A2 = A1 * 0.05 + Ut+2 <= max
	// A3 = A2 * 0.05 + Ut+3 <= max
	// A4 = A3 * 0.05 + Ut+4 <= max
	// 从A4开始回溯

	 if (time == T - 1 ) {
		x0 = max_bd - cache[eid][time] - usage[eid][time].bandwidth;
		ret2 = x0;
	}
	else {

		int RT = min(time + 4, T - 1);

		int r_constrait = max_bd;
		for (int t = RT; t > time; --t) {

			if (local->select_20_edge[t][eid] == 1)
				cache_p = CACHE_P2;
			else
				cache_p = CACHE_P1;

			int r = (r_constrait + 1 - usage[eid][t].bandwidth) * 1.0 / max(cache_p, 1e-3) - 1;

			r_constrait = min(max_bd, r);
		}
		x1 = r_constrait - cache[eid][time] - usage[eid][time].bandwidth;
		x0 = max_bd - cache[eid][time] - usage[eid][time].bandwidth;
		ret2 = min(x0, x1 - 0);
	}

	return max(ret2, 0);
}

/*
* 更新缓存
* 仅更新未来4个时刻以提高效率
*/
void Solution::update_future_cache(int eid, int time, Timeslice usage[][MAX_T_LEN], int32 cache[][MAX_T_LEN], bool update_heap, bool update_all) {
	
	ThreadLocal* local = locals[get_local_id()];

	double cache_p = CACHE_P1;

	if (time < T && local->select_20_edge[time + 1][eid] == 1) 
		cache_p = CACHE_P2;
	else 
		cache_p = CACHE_P1;

	int cur = usage[eid][time].bandwidth + cache[eid][time];
	int next_cache = cur * cache_p;

	// 为了保证精度问题, 迭代4次 100 000 * 0.05^ x < 1 => x < 4.6, 第五次迭代的作用一定是0 但因为取整的关系. 第五次以后也可能+1
	int RT = min(time + 5, T);
	if (update_all) RT = T;

	for (int t = time + 1; t < RT; ++t) {

		if (cache[eid][t] == next_cache) break;
		else {
			cache[eid][t] = next_cache;
			cur = usage[eid][t].bandwidth + cache[eid][t];

			if (time < T && local->select_20_edge[time + 1][eid] == 1)
				cache_p = CACHE_P2;
			else
				cache_p = CACHE_P1;

			next_cache = cur * cache_p;

			if (update_heap) {
				Timeslice temp;
				temp.tid = t;
				temp.bandwidth = cur;
				local->min_heaps[eid].push_topK(temp);
				local->edge_max_used_bandwidth[eid] = max(local->edge_max_used_bandwidth[eid], cur);
			}
		}
	}
}

Timeslice Solution::get_v95(int eid) {

	ThreadLocal* local = locals[get_local_id()];
	return local->min_heaps[eid].get_K();
}

/*
* 结果合法性检查
*/
void Solution::verify() {

	ThreadLocal* local = locals[best_local_id];

	double cache_p = CACHE_P1;
	

	memset(local->edge_usage_bak, 0, sizeof(local->edge_usage_bak));
	// 检查一下 95th是否正确
	for (int i = 0; i < edge_node_count; ++i) {
		for (int t = 0; t < T; ++t) {
			local->edge_usage_bak[i][t].tid = t;
			local->edge_usage_bak[i][t].bandwidth = local->edge_usage[i][t].bandwidth + local->edge_cache[i][t];
		}
	}

	for (int i = 0; i < edge_node_count; ++i) {
		sort(local->edge_usage_bak[i], local->edge_usage_bak[i] + T, Timeslice::sort_by_bandwidth_ascend);
		int index = g95_index;
		if (local->edge_nodes[i].v95 < local->edge_usage_bak[i][index].bandwidth && local->edge_nodes[i].v95 >= base_cost) {
			LOG("v95 not match: " << i << " " << local->edge_nodes[i].v95 << " " << local->edge_usage_bak[i][index].bandwidth);
		}
	}

	// 验证一下是不是全部分配了
	for (int t = 0; t < T; ++t) {
		for (int j = 0; j < client_node_count; ++j) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				if (local->demands_alloc_flag[t][j][s] == 0) {
					ERROR("solution not found");
				}
			}
		}
	}

	// 验证缓存是否正确
	bool has_out = false;
	for (int i = 0; i < edge_node_count; ++i) {
		int max_bd = local->edge_nodes[i].max_bandwidth;

		if (local->select_20_edge[1][i] == 1)
			cache_p = CACHE_P2;
		else
			cache_p = CACHE_P1;

		int cache = (local->edge_usage[i][0].bandwidth + local->edge_cache[i][0]) * cache_p;
		for (int t = 1; t < T; ++t) {
			if (t < T && local->select_20_edge[t + 1][i] == 1)
				cache_p = CACHE_P2;
			else
				cache_p = CACHE_P1;

			// 允许有微小偏差 因为cache更新可能会引起连锁反应 如果更新的非常远, 代价会很大
			// 相应的, 可用带宽 canuse_bandwidth 那边需要做一些限制. 否则可能会引起分配不合法.
			if (abs(cache - local->edge_cache[i][t]) > 1) {
				if (!has_out) {
					LOG("cache is not correct");
					has_out = true;
				}
			}
			local->edge_cache[i][t] = cache;
			cache = (local->edge_usage[i][t].bandwidth + local->edge_cache[i][t]) * cache_p;
		}
	}

	// 验证一下有没有超分配(边缘节点是否过载):
	// 验证一下所有节点是否因缓存机制导致分配不合法.
	int each_sum_demand[135];
	for (int t = 0; t < T; ++t) {
		memset(each_sum_demand, 0, sizeof(each_sum_demand));
		for (int d = 0; d < demands_num[t]; ++d) {
			int s = sorted_demands_desc[t][d].stream_id;
			int j = sorted_demands_desc[t][d].client_id;
			if (local->alloc_mat[t][j][s] != -1) {
				int e = local->alloc_mat[t][j][s];
				if (qos_mat[j][e] >= qos_limit) {
					ERROR("qos overflow");
				}
				each_sum_demand[e] += sorted_demands_desc[t][d].demand;
			}
		}


		for (int i = 0; i < edge_node_count; ++i) {
			if (each_sum_demand[i] + local->edge_cache[i][t] > local->edge_nodes[i].max_bandwidth) {
				//
				LOG(t << " " << i << " overload: " << each_sum_demand[i] + local->edge_cache[i][t] << " " << local->edge_nodes[i].max_bandwidth);
				ERROR("edge node overload");
			}
		}
	}

	// 验证一下所有节点是否因缓存机制导致分配不合法
	//for (int i = 0; i < edge_node_count; ++i) {
	//	int max_bd = local->edge_nodes[i].max_bandwidth;
	//	for (int t = 0; t < T; ++t) {
	//		if (local->edge_usage[i][t].bandwidth + local->edge_cache[i][t] > max_bd) {
	//			ERROR("edge bd overflow!");
	//		}
	//	}
	//}

}

/*
* 为了提高效率 计算可用量以及更新缓存部分只考虑了影响最大的后几个时刻
* 缓存可能会有1的误差
* 为此, 最终需要对缓存进行矫正, 检查是否存在不合法的解, 并尝试修复不合法的解
* 主要还是通过迁移修复, 此时计算可用量从T开始逆推, 更新缓存也更新到T为止, 这样可以保证解的合法性
*/
void Solution::rescue() {

	if (best_local_id == -1) {
		ERROR("best_local_id == -1");
	}

	ThreadLocal* local = locals[best_local_id];

	vector<pair<int, int>> illegal;

	double cache_p = CACHE_P1;

	while (1) {
		bool legal = true;
		illegal.clear();
		for (int i = 0; i < edge_node_count; ++i) {
			int max_bd = local->edge_nodes[i].max_bandwidth;

			if ( local->select_20_edge[1][i] == 1)
				cache_p = CACHE_P2;
			else
				cache_p = CACHE_P1;

			int cache = (local->edge_usage[i][0].bandwidth + local->edge_cache[i][0]) * cache_p;
			for (int t = 1; t < T; ++t) {

				local->edge_cache[i][t] = cache;

				if (t < T && local->select_20_edge[t + 1][i] == 1)
					cache_p = CACHE_P2;
				else
					cache_p = CACHE_P1;


				cache = (local->edge_usage[i][t].bandwidth + local->edge_cache[i][t]) * cache_p;
				if (local->edge_usage[i][t].bandwidth + local->edge_cache[i][t] > max_bd) {
					illegal.push_back(make_pair(i, t));
					legal = false;
				}
			}
		}

		if (legal) break;

		for (auto& p : illegal) {
			int from = p.first;
			int t = p.second;
			int overflow = local->edge_nodes[from].max_bandwidth - (local->edge_usage[from][t].bandwidth + local->edge_cache[from][t]);
			//for (int d = 0; d < demands_num[t]; ++d) {
			// 小的需求往外搬
			bool is_legal = false;
			for (int d = demands_num[t] - 1; d >= 0; --d) {
				int j = sorted_demands_desc[t][d].client_id;
				int s = sorted_demands_desc[t][d].stream_id;
				int batch_demand = sorted_demands_desc[t][d].demand;
				if (batch_demand < overflow) {
					continue;
				}

				if (local->alloc_mat[t][j][s] != from) continue;

				for (int ii = 0; ii < client_nodes[j].degree; ++ii) {
					int i = client_adj[j][ii].nbr_id; // 客户节点连接的边缘节点
					int degree = client_adj[j][ii].nbr_degree;
					if (i == from) continue;
					if (local->edge_use_flag[i] == 0) continue; // 节点完全未启用 (一旦启用 代价会+V)


					int free_use = 0;

					free_use = freeuse_bandwidth(i, t, local->edge_usage, local->edge_cache, true);
					if (free_use >= batch_demand) {

						local->alloc_mat[t][j][s] = i;

						local->edge_usage[from][t].bandwidth -= batch_demand;
						local->edge_usage[i][t].bandwidth += batch_demand;

						update_future_cache(i, t, local->edge_usage, local->edge_cache, false, true);
						update_future_cache(from, t, local->edge_usage, local->edge_cache, false, true);
						is_legal = true;
						break;

					}
				}

				if (is_legal) {
					break;
				}
			}

			if (!is_legal) {
				ERROR("rescue fail!");
			}
		}
	}
}

/*
* 按该重要性顺序进行关机尝试
*/
void Solution::cal_edge_importance(int best_strategy) {

	ThreadLocal* local = locals[get_local_id()];

	int x = 0;
	for (int i = 0; i < edge_node_count; ++i) {
		if (local->edge_nodes[i].degree > 0) {
			local->sorted_edges_importance[x++] = local->edge_nodes[i];
		}
	}
	if (x != conn_edge_count) {
		ERROR("x != conn_edge_count");
	}

	sort(local->sorted_edges_importance, local->sorted_edges_importance + conn_edge_count, EdgeNode::sort_by_degree_ascend);

	//for (int i = 0; i < conn_edge_count; ++i) {
	//	local->sorted_edges_importance[i].score = local->sorted_edges_importance[i].max_bandwidth * 1.0 * (local->sorted_edges_importance[i].degree + 1e-2);
	//}
	//sort(local->sorted_edges_importance, local->sorted_edges_importance + conn_edge_count, EdgeNode::sort_by_score_ascend);

}


/*
* 多线程搜索
* 首先搜索最优的预分配或分配策略
* 在此策略的基础上尝试进行关机
*/
void Solution::gen_solution_mthreads() {

	int local_id = get_local_id();
	ThreadLocal* local = locals[local_id];
	int64 pass_t = 0;
	int count = 0;
	int best_cost = INT_MAX;
	int64 max_use_time = 0;
	int cur_ban_id = -1;
	int ban_cnt = 0;
	int best_ban_cnt = 0;


	int best_strategy = -1;
	for (int s = 0; s < strategy_num[local_id] ; ++s) {
		int64 remain = GLOBAL_MAX_TIME_MS - pass_t;
		if (remain < int64(max_use_time * 1.5)) break; // 时间不足了.. 
		auto t1 = std::chrono::steady_clock::now();

		bool res = std::invoke(ptr_alloc_func[local_id], this, s);

		if (res) {
			int cost = std::invoke(ptr_migra_func[local_id], this);

			if (cost < best_cost) {
				best_cost = cost;
				backup(local->state);
				best_strategy = s;
			}
			LOG("cost: " << cost);
		}
		reset();
		auto t2 = std::chrono::steady_clock::now();
		int64 use_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		pass_t += use_time;
		max_use_time = max(use_time, max_use_time);
		LOG ( "use_time : " << use_time );
	}

	if (best_strategy == -1) {
		//ERROR("NO SOLUTION");
		return;
	}
	LOG("best strategy: " << best_strategy);


	// 不想搜关机 就把这里关掉.
#if TRY_SHUT_DOWN != 1
	recovery(local->state);
	local->cost = best_cost;
	local->ban_cnt = 0;
	return;
#endif 

	// 计算中心节点代价
	// 首先从分配矩阵中统计每天中心节点使用量

	cal_edge_importance(best_strategy);
	memset(local->ban_flag_temp, 0, sizeof(local->ban_flag_temp));

	while (1) {
		int64 remain = GLOBAL_MAX_TIME_MS - pass_t;
		if (remain < int64(max_use_time * 1.5)) break; // 时间不足了.. 

		LOG("============");
		LOG("used time: " << pass_t << " ms");
		LOG("iter num: " << count);

		auto t1 = std::chrono::steady_clock::now();

		while (1) {
			// 没有可以ban的节点了. 跳出
			if (count >= conn_edge_count) break;
			cur_ban_id = local->sorted_edges_importance[count].id;
			local->ban_flag[cur_ban_id] = 1;
			ban_cnt++;
			break;
		}
		if (count >= conn_edge_count) break;

		bool res = std::invoke(ptr_alloc_func[local_id], this, best_strategy);

		if (!res) {
			local->ban_flag[cur_ban_id] = 0; // 注销掉产生无解的ban
			ban_cnt--;
		}
		else { // 有解才进行迁移
			int cost = std::invoke(ptr_migra_func[local_id], this);

			if (cost < best_cost) {
				best_cost = cost;
				backup(local->state);
				best_ban_cnt = ban_cnt;
				memcpy(local->ban_flag_temp, local->ban_flag, sizeof(local->ban_flag));
			}
			else {
				local->ban_flag[cur_ban_id] = 0; // 分数反向. 取消对该节点的ban. 
				ban_cnt--;
			}
			LOG("cost: " << cost);
		}

		reset();

		auto t2 = std::chrono::steady_clock::now();
		int64 use_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

		pass_t += use_time;
		max_use_time = max(use_time, max_use_time);
		LOG("cur used time: " << use_time << " ms");
		LOG("max used time: " << max_use_time << "ms");

		count++;
	}
	LOG("recovery: ");
	recovery(local->state);
	local->cost = best_cost;
	local->ban_cnt = ban_cnt;


	return;
}




// 聚合多个线程的结果
void Solution::reduce() {

	// 确定最优线程
	int32 best_cost = INT_MAX;
	for (int c = 0; c < CPU_CORES; ++c) {
		if (locals[c]->cost == -1) continue; // 无解
		if (locals[c]->cost < best_cost) {
			best_cost = locals[c]->cost;
			best_local_id = c;
		}
	}
	if (best_local_id == -1) {
		ERROR("No Solution Found");
	}
	cout << "ban cnt: " << locals[best_local_id]->ban_cnt << " best cost: " << best_cost << " from thead: " << best_local_id << endl;
}

void Solution::solve() {

	auto t1 = std::chrono::steady_clock::now();
	// 用这句 调整cpu_cores 编排好任务就是多线程
	thread threads[CPU_CORES];
	for (int i = 0; i < CPU_CORES; i++) {
		threads[i] = thread(&Solution::gen_solution_mthreads, this);
	}
	for (auto& t : threads) {
		t.join();
	}

	//gen_solution_mthreads();

	reduce(); // 聚合结果确定最优线程. 随后 get_local_id() 会默认返回最优线程的id
	rescue();
	verify();

	auto t2 = std::chrono::steady_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	cout << "total use time: " << duration.count() << endl;
}


/*
  本文件放置预处理(免费的5%部分)的策略
  =========== pre_allocation ============ 
*/

/*
* 核心思想是找到最佳的连续分块时间长度
* 使得吸收的带宽最多
*/
void Solution::pre_alloc_strategy_002() {

	ThreadLocal* local = locals[get_local_id()];

	int alloc_cnt = 0;
	int64 alloc_sum_demand = 0;

	// 0.预分配5% 这些节点取能拥有最大负载的天数打满. 所以一定是最大的5%
	for (int ii = 0; ii < edge_node_count; ++ii) {
		int i = local->sorted_edge_nodes[ii].id;                  // 按优先级序
		if (local->edge_nodes[i].degree == 0) continue;

		if (local->ban_flag[i] == 1) continue; // 如果被禁用就跳过

		// 计算该节点在所有时间片能够吸收的最大流量
		for (int t = 0; t < T; t++) {
			int32 sum = 0;
			double avg_sum = 0;
			for (int jj = 0; jj < local->edge_nodes[i].degree; ++jj) {
				int j = edge_adj[i][jj].nbr_id;
				for (int s = 0; s < stream_types_num[t]; ++s) {
					if (local->demands_alloc_flag[t][j][s] == 0) {
						int batch_demand = demands_all_types[t][j][s].demand;
						avg_sum += batch_demand * 1.0 / client_nodes[j].degree; // 除以度可以优先考虑度小但需求大的节点
						sum += batch_demand;
					}
				}
			}

			local->edge_max_load[i][t].tid = t;
			local->edge_max_load[i][t].max_load = avg_sum;
		}

		// maxbd * 0.05^x < base_cost
		double pof = log(base_cost * 1.0 / local->edge_nodes[i].max_bandwidth) / log(CACHE_P2);

		int cache_num = int(ceil(pof)); // 白嫖块末尾留给缓存的天数

		// 搜索一下最佳分块数
		vector<pair<int, int64>> slide_vec;
		vector<pair<int, int>> result;
		vector<pair<int, int>> result_temp;
		int best_block_num = 0;
		int64 max_gain = 0;
		for (int b = 2; b <= min(512, max_free_usage_count - cache_num); b = b * 3 / 2) {
			slide_vec.clear();
			result_temp.clear();
			memset(local->tag, 0, sizeof(local->tag));

			int64 slide_sum = 0;
			int slide_size = b;
			for (int t = 0; t < slide_size; ++t) {
				slide_sum += local->edge_max_load[i][t].max_load;
			}

			slide_vec.emplace_back(make_pair(slide_size, slide_sum));

			for (int t = slide_size; t < T; ++t) {
				slide_sum += local->edge_max_load[i][t].max_load;
				slide_sum -= local->edge_max_load[i][t - slide_size].max_load;
				slide_vec.emplace_back(make_pair(t + 1, slide_sum));
			}

			sort(slide_vec.begin(), slide_vec.end(), [](pair<int, int64>& p1, pair<int, int64>& p2) -> bool {
				if (p1.second != p2.second) return p1.second > p2.second;
				else return p1.first < p2.first;
				});

			int canuse = max_free_usage_count;
			int real_use_count = 0;
			int cache_use_count = 0;

			int64 gain = 0;
			for (int it = 0; it < slide_vec.size(); ++it) {
				int R = slide_vec[it].first;
				int L = R - slide_size;
				if (local->tag[L] != 0 || local->tag[R - 1] != 0) continue;
				if (canuse >= (slide_size + cache_num)) {
					gain += slide_vec[it].second;
					memset(local->tag + L, 0xFF, slide_size + cache_num);

					result_temp.push_back(make_pair(L, R));
					canuse -= slide_size + cache_num;

					real_use_count += slide_size;
					cache_use_count += cache_num;
				}
				else {
					R = L + canuse - cache_num;
					int64 temp_gain = 0;
					for (int inner = L; inner < R; ++inner) {
						temp_gain += local->edge_max_load[i][inner].max_load;
					}
					if (temp_gain > 0) {
						gain += temp_gain;
						result_temp.push_back(make_pair(L, R));
						canuse = 0;

						real_use_count += R - L;
						cache_use_count += cache_num;
					}
				}
				if (canuse <= cache_num) break;
			}
			if (gain > max_gain) {
				best_block_num = slide_size;
				result.assign(result_temp.begin(), result_temp.end());
				max_gain = gain;
			}
		}
		// 根据分块结果 进行实际分配
		for (auto& p : result) {
			int L = p.first, R = p.second;
			for (int t = L; t < R; ++t) {
				int free;
				int canuse_flag = (t == L ? 1 : 2);
				free = canuse_bandwidth(i, t, local->edge_usage, local->edge_cache, canuse_flag);
				for (int d = 0; d < demands_num[t]; ++d) {
					int j = sorted_demands_type_desc[t][d].client_id;
					int s = sorted_demands_type_desc[t][d].stream_id;
					int batch_demand = sorted_demands_type_desc[t][d].demand;

					if (qos_mat[j][i] >= qos_limit) continue;
					if (free == 0) break;
					if (local->demands_alloc_flag[t][j][s] != 0) continue; // 已分配过
					if (batch_demand == 0) {
						local->demands_alloc_flag[t][j][s] = 1;
						continue;
					}
					if (free >= batch_demand) {
						alloc(i, j, s, t, false, free);
						local->each_stream_usage[t][i][s] += batch_demand / local->edge_nodes[i].degree;
						local->each_stream_count[t][i][s] ++;
						alloc_cnt++;
						alloc_sum_demand += batch_demand;
						free -= batch_demand;
						local->use_flag[i][t] = 1;
						local->edge_use_flag[i] = 1;
						local->demands_alloc_flag[t][j][s] = 1;
					}
				}
			}
		}
	}
	LOG("pre alloc cnt: " << alloc_cnt);
	LOG("pre alloc sum demand: " << alloc_sum_demand);
}

/*
* 1.对不同时刻 按吸收流量大小降序排序, 记为ranks
* 2.对于某个时刻t, 判断t-1 t+1是否吸收流量也足够大(排名在ranks足够靠前) 如果是, 则加长该块长度
* 3.处理一些边界情况
* 4.白嫖部分可能没有打满, 因此cache_num的每个时刻可能还有一些免费余量(这个免费余量总量是可以推测的) 将这些余量也利用起来
*/
void Solution::pre_alloc_strategy_001() {
	int alloc_cnt = 0;
	int64 alloc_sum_demand = 0;

	ThreadLocal* local = locals[get_local_id()];

	// 0.预分配5% 这些节点取能拥有最大负载的天数打满. 所以一定是最大的5%
	for (int ii = 0; ii < edge_node_count; ++ii) {
		int i = local->sorted_edge_nodes[ii].id;                  // 按优先级序
		int max_bd = local->sorted_edge_nodes[ii].max_bandwidth;
		if (local->edge_nodes[i].degree == 0) continue;

		if (local->ban_flag[i] == 1) continue; // 如果被禁用就跳过

		// 计算该节点在所有时间片能够吸收的最大流量
		for (int t = 0; t < T; t++) {
			int32 sum = 0;
			double avg_sum = 0;
			for (int jj = 0; jj < local->edge_nodes[i].degree; ++jj) {
				int j = edge_adj[i][jj].nbr_id;
				for (int s = 0; s < stream_types_num[t]; ++s) {
					if (local->demands_alloc_flag[t][j][s] == 0) {
						int batch_demand = demands_all_types[t][j][s].demand;
						avg_sum += batch_demand * 1.0 / client_nodes[j].degree; // 算这个可以优先考虑度小但需求大的节点
						sum += batch_demand;
					}
				}
			}

			local->edge_max_load[i][t].tid = t;
			local->edge_max_load[i][t].max_load = avg_sum;
			local->edge_max_load_2[i][t].tid = t;
			local->edge_max_load_2[i][t].sum_demand = sum;
		}
		sort(local->edge_max_load[i], local->edge_max_load[i] + T, Timeslice::sort_by_max_load_descend);

		vector<pair<int, int>> rank_index;
		rank_index.resize(T);
		// 排序索引取 edge_max_load (优先考虑度小但需求大的节点)
		// 实际值取吸收流量的大小 edge_max_load_2 
		// 线上玄学上了一点分
		for (int t = 0; t < T; ++t) {
			int tid = local->edge_max_load[i][t].tid;
			rank_index[tid] = { t , local->edge_max_load_2[i][tid].bandwidth };
		}

		// maxbd * 0.05^x < base_cost
		double pof = log(base_cost * 1.0 / local->edge_nodes[i].max_bandwidth) / log(CACHE_P2);
		int cache_num = int(ceil(pof)); // 白嫖块末尾留给缓存的天数

		vector<pair<int, int>> result;
		vector<pair<int, int>> result2; //缓存

		int64 gain_max = 0;
		double coef = 0.1;
		uint8 tag2[MAX_T_LEN] = { 0 };
		while (coef < 1.01) {
			int64 gain = 0;
			int time_thre = T * coef;
			int canuse_time = max_free_usage_count;
			int all_use_day = 0;
			memset(local->tag, 0, sizeof(local->tag));

			for (int t = 0; t < T; ) {

				int tid = local->edge_max_load[i][t].tid;
				if (local->tag[tid] == 0xFF) {
					++t;
					continue;
				}
				if (t > time_thre) break;
				if (canuse_time <= cache_num) break;
				int t_start = tid;
				int t_end = tid;

				while (true) {
					//如果可用空间不足退出
					if (canuse_time <= cache_num) {
						break;
					}
					bool move_right = ! (t_end >= T || local->tag[t_end] == 0xFF || 
						(rank_index[t_end].first > time_thre && rank_index[t_end].second < max_bd));
					bool move_left = ! (t_start < 0 || local->tag[t_start] == 0xFF || 
						(rank_index[t_start].first > time_thre && rank_index[t_start].second < max_bd));
					
					//如果不能往前走也不能往后走
					if ( !move_right && !move_left) {
						t_end--;
						t_start++;
						break;
					}
					//如果可以往前也可以往后
					else if (move_right && move_left) {
						if (rank_index[t_end].second >= rank_index[t_start].second) {
							all_use_day++;
							gain += rank_index[t_end].second;
							t_end++;
							canuse_time--;
						}
						else {
							all_use_day++;
							gain += rank_index[t_start].second;
							t_start--;
							canuse_time--;
						}
					}
					//如果可以往后但是不能往前
					else if (move_right && !move_left) {
						all_use_day++;
						gain += rank_index[t_end].second;
						t_end++;
						canuse_time--;
					}
					//可以往前但是不能往后
					else if (move_left && !move_right) {
						all_use_day++;
						gain += rank_index[t_start].second;
						t_start--;
						canuse_time--;
					}
				}
				// 边界约束与标记更新
				if (t_end >= T) {
					t_end = T - 1;
				}
				if (t_start < 0) {
					t_start = 0;
				}
				if ((t_end + cache_num) >= T) {
					t_end = T;
					memset(local->tag + t_start, 0xFF, (T - t_start) * sizeof(uint8));
					++t;
				}
				else if (local->tag[t_end + cache_num + 1] == 0xFF) { // 尝试和之前的结果连成片
					t_end = t_end + cache_num + 1;
					memset(local->tag + t_start, 0xFF, (t_end - t_start + cache_num) * sizeof(uint8));
					++t;
				}
				else {
					t_end += 1;
					memset(local->tag + t_start, 0xFF, (t_end - t_start + cache_num) * sizeof(uint8));
					++t;
				}

				// 计算剩余可用, 如果超过canuse的最大值也保护一下
				canuse_time = max_free_usage_count;
				for (int tt = 0; tt < T; tt++) {
					if (local->tag[tt] == 0xFF) {
						if (canuse_time > 0) canuse_time--;
						else local->tag[tt] = 0;
					}
				}
			}

			int all_cache_day = 0;
			for (int tt = 0; tt < T; tt++) {
				if (local->tag[tt] == 0xFF) {
					all_cache_day++;
				}
			}
			all_cache_day = all_cache_day - all_use_day;

			// 预估可吸收的流量
			if ((gain + all_cache_day * base_cost) > gain_max) {
				gain_max = (gain + all_cache_day * base_cost);
				memcpy(tag2, local->tag, sizeof(local->tag));
			}
			coef += 0.025;
		}

		// 统计白嫖块的区间 rbegin 和 rend一一对应
		vector<int> rbegin;
		vector<int> rend;
		if (tag2[0] == 0xFF) {
			rbegin.push_back(0);
		}
		for (int tt = 0; tt < T - 1; tt++) {
			if (tag2[tt] == 0x00 && tag2[tt + 1] == 0xFF) {
				rbegin.push_back(tt + 1);
			}
		}
		for (int tt = 0; tt < T - 1; tt++) {
			if (tag2[tt] == 0xFF && tag2[tt + 1] == 0x00) {
				rend.push_back(tt + 1);
			}
		}
		if (tag2[T - 1] == 0xFF) {
			rend.push_back(T);
		}
		if (rbegin.size() != rend.size()) {
			ERROR("ERROR");
		}

		// 统计区分白嫖块中实际打满的时刻和作为缓存的时刻
		// 并记录到result/result2中
		for (int i = 0; i < rbegin.size(); i++) {
			pair<int, int> temp_pair;
			temp_pair.first = rbegin[i];
			if (rend[i] == T) {
				temp_pair.second = T;
			}
			else {
				temp_pair.second = rend[i] - cache_num;
				// space 是预估不会引起整体v95上涨的最大可用量, 从base_cost开始逆推
				int space = base_cost;
				for (int c = 0; c < cache_num; c++) {
					pair<int, int> temp_pair2;
					temp_pair2.first = rend[i] - (c + 1);
					temp_pair2.second = space;
					space = space / CACHE_P2;
					result2.push_back(temp_pair2);
				}
			}
			result.push_back(temp_pair);
		}

		sort(result.begin(), result.end(), [](pair<int, int >& a1, pair<int, int >& a2) { return a1.first < a2.first;});
		sort(result2.begin(), result2.end(), [](pair<int, int >& a1, pair<int, int >& a2) { return a1.first < a2.first;});

		for (auto& p : result) {
			int L = p.first, R = p.second;
			for (int t = L; t < R; ++t) {
				int free = canuse_bandwidth(i, t, local->edge_usage, local->edge_cache, 1);
				for (int d = 0; d < demands_num[t]; ++d) {
					int j = sorted_demands_type_desc[t][d].client_id;
					int s = sorted_demands_type_desc[t][d].stream_id;
					int batch_demand = sorted_demands_type_desc[t][d].demand;

					if (qos_mat[j][i] >= qos_limit) continue;
					if (free == 0) break;
					if (local->demands_alloc_flag[t][j][s] != 0) continue; // 已分配过
					if (batch_demand == 0) {
						local->demands_alloc_flag[t][j][s] = 1;
						continue;
					}
					if (free >= batch_demand) {
						alloc(i, j, s, t, false, free);
						local->each_stream_usage[t][i][s] += batch_demand / local->edge_nodes[i].degree;
						local->each_stream_count[t][i][s] ++;
						alloc_cnt++;
						alloc_sum_demand += batch_demand;
						free -= batch_demand;
						local->use_flag[i][t] = 1;
						local->edge_use_flag[i] = 1;
						local->demands_alloc_flag[t][j][s] = 1;
					}
				}
			}
		}
		// 可打满的时刻打满后, 重新考虑这些白嫖块末尾的缓存是否完全利用上了
		for (auto& p : result2) {
			int t = p.first;
			int free;
			int canuse_flag = 1;
			free = min(canuse_bandwidth(i, t, local->edge_usage, local->edge_cache, 1), p.second - local->edge_usage[i][t].bandwidth - local->edge_cache[i][t]);

			for (int d = 0; d < demands_num[t]; ++d) {
				int j = sorted_demands_type_desc[t][d].client_id;
				int s = sorted_demands_type_desc[t][d].stream_id;
				int batch_demand = sorted_demands_type_desc[t][d].demand;

				if (qos_mat[j][i] >= qos_limit) continue;
				if (free == 0) break;
				if (local->demands_alloc_flag[t][j][s] != 0) continue; // 已分配过
				if (batch_demand == 0) {
					local->demands_alloc_flag[t][j][s] = 1;
					continue;
				}

				if (free >= batch_demand) {
					alloc(i, j, s, t, false, free);
					local->each_stream_usage[t][i][s] += batch_demand / local->edge_nodes[i].degree;
					local->each_stream_count[t][i][s] ++;
					alloc_cnt++;
					alloc_sum_demand += batch_demand;
					free -= batch_demand;
					local->edge_use_flag[i] = 1;
					local->demands_alloc_flag[t][j][s] = 1;
				}
			}
		}
	}
	LOG("pre alloc cnt: " << alloc_cnt);
	LOG("pre alloc sum demand: " << alloc_sum_demand);

}


/*
  本文件放置分配策略

  ====== allocation =======
*/

bool Solution::alloc_strategy_001() {

	ThreadLocal* local = locals[get_local_id()];

	int cannot_alloc = 0;
	int alloc_cnt = 0;
	int64 alloc_sum_demand = 0;
	for (int t = 0; t < T; ++t) {
		if (t > 0) {
			select_today_20(t);
			recal_cache(t - 1);
		}

		int edge_remains[MAX_EDGENODE_NUM];
		for (int i = 0; i < edge_node_count; ++i) {
			edge_remains[i] = canuse_bandwidth(i, t, local->edge_usage, local->edge_cache);
		}

		for (int d = 0; d < demands_num[t]; ++d) {
			int batch_demand;
			int j;
			int s;
			j = sorted_demands_desc[t][d].client_id;
			s = sorted_demands_desc[t][d].stream_id;
			batch_demand = sorted_demands_desc[t][d].demand;
			
			if (local->demands_alloc_flag[t][j][s] == 1) continue; // 跳过已分配
			if (batch_demand == 0) {
				local->demands_alloc_flag[t][j][s] = 1;
				continue;
			}

			// 找最合适的边缘
			int sel_node = -1;
			int state = INT_MAX;
			int max_score_state[6];
			memset(max_score_state, 0, sizeof(max_score_state));
			max_score_state[5] = INT_MIN;

			for (int ii = 0; ii < client_nodes[j].degree; ++ii) {
				int i = client_adj[j][ii].nbr_id;				//边缘的结点
				int degree = client_adj[j][ii].nbr_degree;		//边缘的度
				if (degree == 0) continue;
				if (local->edge_use_flag[i] == 0) continue; // 未启用就再也不启用
				int remain = edge_remains[i];
				if (remain < batch_demand) continue;
				int free_use = 0;

				if (local->use_flag[i][t] == 1) {  // 属于白嫖且可以打满的时刻
					free_use = canuse_bandwidth(i, t, local->edge_usage, local->edge_cache);
					if (free_use >= batch_demand) {
						if (local->each_stream_usage[t][i][s] > 0 && local->each_stream_usage[t][i][s] > max_score_state[0]) {
							state = 1;	//状态1为可白嫖 且有该类型流的
							sel_node = i;
							max_score_state[0] = local->each_stream_usage[t][i][s];
						}
						if (state == 1) continue;
						if (free_use / degree > max_score_state[1]) { //状态2为可白嫖 且无该类型流的,找剩余空间最多的,且度最小的
							state = 2;
							sel_node = i;
							max_score_state[1] = free_use / degree;
						}
					}
				}
				if (state <= 2) continue;
			
				free_use = min(max((get_v95(i).bandwidth), base_cost) 
					- local->edge_usage[i][t].bandwidth - local->edge_cache[i][t], remain);

				if (free_use >= batch_demand) {
					if (local->each_stream_usage[t][i][s] > 0 && local->each_stream_usage[t][i][s] > max_score_state[2]) {
						state = 3;	//状态3为存在免费额度 且有该类型流的
						sel_node = i;
						max_score_state[2] = local->each_stream_usage[t][i][s];
					}
					if (state <= 3) continue;

					if (free_use / degree > max_score_state[3]) { //状态4为存在免费额度, 找剩余空间最多的, 度最小的
						state = 4;
						sel_node = i;
						max_score_state[3] = free_use / degree;
					}
				}
				if (state <= 4) continue;

				if (free_use > 0) { 
					int score = sqrt(degree) * free_use;
					if (score > max_score_state[4]) { // 状态5为普通扩容, 免费额度不足以容纳当前需求, 找度大剩余空间大的(玄)
						max_score_state[4] = score;
						sel_node = i;
						state = 5;
					}
				}
				if (state <= 5) continue;
				
				// 状态6为特殊扩容, 在处理前几个时刻(t=0附近) v95值还很小, 但存在一些白嫖块中的缓存(这些时刻并不能打满, 但由于v95还很小 故当前又属于后5%)
				// 不能让这些时刻打得过多, 否则它们对下一天的缓存影响会增大
				// 如果非要在这些节点选, 我们尽量选使用量小, 度大的
				// 这里的扩容选择影响还是挺大的
				state = 6;
				int score;
				score = -(local->edge_usage[i][t].bandwidth + local->edge_cache[i][t]) / (degree + 1e-3);
				if (score > max_score_state[5]) {
					max_score_state[5] = score;
					sel_node = i;
				}
			}

			if (sel_node == -1) {
				LOG("not solution found");
				cannot_alloc++;
				return false;
			}
			else {
				alloc(sel_node, j, s, t, true, edge_remains[sel_node]);
				local->edge_use_flag[sel_node] = 1;
				local->demands_alloc_flag[t][j][s] = 1;
				edge_remains[sel_node] -= batch_demand;
				local->each_stream_usage[t][sel_node][s] += batch_demand / local->edge_nodes[sel_node].degree;
				local->each_stream_count[t][sel_node][s] ++;
				alloc_cnt++;
				alloc_sum_demand += batch_demand;
			}
		}
	}
	LOG("alloc: " << alloc_cnt);
	LOG("alloc sum demand: " << alloc_sum_demand);
	return true;
}


bool Solution::alloc_strategy_002() {

	ThreadLocal* local = locals[get_local_id()];

	int cannot_alloc = 0;
	int alloc_cnt = 0;
	int64 alloc_sum_demand = 0;
	int edge_remains[MAX_EDGENODE_NUM];

	for (int t = 0; t < T; ++t) {
		if (t > 0) {
			select_today_20(t);
			recal_cache(t - 1);
		}
		for (int i = 0; i < edge_node_count; ++i) {
			edge_remains[i] = canuse_bandwidth(i, t, local->edge_usage, local->edge_cache);
		}
		for (int d = 0; d < demands_num[t]; ++d) {
			int batch_demand;
			int j;
			int s;

			j = sorted_demands_desc[t][d].client_id;
			s = sorted_demands_desc[t][d].stream_id;
			batch_demand = sorted_demands_desc[t][d].demand;
			if (local->demands_alloc_flag[t][j][s] == 1) continue; // 跳过已分配
			if (batch_demand == 0) {
				local->demands_alloc_flag[t][j][s] = 1;
				continue;
			}

			Package state;
			state.level = INT_MAX;
			state.score = -1;
			state.score2 = -1;
			state.select_node = -1;

			for (int ii = 0; ii < client_nodes[j].degree; ++ii) {
				int i = client_adj[j][ii].nbr_id;				//边缘的结点
				int degree = client_adj[j][ii].nbr_degree;		//边缘的度
				if (degree == 0) continue;
				if (local->edge_use_flag[i] == 0) continue; // 未启用就再也不启用
				int remain = edge_remains[i];
				if (remain < batch_demand) continue;

				int free_use = 0;
				if (local->use_flag[i][t] == 1) {  // 属于白嫖且能够打满的时刻
					free_use = remain;
					double score = free_use * 1.0 / degree;
					if (free_use >= batch_demand) {
						if ((state.level > 1 && local->each_stream_usage[t][i][s] > 0)
							|| (state.level == 1 && local->each_stream_usage[t][i][s] > state.score)
							|| (state.level == 1 && local->each_stream_usage[t][i][s] == state.score &&
								local->edge_nodes[i].max_bandwidth > state.score2)) {
							// 状态1: 可白嫖且有同流 找最大同流
							state.level = 1;
							state.select_node = i;
							state.score = local->each_stream_usage[t][i][s];
							state.score2 = local->edge_nodes[i].max_bandwidth;
						}
						else if (state.level > 2 || (state.level == 2 && score > state.score)
							|| (state.level == 2 && score == state.score && local->edge_nodes[i].max_bandwidth > state.score2)) {
							//状态2: 可白嫖但无同流 找最大可用(优先度小)
							state.level = 2;
							state.select_node = i;
							state.score = score;
							state.score2 = local->edge_nodes[i].max_bandwidth;
						}
					}
				}
				else {
					int usage = local->edge_usage[i][t].bandwidth + local->edge_cache[i][t];
					free_use = max((get_v95(i).bandwidth), base_cost) - usage;
					free_use = min(free_use, remain);
					double score = free_use * 1.0 / degree;
					if (free_use >= batch_demand) {
						if ((state.level > 3 && local->each_stream_usage[t][i][s] > 0)
							|| (state.level == 3 && local->each_stream_usage[t][i][s] > state.score)
							|| (state.level == 3 && local->each_stream_usage[t][i][s] == state.score &&
								local->edge_nodes[i].max_bandwidth > state.score2)) {
							//状态3: 不扩容且有同流 找最大同流
							state.level = 3;
							state.select_node = i;
							state.score = local->each_stream_usage[t][i][s];
							state.score2 = edge_nodes[i].max_bandwidth;
						}

						else if (state.level > 4 || (state.level == 4 && score > state.score) ||
							(state.level == 4 && score == state.score && local->edge_nodes[i].max_bandwidth > state.score2)) {
							//状态4: 不扩容且无同流 找最大可用(优先度小)
							state.level = 4;
							state.select_node = i;
							state.score = score;
							state.score2 = local->edge_nodes[i].max_bandwidth;
						}
					}

					else if (free_use > 0 ) { // 普通扩容
						double score = sqrt(degree) * free_use;
						if (state.level > 5 || (score > state.score || (score == state.score && local->edge_nodes[i].max_bandwidth > state.score2))) {
							state.level = 5;
							state.select_node = i;
							state.score = score;
							state.score2 = local->edge_nodes[i].max_bandwidth;
						}
					}
					else  { // 缓存扩容
						double score = -(local->edge_usage[i][t].bandwidth + local->edge_cache[i][t]) * 1.0 / (degree );
						if (state.level > 6 || (state.level == 6 && score > state.score)
							|| (state.level == 6 && score == state.score && local->edge_nodes[i].max_bandwidth > state.score2)) {
							state.level = 6;
							state.select_node = i;
							state.score = score;
							state.score2 = edge_nodes[i].max_bandwidth;
						}

					}
				}
			}
			int sel_node = state.select_node;
			if (sel_node == -1) {
				LOG("not solution found");
				cannot_alloc++;
				return false;
			}
			else {
				alloc(sel_node, j, s, t, true, edge_remains[sel_node]);
				local->edge_use_flag[sel_node] = 1;
				local->demands_alloc_flag[t][j][s] = 1;
				edge_remains[sel_node] -= batch_demand;
				local->each_stream_usage[t][sel_node][s] += batch_demand / local->edge_nodes[sel_node].degree;
				local->each_stream_count[t][sel_node][s] ++;
				alloc_cnt++;
				alloc_sum_demand += batch_demand;

			}
		}
	}
	LOG("alloc: " << alloc_cnt);
	LOG("alloc sum demand: " << alloc_sum_demand);
	return true;
}


/*
  本文件放置迁移策略

  ============ migration ===============
*/


/*
* migration_edge 降低边缘的成本
* migration_center 降低中心的成本
* 随后确定了中心的后5%, 进而再次尝试降低一些边缘成本
* 最后一轮迁移边缘节点, 在中心节点的5%确定后, 没有严格锁住中心成本的上升, 这样线上效果更好一些
*/
int Solution::migration_process() {

	ThreadLocal* local = locals[get_local_id()];

	// 根据堆更新v95. 此后固定该v95为基准
	for (int i = 0; i < edge_node_count; ++i) {
		auto temp = get_v95(i);
		local->edge_nodes[i].v95_tid = temp.tid;
		local->edge_nodes[i].v95 = temp.bandwidth;
	}

	int prev_cost = 0, cost = -1;


	memset(local->center_use, 0, sizeof(local->center_use));

	migration_edge();
	migration_center();
	migration_edge();

	//cost = cal_cost(false, false);
	cost = resort_cal_cost();

	cout << "center " << int(local->center_v95 * center_cost) << endl;
	cout << "after migration cost:" << cost << endl;

	return cost;
}

/*
* 计算代价更大
* 不断执行迁移1和迁移2
* 直到分数不变或更差
* 线上效果不如process_1
*/
int Solution::migration_process_2() {

	ThreadLocal* local = locals[get_local_id()];

	// 根据堆更新v95. 此后固定该v95为基准
	for (int i = 0; i < edge_node_count; ++i) {
		auto temp = get_v95(i);
		local->edge_nodes[i].v95_tid = temp.tid;
		local->edge_nodes[i].v95 = temp.bandwidth;
	}

	int prev_cost = 0, cost = -1;

	memset(local->center_use, 0, sizeof(local->center_use));
	
	backup(local->state2);
	prev_cost = resort_cal_cost();
	while (1) {
		migration_edge();
		migration_center();
		cost = resort_cal_cost();
		if (cost < prev_cost) {
			backup(local->state2);
			prev_cost = cost;
		}
		else break;
	}
	recovery(local->state2);

	return prev_cost;
}



/*
* 采用对顶堆维护迁移对象
* 核心思路: 选择一个利用率低的节点, 选择其95分位的时刻, 不断搬出该时刻的需求; 然后选择新的95分位时刻, 不断往复.
*/
void Solution::migration_edge() {

	ThreadLocal* local = locals[get_local_id()];

	unordered_set<int> record;


	local->migra_min_heap.set_K(max_free_usage_count + 1);
	local->migra_max_heap.set_K(T - (max_free_usage_count + 1));

	while (1) {

		int min_score = INT_MAX, from_id = -1;
		// 寻找利用率最低的节点
		for (int i = 0; i < edge_node_count; ++i) {
			int sum_bd = 0;
			if (local->edge_nodes[i].degree == 0) continue;
			if (record.find(i) != record.end()) continue;

			for (int t = 0; t < T; ++t) {
				sum_bd += local->edge_usage[i][t].bandwidth + local->edge_cache[i][t];
			}
			double score = sum_bd * 1.0 / local->edge_nodes[i].v95;
			if (min_score > score) {
				min_score = score;
				from_id = i;
			}
		}
		if (from_id == -1) break; // 说明所有已启用的节点都已经尝试迁移过了
		record.insert(from_id);

		int cur_v95th = local->edge_nodes[from_id].v95, cur_v95th_tid = local->edge_nodes[from_id].v95_tid;
		local->migra_min_heap.reset();
		local->migra_max_heap.reset();

		// 维护初始对顶堆
		for (int tt = 0; tt < T; tt++) {
			int cur_usage = local->edge_usage[from_id][tt].bandwidth + local->edge_cache[from_id][tt];
			if (local->migra_min_heap.get_K().bandwidth < cur_usage) {
				if (local->migra_min_heap.n == local->migra_min_heap.K) {
					local->migra_max_heap.push_topK(local->migra_min_heap.get_K());
				}
				local->migra_min_heap.push_topK(Timeslice{ tt, cur_usage });
			}
			else {
				local->migra_max_heap.push_topK(Timeslice{ tt, cur_usage });
			}
		}

		while (true)
		{
			int t = local->migra_min_heap.get_K().tid;
			int init_cost = local->migra_min_heap.get_K().bandwidth;

			for (int d = 0; d < demands_num[t]; ++d) {
				//for (int d = demands_num[t] - 1; d >= 0; --d) {
				int j = sorted_demands_desc[t][d].client_id;
				int s = sorted_demands_desc[t][d].stream_id;
				if (local->alloc_mat[t][j][s] != from_id) continue;
				int batch_demand = sorted_demands_desc[t][d].demand;
				if (local->edge_usage[from_id][t].bandwidth + local->edge_cache[from_id][t] <= base_cost) {
					break;
				}

				int score = 0;
				int best_to = -1;

				for (int ii = 0; ii < client_nodes[j].degree; ++ii) {
					int i = client_adj[j][ii].nbr_id; // 客户节点连接的边缘节点
					int degree = client_adj[j][ii].nbr_degree;
					if (i == from_id) continue;
					if (local->edge_use_flag[i] == 0) continue; // 节点完全未启用 (一旦启用 代价会+V)


					int free_use = 0;

					// 使用下面的约束条件 且 外面标记center_use[t] = 1后, 可以锁住中心上涨 (但线上效果不好..)
					//if (local->center_use[t] == 0) continue;
					//if (local->each_stream_count[t][i][s] == 0 && local->each_stream_count[t][from_id][s] > 1 ) {
					//	continue;
					//}

					if (local->each_stream_count[t][i][s] == 0 && local->center_use[t] != 1 && local->each_stream_count[t][from_id][s] > 1) {
						continue;
					}

					free_use = freeuse_bandwidth(i, t, local->edge_usage, local->edge_cache);

					if (free_use < batch_demand) {
						continue;
					}
					//优先度小的
					if (free_use / degree > score) {
						best_to = i;
						score = free_use / degree;
					}
				}
				if (best_to != -1) {
					move(from_id, best_to, j, s, t, false);
					local->each_stream_count[t][from_id][s]--;
					local->each_stream_count[t][best_to][s]++;
				}
			}

			// 更新堆. 不同于中心节点对顶堆的维护, 边缘的维护需要考虑缓存
			for (int t0 = t; t0 < min(T, t + 5); ++t0) {
				int cur_usage = local->edge_usage[from_id][t0].bandwidth + local->edge_cache[from_id][t0];

				// 原数据在小顶堆中
				if (local->migra_min_heap.mp[t0] != -1) {
					local->migra_min_heap.push_topK(Timeslice{ t0, cur_usage });
				}
				else { // 在大顶堆中
					local->migra_max_heap.push_topK(Timeslice{ t0, cur_usage });
				}

				// 判断是否要交换
				if (local->migra_max_heap.get_K().bandwidth > local->migra_min_heap.get_K().bandwidth) {
					auto temp = local->migra_min_heap.get_K();
					local->migra_min_heap.push_topK(local->migra_max_heap.get_K()); // 一定会成功 (小顶堆 顶部 小于 大顶堆顶部)
					local->migra_max_heap.push_topK(temp);
				}
			}

			int cur_usage = local->edge_usage[from_id][t].bandwidth + local->edge_cache[from_id][t];
			if (cur_usage == init_cost || cur_usage >= local->migra_max_heap.get_K().bandwidth) { // 没迁移动 或 比大顶堆首个还要大
				break;
			}
		}

		cur_v95th = local->migra_min_heap.get_K().bandwidth;
		cur_v95th_tid = local->migra_min_heap.get_K().tid;

		// 在决赛的新需求下, 如果继续更换use_flag(也就是后5%的选择) 会影响之前选择的保留率为0.01的20个节点
		// 但影响可能不是很大-w- .. 就不修了
		memset(local->use_flag[from_id], 0, sizeof(local->use_flag[0]));
		// 更新后5%的标记
		for (int tt = 0; tt < T; tt++) {
			if (local->migra_min_heap.mp[tt] != -1 && tt != cur_v95th_tid) local->use_flag[from_id][tt] = 1;
		}

		local->edge_nodes[from_id].v95 = cur_v95th;
		local->edge_nodes[from_id].v95_tid = cur_v95th_tid;
	}
	return;
}

/*
* 采用对顶堆维护迁移对象
* 核心思路: 选择其中心节点95分位的时刻, 不断搬出在边缘节点中同类型数目只有1种的客户需求
*/
void Solution::migration_center()
{
	ThreadLocal* local = locals[get_local_id()];

	int preday95 = -1;
	int postday95 = -1;

	local->migra_min_heap.reset();
	local->migra_max_heap.reset();
	local->migra_min_heap.set_K(max_free_usage_count + 1);
	local->migra_max_heap.set_K(T - (max_free_usage_count + 1));

	//计算中心节点5%是哪几天
	//首先从分配矩阵中统计每天中心节点使用量
	for (int t = 0; t < T; ++t) {
		memset(local->temp_each_stream_maxval, 0, sizeof(local->temp_each_stream_maxval));

		for (int j = 0; j < client_node_count; ++j) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				int i = local->alloc_mat[t][j][s];
				if (i != NOT_ALLOCATION) {
					int batch_demand = demands_all_types[t][j][s].demand;
					local->temp_each_stream_maxval[i][s] = max(local->temp_each_stream_maxval[i][s], batch_demand);
				}
			}
		}
		int32 sum_bd = 0;
		for (int i = 0; i < edge_node_count; ++i) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				sum_bd += local->temp_each_stream_maxval[i][s];
			}
		}

		if (local->migra_min_heap.get_K().bandwidth < sum_bd) {
			if (local->migra_min_heap.n == local->migra_min_heap.K) {
				local->migra_max_heap.push_topK(local->migra_min_heap.get_K());
			}
			local->migra_min_heap.push_topK(Timeslice{ t, sum_bd });
		}
		else {
			local->migra_max_heap.push_topK(Timeslice{ t, sum_bd });
		}
	}

	local->center_v95 = local->migra_min_heap.get_K().bandwidth;
	preday95 = local->migra_min_heap.get_K().tid;

	while (true) {
		//迁移
		int t = preday95;
		for (int d = 0; d < demands_num[t]; ++d) {
			//for (int d = demands_num[t] - 1; d >= 0; --d) {
			int j = sorted_demands_desc[t][d].client_id;
			int s = sorted_demands_desc[t][d].stream_id;
			int batch_demand = sorted_demands_desc[t][d].demand;
			int from_id = local->alloc_mat[t][j][s];
			if (from_id == NOT_ALLOCATION) continue;
			if (local->each_stream_count[t][from_id][s] == 1) {
				int score = 0;
				int best_to = -1;
				for (int ii = 0; ii < client_nodes[j].degree; ++ii) {
					int i = client_adj[j][ii].nbr_id; // 客户节点连接的边缘节点
					int degree = client_adj[j][ii].nbr_degree;
					if (i == from_id) continue;
					if (local->edge_use_flag[i] == 0) continue; // 节点完全未启用 (一旦启用 代价会+V)

					if (local->each_stream_count[t][i][s] == 0 ) {
						continue;
					}

					int free_use = freeuse_bandwidth(i, t, local->edge_usage, local->edge_cache);

					if (free_use < batch_demand) continue;
					//优先剩余空间大 & 度小的
					if (free_use / degree > score) {
						best_to = i;
						score = free_use / degree;
					}
				}
				if (best_to != -1) {
					move(from_id, best_to, j, s, t, false);
					local->each_stream_count[t][from_id][s]--;
					local->each_stream_count[t][best_to][s]++;
				}
			}

		}

		// 重新统计迁移后该时刻的中心节点代价
		memset(local->temp_each_stream_maxval, 0, sizeof(local->temp_each_stream_maxval));
		for (int j = 0; j < client_node_count; ++j) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				int i = local->alloc_mat[t][j][s];
				if (i != NOT_ALLOCATION) {
					int batch_demand = demands_all_types[t][j][s].demand;
					local->temp_each_stream_maxval[i][s] = max(local->temp_each_stream_maxval[i][s], batch_demand);
				}
			}
		}
		int32 sum_bd = 0;
		for (int i = 0; i < edge_node_count; ++i) {
			for (int s = 0; s < stream_types_num[t]; ++s) {
				sum_bd += local->temp_each_stream_maxval[i][s];
			}
		}

		// 维护对顶堆
		local->migra_min_heap.push_topK(Timeslice{ t, sum_bd });
		if (local->migra_max_heap.get_K().bandwidth > sum_bd) {
			local->migra_min_heap.push_topK(local->migra_max_heap.get_K());
			local->migra_max_heap.push_topK(Timeslice{ t, sum_bd });
		}

		postday95 = local->migra_min_heap.get_K().tid;
		if (postday95 == preday95) {
			memset(local->center_use, 0, sizeof(local->center_use));
			local->center_v95 = local->migra_min_heap.get_K().bandwidth;
			local->center_v95_t = postday95;
			// 标记后5%的时刻
			for (int tt = 0; tt < T; tt++) {
				if (local->migra_min_heap.mp[tt] != -1 && tt != postday95) local->center_use[tt] = 1;
			}
			break;
		}
		preday95 = postday95;
	}
}
