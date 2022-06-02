#include "CodeCraft-2022.h"
#include "reader.h"
#include <fstream>
#include <sstream>

int32 str2Int(string str) {
	int32 val = 0;
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == ' ') continue;
		if (str[i] < '0' || str[i] > '9') {
			perror("str[i] <'0' || str[i] > '9'");
			abort();
		}
			
		val = val * 10 + str[i] - '0';
	}
	return val;
}

void read_textstream(string path, string& result) {
	ifstream in(path);
	stringstream buffer;
	buffer << in.rdbuf();
	result = buffer.str();
	in.close();
	buffer.clear();
}

int split(string& s, int offset, vector<string>& tokens, char delimiters = ' ', char endChar = '\n')
{
	int32 lastPos = s.find_first_not_of(delimiters, offset);
	int32 pos = s.find(delimiters, lastPos);
	while (-1 != pos || -1 != lastPos) {
		tokens.emplace_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delimiters, pos);
		//pos = s.find(delimiters, lastPos);
		pos = lastPos;
		while (pos < s.length() && s[pos] != delimiters) {
			if (s[pos] == endChar) {
				tokens.emplace_back(s.substr(lastPos, pos - lastPos));
				return pos;
			}
			++pos;
		}
	}
	return pos;
}

int split_line(string& s, int offset, vector<string>& tokens, char delimiters = ' ')
{
	int32 lastPos = s.find_first_not_of(delimiters, offset);
	int32 pos = s.find(delimiters, lastPos);
	while (-1 != pos || -1 != lastPos) {
		tokens.emplace_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delimiters, pos);
		//pos = s.find(delimiters, lastPos);
		pos = lastPos;
		while (pos < s.length() && s[pos] != delimiters) {
			if (s[pos] == '\r' || s[pos] == '\n') {
				tokens.emplace_back(s.substr(lastPos, pos - lastPos));
				// 不同系统上可能换行符为\r\n. 
				return pos + 1 < s.length() && s[pos + 1] == '\n' ? pos + 1 : pos;
			}
			++pos;
		}
	}
	return pos;
}

void Reader::read_config(string path) {
	string buffer;
	if (pState == nullptr) {
		perror("pState == nullptr");
		abort();
	}
		
	read_textstream(path, buffer);
	const string keystr = "qos_constraint";
	
	int32 pos = buffer.find(keystr, 0) + keystr.length() + 1;
	if (pos == buffer.npos) {
		perror("pos == buffer.npos");
		abort();
	}
	int32 value = 0;
	for (int i = pos; i <buffer.size() && buffer[i] != '\n'; ++i) {
		if (buffer[i] == ' ' || buffer[i] == '\r') continue;
		if (buffer[i] < '0' || buffer[i] > '9') {
			perror("buffer[i] <'0' || buffer[i] > '9'");
			abort();
		}
		value = value * 10 + buffer[i] - '0';
	}
	pState->qos_limit = value;

	const string keystr2 = "base_cost";
	pos = buffer.find(keystr2, 0) + keystr2.length() + 1;

	value = 0;
	for (int i = pos; i < buffer.size() && buffer[i] != '\n'; ++i) {
		if (buffer[i] == ' ' || buffer[i] == '\r') continue;
		if (buffer[i] < '0' || buffer[i] > '9') {
			perror("buffer[i] <'0' || buffer[i] > '9'");
			abort();
		}
		value = value * 10 + buffer[i] - '0';
	}
	pState->base_cost = value;


	const string keystr3 = "center_cost";
	pos = buffer.find(keystr3, 0) + keystr3.length() + 1;

	int str_end = pos;
	while (str_end < buffer.size()) {
		if ((buffer[str_end] >= '0' && buffer[str_end] <= '9') || buffer[str_end] == '.') {
			str_end++;
		}
		else break;
	}
	string str_center_cost = buffer.substr(pos, str_end - pos);
	double v = stold(str_center_cost);
	pState->center_cost = v;
	//LOG(str_center_cost);


}

void Reader::read_demand(string path) {
	string buffer;
	if (pState == nullptr) {
		perror("pState == nullptr");
		abort();
	}

	read_textstream(path, buffer);

	int32 j = 0 ;

	vector<string> temp;
	j = split_line(buffer, 0, temp, ',');
	for (int k = 2; k < temp.size(); ++k) {     // 跳过两列
		pState->client_id_str[k - 2] = temp[k]; // 创建客户节点映射
		pState->client_str_id[temp[k]] = k - 2;
	}
	pState->client_node_count = temp.size() - 2; // 设置客户节点数量
	
	int32 t = 0;
	// read reamin lines
	vector<string> stream_name;
	string prev_time = "";
	int32 stream_id = -1;
	while (j + 1 < buffer.length()) {
		temp.clear();
		j = split_line(buffer, j + 1, temp, ',');
		if (temp[0] != prev_time) {
			if (t != 0) pState->stream_types_name.emplace_back(stream_name);
			if (t != 0) pState->stream_types_num[t - 1] = stream_id + 1;

			stream_name.clear();
			stream_id = 0;
			stream_name.emplace_back(temp[1]);

			++t;
			prev_time = temp[0];
		}
		else {
			stream_id++;
			stream_name.emplace_back(temp[1]);
		}
			
		for (int k = 2; k < temp.size(); ++k) {
			int32 val = str2Int(temp[k]);
			//if (t - 1) {
			//	int bbb = 0;
			//}
			pState->demands_all_types[t - 1][k - 2][stream_id].stream_id = stream_id;
			//pState->demands_all_types[t - 1][k - 2][stream_id].client_id = k - 2;
			pState->demands_all_types[t - 1][k - 2][stream_id].demand = val; // 设置某时刻某客户的需求流量
		}
	}
	// 尾部数据
	if (stream_name.size() > 0) {
		pState->stream_types_name.emplace_back(stream_name);
		pState->stream_types_num[t - 1] = stream_id + 1;
	}
	
	if (pState->T == 0) pState->T = t; // 设置时间长度
}

void Reader::read_site_bandWidth(string path) {
	string buffer;
	if (pState == nullptr) {
		perror("pState == nullptr");
		abort();
	}
	
	read_textstream(path, buffer);

	//LOG(buffer);
	int32 j = 0, c = 0;
	vector<string> temp;

	while (buffer[j] != '\n') ++j;
	// 读入边缘节点最大带宽并创建结构.
	while (j + 1 < buffer.length()) {
		temp.clear();
		j = split_line(buffer, j + 1, temp, ',');
		pState->edge_str_id[temp[0]] = c;
		pState->edge_id_str[c] = temp[0];
		int32 val = str2Int(temp[1]);
		pState->edge_nodes[c].max_bandwidth = val;
		//pState->edge_nodes[c].freeBandwidth = val;
		c++;
	}

	pState->edge_node_count = c;

}

void Reader::read_qos(string path) {
	string buffer;
	if (pState == nullptr) {
		perror("pState == nullptr");
		abort();
	}
	if (pState->edge_node_count == 0 || pState->client_node_count == 0) {
		perror("readDemand and readBandwidth should be called first");
		abort();
	}

	read_textstream(path, buffer);
	int32 j = 0, c = 0;
	vector<string> client_name;
	vector<string> temp;
	//while (buffer[j] != '\n') ++j;
	j = split_line(buffer, 0, client_name, ',');
	// 创建Qos矩阵. 假设顺序与 demand文件和site_bandwidth文件顺序一致，不作映射..
	while (j + 1 < buffer.length()) {
		temp.clear();
		j = split_line(buffer, j + 1, temp, ',');
		int32 edge_id = pState->edge_str_id[temp[0]];
		for (int k = 1; k < temp.size(); ++k) {
			int32 val = str2Int(temp[k]); 
			int32 client_id = pState->client_str_id[client_name[k]];
			pState->qos_mat[client_id][edge_id] = val; // k - 1: 客户节点序号; c : 边缘节点序号
			//pState->qos_mat[k-1][c] = val; // k - 1: 客户节点序号; c : 边缘节点序号
		}
		c++;
	}
}