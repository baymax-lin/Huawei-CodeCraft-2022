#include "writer.h"
#include "solution.h"
#include "CodeCraft-2022.h"


char linebuf[1024 * 32];
void Writer::write_result(string path) {
	if (pState == nullptr) {
		perror("pState == nullptr");
		abort();
	}
	if (pState->T == 0) {
		perror("pState->T == 0");
		abort();
	}


	ofstream fout(path);
	streambuf* coutbackup = cout.rdbuf(fout.rdbuf());


	char* pAddr = linebuf;
	//char buf[100];
	// 先输出90点的十台
	//int count = 0;
	//bool _dot = false;
	//for (int i = 0; i < pState->edge_node_count; ++i) {
	//	if (pState->is_90[i] == 1) {
	//		count++;
	//		if (_dot) *pAddr++ = ',';
	//		memcpy(pAddr, pState->edge_id_str[i].c_str(), pState->edge_id_str[i].length());
	//		pAddr += pState->edge_id_str[i].length();
	//		_dot = true;
	//	}
	//}
	//if (count > 0) {
	//	*pAddr++ = '\n';
	//	cout << linebuf;
	//}
		
	
	
	for (int i = 0; i < pState->T; ++i) {

		// 决赛变更点的输出
		
		if (i > 0) {
			pAddr = linebuf;
			int cnt = 0;
			bool dot = false;
			for (int kk = 0; kk < pState->edge_node_count; ++kk) {
				if (pState->locals[pState->best_local_id]->select_20_edge[i][kk] == 0) continue;
				if (dot) *pAddr++ = ',';
				cnt++;
				int eid = kk;

				memcpy(pAddr, pState->edge_id_str[eid].c_str(), pState->edge_id_str[eid].length());
				pAddr += pState->edge_id_str[eid].length();

				dot = true;
			}
			if (cnt != 20) {
				LOG("wrong cnt != 20  ");
			}
			*pAddr++ = '\n';
			*pAddr++ = 0;
			cout << linebuf;
		}

		for (int j = 0; j < pState->client_node_count; ++j) {
			pAddr = linebuf;
			memcpy(pAddr, pState->client_id_str[j].c_str(), pState->client_id_str[j].length());
			pAddr += pState->client_id_str[j].length();
			*pAddr++ = ':';
			bool dot = false;
			//if (j == 1) {
			//	int bb = 0;
			//}
			// pState->alloc_mat 修订.
			for (int k = 0; k < pState->stream_types_num[i]; ++k) {
				if (pState->locals[pState->best_local_id]->alloc_mat[i][j][k] != NOT_ALLOCATION) {
					if (dot) *pAddr++ = ',';
					*pAddr++ = '<';
					int eid = pState->locals[pState->best_local_id]->alloc_mat[i][j][k];  // 该段流量分配去向
					memcpy(pAddr, pState->edge_id_str[eid].c_str(), pState->edge_id_str[eid].length());
					pAddr += pState->edge_id_str[eid].length();
					*pAddr++ = ',';

					//string str = to_string(pState->alloc_mat[i][k][j]);
					string str = pState->stream_types_name[i][k];
					memcpy(pAddr, str.c_str(), str.length());
					pAddr += str.length();
					*pAddr++ = '>';
					dot = true;
				}
				else {
					// 若为NOT_ALLOCATION, 说明该节点未分配或本身就没有流量?
				}
			}
			// 手动补换行符
			if ( ! (i == pState->T - 1 && j == pState->client_node_count - 1) ) *pAddr++ = '\n';
			*pAddr = 0;
			cout << linebuf;
		}
	}

	cout.rdbuf(coutbackup);
}