#pragma once

#include "CodeCraft-2022.h"

#define MAX_HEAP_LEN (int(MAX_T_LEN * 0.05) + 100)

template<class type, class _cmp>
class Heap {
public:

	//type arr[MAX_T_LEN];
	type* arr;
	short mp[MAX_T_LEN];

	int n = 0;
	int K = 0;
	_cmp cmp;
	
	Heap(int size = MAX_HEAP_LEN) {
		n = 0;
		K = 0;
		memset(mp, -1, sizeof(mp));
		arr = new type[size];
	}
	~Heap() {
		delete[] arr;
	}

	void copy_from(const Heap & hp) {
		n = hp.n;
		K = hp.K;
		memcpy(arr, hp.arr, sizeof(arr));
		memcpy(mp, hp.mp, sizeof(mp));
	}

	void set_K(int _K) {
		K = _K;
		if (K == 0) K = 1;
	}
	void reset() {
		n = 0;
		memset(mp, -1, sizeof(mp));
	}
	bool empty() {
		return n == 0;
	}

	type top() {
		if (n == 0) return type(); // 若该节点没有被使用, 返回默认构造
		return arr[0];
	}

	bool fill_K() {
		return n >= K;
	}

	bool exist_but_not_top(int t) {
		if (n == 0) return false;
		return (mp[t] != -1 && t != arr[0].getkey());
	}

	type get_K() {
		if (n >= K) return arr[0];
		else return type();
	}
	type get_K_next() {
		if (n >= K)  return second();
		else if (n == K - 1) return arr[0];
		else return type();
	}

	type second() {
		int a = 1, b = 2;
		if (n <= 1) return type();
		else if (n == 2) return arr[1];
		if (cmp(arr[a], arr[b])) return arr[a];
		else return arr[b];
	}

	type get_by_key(int key) {
		if (mp[key] == -1) return type();
		return arr[mp[key]];
	}
	void push_topK(type x) {
		if (K == 0) return;
		int key = x.getkey();
		if (mp[key] == -1) {
			if (n < K) {
				arr[n++] = x;
				up(n - 1);
			}
			else if (cmp(arr[0], x)) { // arr[0] < x # arr[0].t > x.t -> 替换 
				mp[arr[0].getkey()] = -1;
				arr[0] = x;
				down(0);
			}
		}
		else {
			int i = mp[key];
			auto temp = arr[i];
			arr[i] = x;

			if (cmp(x, temp)) up(i);  // x < temp # x.t < temp.t
			else down(i);
		}
	}

	void push(type x) {
		int key = x.getkey();
		if (mp[key] == -1) {
			arr[n++] = x;
			up(n - 1);
		}
		else {
			int i = mp[key];
			auto temp = arr[i];
			arr[i] = x;

			if (cmp(x, temp)) up(i);
			else down(i);
		}

	}

	void pop() {
		if (n == 0) return;
		mp[arr[0].getkey()] = -1;
		if (n == 1) {
			n--;
			return;
		}
		swap(arr[0], arr[n - 1]);
		n--;
		down(0);
	}


	void up(int i) {
		type v = arr[i];
		int j = (i - 1) / 2;    // 父结点

		while (j >= 0 && i != 0) {
			if (cmp(arr[j], v)) // arr[j] < v
				break;
			arr[i] = arr[j];
			mp[arr[i].getkey()] = i;
			i = j;
			j = (i - 1) / 2;
		}
		arr[i] = v;
		mp[arr[i].getkey()] = i;
	}

	void down(int i) {
		int j = 2 * i + 1;     //左孩子
		type v = arr[i];
		while (j < n) {
			if (j + 1 < n && cmp(arr[j + 1], arr[j]))
				j++;           // 寻找左右孩子中较小的
			if (cmp(v, arr[j])) // v < arr[j]
				break;
			arr[i] = arr[j];
			mp[arr[i].getkey()] = i;
			i = j;
			j = 2 * i + 1;
		}
		arr[i] = v;
		mp[arr[i].getkey()] = i;
	}

};