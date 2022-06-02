#pragma once

#include "CodeCraft-2022.h"
#include "solution.h"



class Writer {
	Solution* pState;
public:
	void init_solution(Solution* _pState) {
		if (nullptr == _pState) {
			perror("nullptr == _pState");
			abort();
		}
		pState = _pState;
	}
	void write_result(string path);
	
};