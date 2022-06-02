#pragma once

#include <stdio.h>
#include <string>
#include "solution.h"
#include "CodeCraft-2022.h"


class Reader {

	Solution* pState;
public:
	void init_solution(Solution* _pState) {
		if (nullptr == _pState) {
			perror("nullptr == _pState");
			abort();
		}
		pState = _pState;
	}

	void read_config(string path);
	void read_demand(string path);
	void read_qos(string path);
	void read_site_bandWidth(string path);

};