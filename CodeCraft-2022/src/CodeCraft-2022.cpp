#include <iostream>
#include "CodeCraft-2022.h"
#include "writer.h"
#include "reader.h"
#include "solution.h"

#if __ONLINE__ == 1
 string configPath = "/data/config.ini";
 string demandPath = "/data/demand.csv";
 string siteBandwidthPath = "/data/site_bandwidth.csv";
 string qosPath = "/data/qos.csv";
 string resultPath = "/output/solution.txt";
#else
 string configPath = "../data/config.ini";
 string demandPath = "../data/demand.csv";
 string siteBandwidthPath = "../data/site_bandwidth.csv";
 string qosPath = "../data/qos.csv";
 string resultPath = "../data/solution.txt";
#endif


Writer wr;
Reader rd;
Solution *worker = new Solution();

int main() {
	ios::sync_with_stdio(false);

	rd.init_solution(worker);
	rd.read_config(configPath);
	rd.read_demand(demandPath);
	rd.read_site_bandWidth(siteBandwidthPath);
	rd.read_qos(qosPath);

	worker->init();
	worker->solve();

	worker->evaluate();

	wr.init_solution(worker);
	wr.write_result(resultPath);
	
    //std::cout << "Hello world!"<<std::endl;

	delete worker;

	return 0;
}