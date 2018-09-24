//#include<Python.h>
#include "baseframework.h"
#include "MyGlobalParameters.h"

using namespace std;

int main(int argc, const char * argv[])
{
	/*
		命令行参数说明：
		argv[1] : m_inputpath   数据集文件名 {dblp_8w.txt, flickr_9w.txt 等}
		argv[2] : g_epsilon     ppr 误差值 {0.001}
		argv[3] : g_delta       dbscan 阈值 {0.005}
		argv[4] : m_minPts      dbscan minPts {4}
		argv[5] : dataset_id    数据集类型 {1 = dblp_8w, 2 = flickr_9w, 3 = dblp_3k, 4 = flickr_4k, 5 = dblp_cs, 6 = football}
		argv[6] : schemeid      算法方案编号 {1 = base, 2 = approximate, 3 = partial, 4 = game theory}
		argv[7] ：pre_flag      是否进行预处理 {0 = no, 1 = yes}
		argv[8] : g_gamma		博弈论代价函数参数 {0.5}
	*/

	// 选择指定方案
	g_schemeid = atoi(argv[6]);
	switch (g_schemeid)
	{
	case 1:{
		// basic 
		BaseReservePush brp(argv);
		brp.execute();
	}
		   break;
	case 2:{
		// approximate
		ApproximateReservePush arp(argv);
		arp.execute();
	}
		   break;
	case 3:{
		// partial
		PartialReservePush p_pr(argv);
		p_pr.execute();
	}
		   break;
	case 4:{
		// game theory
		GameTheory gt(argv);
		gt.execute();
	}
		   break;
	default: 
		cout << "Get the wrong scheme id !" << endl;
		break;
	}
	return 0;
}
