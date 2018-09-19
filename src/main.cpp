//#include<Python.h>
#include "baseframework.h"
#include "MyGlobalParameters.h"

using namespace std;


// 原始main函数
// ===================================================================================================================================
#define MAIN
#ifdef MAIN
int main(int argc, const char * argv[])
{
	/*
		命令行参数说明：
		argv[1] : m_inputpath   数据集文件名 {dblp_8w.txt, flickr_9w.txt 等}
		argv[2] : g_epsilon     ppr 误差值 {0.001}
		argv[3] : g_delta       dbscan 阈值 {0.005}
		argv[4] : m_minPts      dbscan minPts {4}
		argv[5] : dataset_id    数据集类型 {1 = dblp_8w, 2 = flickr_9w, 3 = dblp_3k, 4 = flickr_4k, 5 = case study - dblp_cs, 6 = football}
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
		BaseReversePush brp(argv);
		brp.execute();
	}
		   break;
	case 2:{
		// approximate
		ApproximateReversePush arp(argv);
		arp.execute();
	}
		   break;
	case 3:{
		// partial
		PartialReversePush p_pr(argv);
		p_pr.execute();
	}
		   break;
	case 4:{
		// game theory
		GameTheoty gt(argv);
		gt.execute();
	}
		   break;
	default: 
		cout << "Get the wrong scheme id !" << endl;
		break;
	}
	return 0;
}
#endif


// 选择符合要求的参数
// ===================================================================================================================================
//#define PRA
#ifdef PRA

int main(int argc, const char * argv[])
{
	/*
	命令行参数说明：
	argv[1] : m_inputpath   数据集文件名 {dblp_8w.txt, flickr_9w.txt 等}
	argv[2] : g_epsilon     ppr 误差值 {0.001}
	argv[3] : g_delta       dbscan 阈值 {0.005}
	argv[4] : m_minPts      dbscan minPts {4}
	argv[5] : dataset_id    数据集类型 {1 = dblp_8w, 2 = flickr_9w, 3 = dblp_3k, 4 = flickr_4k, 5 = case study - dblp_cs, 6 = football}
	argv[6] : schemeid      算法方案编号 {1 = base, 2 = approximate, 3 = partial, 4 = game theory}
	argv[7] ：pre_flag      是否进行预处理 {0 = no, 1 = yes}
	argv[8] : g_gamma		博弈论代价函数参数 {0.5}
	*/

	ofstream nmi_ou;
	nmi_ou.open("F:\\WorkSpace\\GraphClustering\\GC_ApproximateReversePush\\x64\\Release\\football_nmi_result.txt", ios::app);
	Py_Initialize();      // 初始化python解释器,告诉编译器要用的python编译器
	PyRun_SimpleString("import Py_NMI");         // 调用python文件

	// 设定参数区间
	// ===========
	// 需要在baseframework.cpp中注释语句
	// ===========

	// ApproximateReversePush.exe football.txt 0.001 0.005 4 6 1 0 0.5
	g_epsilon = 0.001f;

	int ptsNum = 6;
	int array_minPts[6] = { 3, 4, 5, 6, 7, 8};
	//int array_minPts[1] = { 4 };

	int scheNum = 1;
	int array_schemeid[2] = { 4 };
	//int array_schemeid[2] = { 1, 4 };

	float delta_start = 0.001f;
	float delta_end = 0.008f;
	float delta_rate = 0.0002f;

	for (float delta_cur = delta_start; delta_cur <= delta_end; delta_cur += delta_rate)
	{
		// 1. 设置两个参数
		g_delta = delta_cur;  // 参数1

		for (int i = 0; i < ptsNum; i++)
		{
			g_minPts = array_minPts[i];  // 参数2

			cout << "test : " << "delta = " << g_delta << "\t" << "minPts = " << g_minPts << endl;

			// 2. 计算两个方案的结果
			for (int j = 0; j < scheNum; j++)
			{
				g_schemeid = array_schemeid[j];  // 选择指定方案

				// 方案选择执行
				switch (g_schemeid)
				{
				case 1:{
					// basic 
					BaseReversePush brp(argv);
					brp.execute();
				}
					   break;
				case 2:{
					// approximate
					ApproximateReversePush arp(argv);
					arp.execute();
				}
					   break;
				case 3:{
					// partial
					PartialReversePush p_pr(argv);
					p_pr.execute();
				}
					   break;
				case 4:{
					// game theory
					GameTheoty gt(argv);
					gt.execute();
				}
					   break;
				default:
					cout << "Get the wrong scheme id !" << endl;
					break;
				}
			}

			// 3. 效果评价
			PyRun_SimpleString("Py_NMI.printNMI()");       // 调用python文件中的函数
			
			// 4. 结果分析
			ifstream nmi_in;
			string nmi_r_path = "nmi_result";
			nmi_in.open(nmi_r_path, ios::in);
			int nmi_res;
			float base;
			float game;
			nmi_in >> nmi_res;
			nmi_in >> base;
			nmi_in >> game;
			nmi_in.close();

			if (nmi_res)
			{
				nmi_ou << "YES" << "\tdelta = " << g_delta << "\tminPts = " << g_minPts << "\tbase = " << base <<  "\tgame = " << game << endl;
			}
			else
			{
				nmi_ou << "NO" << "\tdelta = " << g_delta << "\tminPts = " << g_minPts << "\tbase = " << base << "\tgame = " << game << endl;
			}
		}
	}
	nmi_ou.close();
	Py_Finalize();       // 结束python解释器，释放资源
	system("pause");
	return 0;
}
#endif