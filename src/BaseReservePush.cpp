#include <fstream>
#include <time.h>

#include "baseframework.h"
#include "MyGlobalParameters.h"

using namespace std;

// 计算转移概率: 可处理非完整图
inline float getTransprob(Vertex * _u, Vertex * _v)
{
	float real_total_edgeweight = 0.0f;
	for (auto n : _u->edgetype_outdegree)
	{
		if (n != 0)
			real_total_edgeweight += 1.0f;
	}
	return (g_edgeweight[_u->vertextype][_v->vertextype] / (_u->edgetype_outdegree[_v->vertextype] * real_total_edgeweight));
}


void BaseReservePush::baseReservePush()
{
	g_dbscanneighborsets.clear();  // 初始化

	for (int targetID = g_cluster_startid; targetID <= g_cluster_endid; targetID++)
	{
		vector<float> p(g_vertexnum, 0.0f);
		vector<float> r(g_vertexnum, 0.0f);

		set<int> pushback_queue;          // 存放带pushback的节点

		r[targetID] = 1.0f;               // target point
		pushback_queue.insert(targetID);

		while (pushback_queue.size() > 0)
		{
			int uID = *pushback_queue.begin();
			pushback_queue.erase(pushback_queue.begin());

			Vertex * u = g_vertices[uID];   // 读取带pushback节点信息

			p[uID] += g_alpha * r[uID];     // estimated value
			g_pushbackcount++;

			//遍历u能够到达的点（reserve push）
			for (auto & wID : u->neighborvertex)
			{
				Vertex * w = g_vertices[wID];
				r[wID] += (1 - g_alpha) * r[uID] * getTransprob(w, u);    // residual value

				if (r[wID] > g_epsilon)
				{
					pushback_queue.insert(wID);
				}
			}

			r[uID] = 0;
		}
		pushback_queue.clear();

		// 获取邻居
		set<int> n_set;
		for (int i = g_cluster_startid; i <= g_cluster_endid; i++)
		{
			if (p[i] > g_delta)
			{
				n_set.insert(i);
			}
		}
		g_dbscanneighborsets.insert(pair<int, set<int>>(targetID, n_set));
	}
}


void BaseReservePush::baseReservePushWithMemory()
{
	// 初始化
	m_pprDistances.clear();           
	g_dbscanneighborsets.clear();

	for (int targetID = g_cluster_startid; targetID <= g_cluster_endid; targetID++)
	{
		vector<float> p(g_vertexnum, 0.0);
		vector<float> r(g_vertexnum, 0.0);

		set<int> pushback_queue;      // 存放带pushback的节点

		r[targetID] = 1.0f;           // target point
		pushback_queue.insert(targetID);

		while (!pushback_queue.empty())
		{
			int uID = *pushback_queue.begin();
			pushback_queue.erase(pushback_queue.begin());

			Vertex * u = g_vertices[uID];   // 读取带pushback节点信息

			p[uID] += g_alpha * r[uID];     // estimated value
			g_pushbackcount++;

			//遍历u能够到达的点（reserve push）
			for (auto & wID : u->neighborvertex)
			{
				Vertex * w = g_vertices[wID];
				r[wID] += (1 - g_alpha) * r[uID] * getTransprob(w, u);    // residual value

				if (r[wID] > g_epsilon)
				{
					pushback_queue.insert(wID);
				}
			}
			r[uID] = 0.0;
		}

		// 保存距离
		// ===============================================
		map<int, float> m_map;
		for (int i = g_cluster_startid; i <= g_cluster_endid; i++)
		{
			if (p[i] > g_epsilon)     // ppr < g_epsilon 的点，存在误差，将其视为0
			{
				m_map.insert(pair<int, float>(i, p[i]));
			}
		}
		m_pprDistances.insert(pair<int, map<int, float>>(targetID, m_map));
		// ===============================================
	}
}


void BaseReservePush::execute()
{
	string result_output = g_resultpath + "result_" + to_string(g_datasetid) + "_" + to_string(g_delta) 
		+ "_" + to_string(m_minPts) + "_" + to_string(g_gamma) + "_" + to_string(g_epsilon) + ".txt";
	string cluster_output = g_resultpath + "cluster_result_" + to_string(g_datasetid) + "_" + to_string(g_schemeid) + "_"
		+ to_string(g_delta) + "_" + to_string(m_minPts) + "_" + to_string(g_gamma) + "_" + to_string(g_epsilon) + ".txt";

	ofstream log_ou;
	log_ou.open(result_output, ios::app);

	float diff = 1e10;
	int iterTimes = 0;
	long long total_pushback_times = 0;
	clock_t total_start, total_end;

	log_ou << endl << "********************************" << endl << "Base Approach: " << endl;
	// ====================================== 读图 ======================================
	readGraph();
	g_vertexnum = (int)g_vertices.size();

	// ====================================== 迭代计算 ======================================
	// 运行时间重复20次取平均值
	int runTimes = 1;  // default = 20
	long long total_running_time = 0;
	for (int i = 0; i < runTimes; i++)
	{
		// Initialize
		diff = 1e10;
		iterTimes = 0;
		total_pushback_times = 0;

		total_start = clock();
		while (diff > 1e-2)
		{
			iterTimes++;
			cout << "iterTimes = " << iterTimes << endl;
			if (iterTimes > 20)
			{
				log_ou << "Can not converge!" << endl;
				return;
			}

			// Compute ppr score
			g_pushbackcount = 0;
			baseReservePush();
			//baseReservePushWithMemory();
			total_pushback_times += g_pushbackcount;

			// Symmetrization
			symmetrization();
			//symmetrizationWithMemory();

			// DBSCAN
			dbscan();

			// Weight update
			log_ou << "Current weight: ";
			for (int i = ATTRIBUTE_1; i <= ATTRIBUTE_NUM; i++)
			{
				log_ou << g_edgeweight[STRUCTURE][i] << "\t";
			}
			log_ou << endl;
			diff = weightUpdate_Entropy();        // 熵（只实现按主类判断）
			//diff = weightUpdate_Vote();         // 投票（只实现按主类判断）
		}
		total_end = clock();

		total_running_time += total_end - total_start;
	}
	
	// ==================================== 统计结果 ===================================
	log_ou << "Total runningtime: " << total_running_time / runTimes << endl;
	log_ou << "Iteration Times: " << iterTimes << endl;
	log_ou << "Total Pushback Times: " << total_pushback_times << endl;
	log_ou << "Current weight: ";
	for (int i = ATTRIBUTE_1; i <= ATTRIBUTE_NUM; i++)
	{
		log_ou << g_edgeweight[STRUCTURE][i] << "\t";
	}
	log_ou << endl;
	
	// ----- 聚类分析 -----
	// 聚类数目
	log_ou << "Cluster_Amount: " << m_clusters.size() << endl;
	log_ou << "Valid_Cluster_points: " << m_valid_cluster_points.size() << endl;
	// density
	log_ou << "Cluster_Density: " << clusterEvaluate_Density() << endl;
	// entropy
	log_ou << "Cluster_Entropy: " << clusterEvaluate_NormalEntropy() << endl;

	// SToC
	// --------------------------------------------------------
	//string stoc_cluster_result_path = "F:\\WorkSpace\\GraphClustering\\GC_Result\\GT_vs_SToC\\stoc_" + m_stoc_file + "_" + to_string(m_minPts) + "_" + to_string(g_delta);
	//cout << "stoc_cluster_result_path = " << stoc_cluster_result_path << endl;
	//getClusterResult(stoc_cluster_result_path);
	//log_ou << "-----------------------------------------------------" << endl;
	//log_ou << "stoc_cluster_result_path = " << stoc_cluster_result_path << endl
	//       << "length = " << m_minPts << endl
	//	   << "tau = " << g_delta << endl  
	//	   << "cluster_amount = " << m_clusters.size() << endl
	//	   << "valid_cluster_points = " << m_valid_cluster_points.size() << endl
	//	   << "density = " << clusterEvaluate_Density() << endl
	//	   << "entropy = " << clusterEvaluate_NormalEntropy() << endl;
	//log_ou << "-----------------------------------------------------" << endl;
	// --------------------------------------------------------
	log_ou.close();

	//// 存储聚类结果
	//storeClusterResult(cluster_output);
	//// 存储聚类结果 == 对比SToC
	//string cluster_output_path = "F:\\WorkSpace\\GraphClustering\\GC_Result\\GT_vs_SToC\\cluster_result_1";
	//storeClusterResultForComparison(cluster_output_path);
}
