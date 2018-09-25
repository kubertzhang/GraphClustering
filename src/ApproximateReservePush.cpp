#include <fstream>
#include <time.h>
#include <process.h>

#include "baseframework.h"
#include "MyGlobalParameters.h"

using namespace std;


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


void ApproximateReservePush::reservePush_Pre_Encode()
{
	static float * p = new float[g_vertexnum];
	static float * r = new float[g_vertexnum];
	static float * ur = new float[g_vertexnum * ATTRIBUTE_NUM];

	const int encode_buffer_size = 5 + (5 + ATTRIBUTE_NUM) * g_vertexnum;     // data block 的最大占用空间(包含最后结束判断位)
	float * encode_buffer = new float[encode_buffer_size];

	unsigned int write_buffer_size = 0;   // 实际 data block 大小
	float flag[1];                        // 标记该块后面是否还有数据块

	// 内存映射相关
	MMF * mmf_ptr = NULL;
	float * mmfm_base_address = NULL;
	float * write_dst_ptr;
	int mapfileid = 0;

	// 第一个mmf文件
	mmf_ptr = new MMF();
	char sharedname[1024] = "\0";
	char filename[1024] = "\0";
	sprintf_s(sharedname, "%d%s%d%s%d%s%d", g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", mapfileid);
	sprintf_s(filename, "%s%d%s%d%s%d%s%d", g_mmfPath.data(), g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", mapfileid, ".mmf");
	mmf_ptr->mmfCreate(sharedname, filename, a_mmfsizehigh, a_mmfsizelow, mmfm_base_address);
	write_dst_ptr = mmfm_base_address;
	mapfileid++;

	for (int targetID = g_cluster_startid; targetID <= g_cluster_endid; targetID++)
	{
		if (targetID % 1000 == 0)
			cout << targetID << endl;

		memset(p, 0, g_vertexnum * sizeof(float));
		memset(r, 0, g_vertexnum * sizeof(float));
		memset(ur, 0, g_vertexnum * ATTRIBUTE_NUM * sizeof(float));

		memset(encode_buffer, 0, encode_buffer_size * sizeof(float));

		// ---------- pushback ----------
		set<int> pushback_queue;

		r[targetID] = 1.0f;  // target point

		pushback_queue.insert(targetID);

		while (pushback_queue.size() > 0)
		{
			int uID = *pushback_queue.begin();
			pushback_queue.erase(pushback_queue.begin());

			Vertex * u = g_vertices[uID];    // 读取带pushback节点信息

			p[uID] += g_alpha * r[uID];      // estimated value

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

		// ---------- 计算ur ----------
		int u_ptr = 0;
		for (int i = 0; i < g_vertexnum; i++)
		{
			Vertex * u = g_vertices[i];
			int checkhop = u->edgetype_outdegree[STRUCTURE];       // 标记属性节点开始位置

			for (int j = ATTRIBUTE_1; j <= ATTRIBUTE_NUM; j++)   // 对应属性编号，从1开始
			{
				float sum_atribute_p = 0.0;
				for (int k = 0; k < u->edgetype_outdegree[j]; k++)
				{
					sum_atribute_p += p[u->neighborvertex[checkhop]];
					checkhop++;
				}

				float r_per_attribute;
				if (u->edgetype_outdegree[j] == 0)
					r_per_attribute = 0.0;
				else
					r_per_attribute = (1 - g_alpha) / (g_alpha * g_total_edgeweight * u->edgetype_outdegree[j]) * sum_atribute_p;

				ur[u_ptr] = r_per_attribute;
				u_ptr++;
			}
		}

		// ---------- 数据压缩并写入文件 ----------
		int buffer_index = 4;    // [0] - [3] 为数据块头 
		int b_p_count = 0;
		int b_r_count = 0;
		int b_ur_count = 0;

		// estimated value
		for (int i = 0; i < g_vertexnum; i++)
		{
			if (p[i] > 0)
			{
				b_p_count++;
				encode_buffer[buffer_index++] = (float)i;
				encode_buffer[buffer_index++] = p[i];
			}
		}
		// residual value
		for (int i = 0; i < g_vertexnum; i++)
		{
			if (r[i] > 0)
			{
				b_r_count++;
				encode_buffer[buffer_index++] = (float)i;
				encode_buffer[buffer_index++] = r[i];
			}
		}
		// update residual value
		for (int i = 0; i < g_vertexnum; i++)
		{
			bool c_flag = true;
			for (int tt = 0; tt < ATTRIBUTE_NUM; tt++)    // 三个属性对应值全为0的话，去掉
			{
				if (ur[i * ATTRIBUTE_NUM + tt] > 0)
				{
					c_flag = false;
					break;
				}
			}
			if (c_flag)
				continue;

			b_ur_count++;
			encode_buffer[buffer_index++] = (float)i;
			for (int j = i * ATTRIBUTE_NUM; j < i * ATTRIBUTE_NUM + ATTRIBUTE_NUM; j++)
			{
				encode_buffer[buffer_index++] = ur[j];
			}
		}

		unsigned int block_size = 4 + b_p_count * 2 + b_r_count * 2 + b_ur_count * (1 + ATTRIBUTE_NUM);   // 数据块大小

		// buffer块头部统计
		encode_buffer[0] = (float)targetID;     // 节点id
		encode_buffer[1] = (float)b_p_count;
		encode_buffer[2] = (float)b_r_count;
		encode_buffer[3] = (float)b_ur_count;

		write_buffer_size += (block_size + 1) * sizeof(float);   // 包含最后结束判断位

		if (write_buffer_size < a_mmf_buffer_size)
		{
			if (targetID > g_cluster_startid)  // 第一个数据块不需要在前面加标志位
			{
				flag[0] = 1.0f;     // 后面还有数据块
				memcpy(write_dst_ptr, flag, sizeof(float));
				write_dst_ptr = write_dst_ptr + 1;
			}

			// 数据写入内存映射文件
			memcpy(write_dst_ptr, encode_buffer, block_size * sizeof(float));
			write_dst_ptr = write_dst_ptr + block_size;
		}
		else     // 上一个文件大小已经达到上限
		{
			flag[0] = -1.0f;
			memcpy(write_dst_ptr, flag, sizeof(float));
			write_dst_ptr = write_dst_ptr + 1;

			mmf_ptr->mmfClose(mmfm_base_address);
			delete mmf_ptr;

			write_buffer_size = block_size * sizeof(float);    // 更新统计块

			// 创建一个新的文件
			mmf_ptr = new MMF();
			char sharedname[1024] = "\0";
			char filename[1024] = "\0";
			sprintf_s(sharedname, "%d%s%d%s%d%s%d", g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", mapfileid);
			sprintf_s(filename, "%s%d%s%d%s%d%s%d", g_mmfPath.data(), g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", mapfileid, ".mmf");
			mmf_ptr->mmfCreate(sharedname, filename, a_mmfsizehigh, a_mmfsizelow, mmfm_base_address);
			write_dst_ptr = mmfm_base_address;
			mapfileid++;

			// 数据写入内存映射文件
			memcpy(write_dst_ptr, encode_buffer, block_size * sizeof(float));
			write_dst_ptr = write_dst_ptr + block_size;
		}
	}
	// 关闭最后一个内存映射文件
	flag[0] = -1.0f;
	memcpy(write_dst_ptr, flag, sizeof(float));
	mmf_ptr->mmfClose(mmfm_base_address);
	delete mmf_ptr;

	// 保存预存文件的数量
	int pre_file_num = mapfileid;
	ofstream pre_file_ou;
	string pre_file_outputpath = g_mmfPath + to_string(g_datasetid) + "_" + to_string(g_clustertype) + "_" + to_string(g_schemeid) + ".num";
	pre_file_ou.open(pre_file_outputpath, ios::out);
	pre_file_ou << pre_file_num << endl;
	pre_file_ou.close();

	// 释放缓冲区
	delete p;
	delete r;
	delete ur;
	delete encode_buffer;
}


// 单线程，用于测试
void approximateReservePushSingleThreadUpdate(ofstream & log_ou)
{
	int taskid = 0;  // 人为指定文件名称编号

	int pr_id = 0;
	int ud_id = 0;
	int checkID = -1;

	g_mmfpool[taskid] = new MMF();
	char sharedname[1024] = "\0";
	char filename[1024] = "\0";
	sprintf_s(sharedname, "%d%s%d%s%d%s%d", g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", taskid);
	sprintf_s(filename, "%s%d%s%d%s%d%s%d", g_mmfPath.data(), g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", taskid, ".mmf");
	g_mmfpool[taskid]->mmfCreate(sharedname, filename, a_mmfsizehigh, a_mmfsizelow, g_mmfm_address_pool[taskid]);
	memset(g_bufferpool[taskid], 0, a_mmf_buffer_size);
	memcpy(g_bufferpool[taskid], g_mmfm_address_pool[taskid], a_mmf_buffer_size);
	g_mmfpool[taskid]->mmfClose(g_mmfm_address_pool[taskid]);
	delete g_mmfpool[taskid];
	// --------------------------------------------------------------------------------------

	float * buffersection = g_bufferpool[taskid];
	int b_index = 0;
	int b_p_count;
	int b_r_count;
	int b_ur_count;
	float continue_flag;

	int targetID;

	while (true)
	{
		memset(g_pr_pool[pr_id], 0, 2 * g_vertexnum * sizeof(float));   // 最后一位保存targetid

		targetID = (int)round(buffersection[b_index++]);                // _pr 最后一位保存 target_id
		b_p_count = (int)round(buffersection[b_index++]);
		b_r_count = (int)round(buffersection[b_index++]);
		b_ur_count = (int)round(buffersection[b_index++]);

		cout << targetID << endl;

		// decode estimate value
		for (int i = 0; i < b_p_count; i++)
		{
			int p_id = (int)round(buffersection[b_index++]);
			g_pr_pool[pr_id][p_id] = buffersection[b_index++];
		}

		// decode residual value
		for (int i = 0; i < b_r_count; i++)
		{
			int r_id = g_vertexnum + (int)round(buffersection[b_index++]);
			g_pr_pool[pr_id][r_id] = buffersection[b_index++];
		}

		// decode update residual value
		for (int i = 0; i < b_ur_count; i++)
		{
			int ur_id = g_vertexnum + (int)round(buffersection[b_index++]);
			for (int j = 1; j <= ATTRIBUTE_NUM; j++)
			{
				g_pr_pool[pr_id][ur_id] += (g_edgeweight[STRUCTURE][j] - g_base_weight[j]) * buffersection[b_index++];
			}
		}

		// 更新
		set<int> pushback_queue;

		for (int i = g_cluster_startid; i <= g_cluster_endid; i++)
		{
			if (g_pr_pool[ud_id][g_vertexnum + i] > g_epsilon)
			{
				pushback_queue.insert(i);
			}
		}

		// 根据更新结果再次进行 pushback
		while (pushback_queue.size() > 0)
		{
			//push_back u
			int uID = *pushback_queue.begin();
			pushback_queue.erase(pushback_queue.begin());

			Vertex * u = g_vertices[uID];                                            // 读取带pushback节点信息

			g_pr_pool[ud_id][uID] += g_alpha * g_pr_pool[ud_id][g_vertexnum + uID];      // estimated value

			//遍历u能够到达的点（reserve push）
			for (auto & wID : u->neighborvertex)
			{
				Vertex * w = g_vertices[wID];
				g_pr_pool[ud_id][g_vertexnum + wID] += (1 - g_alpha) * g_pr_pool[ud_id][g_vertexnum + uID] * getTransprob(w, u);    // residual value

				if (g_pr_pool[ud_id][g_vertexnum + wID] > g_epsilon)
				{
					pushback_queue.insert(wID);
				}
			}
			g_pr_pool[ud_id][g_vertexnum + uID] = 0;
		}
		pushback_queue.clear();

		// 获取邻居
		set<int> n_set;
		for (int i = 0; i < g_vertexnum; i++)
		{
			if (g_pr_pool[ud_id][i] > g_delta)
			{
				n_set.insert(i);
			}
		}

		g_dbscanneighborsets.insert(pair<int, set<int>>(targetID, n_set));
		n_set.clear();

		// 结束标志
		continue_flag = buffersection[b_index++];
		if (continue_flag < 0)
		{
			break;
		}
	}
}


unsigned int __stdcall approximateReservePushMultiThreadUpdateFunc(PVOID pM)
{
	int taskid = *(int*)pM;
	SetEvent(g_hThreadEvent);      //触发事件

	//互斥获取缓冲区
	EnterCriticalSection(&g_cs);

	if (g_buffer_queue.size() == 0)
	{
		cout << "Error: the buffer queue is empty!" << endl;
	}
	int buffer_id = *g_buffer_queue.begin();
	g_buffer_queue.erase(g_buffer_queue.begin());
	g_mmfpool[buffer_id] = new MMF();
	char sharedname[1024] = "\0";
	char filename[1024] = "\0";
	sprintf_s(sharedname, "%d%s%d%s%d%s%d", g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", taskid);
	sprintf_s(filename, "%s%d%s%d%s%d%s%d", g_mmfPath.data(), g_datasetid, "_", g_clustertype, "_", g_schemeid, "_", taskid, ".mmf");
	g_mmfpool[buffer_id]->mmfCreate(sharedname, filename, a_mmfsizehigh, a_mmfsizelow, g_mmfm_address_pool[buffer_id]);
	memset(g_bufferpool[buffer_id], 0, a_mmf_buffer_size);
	memcpy(g_bufferpool[buffer_id], g_mmfm_address_pool[buffer_id], a_mmf_buffer_size);
	g_mmfpool[buffer_id]->mmfClose(g_mmfm_address_pool[buffer_id]);
	delete g_mmfpool[buffer_id];
	LeaveCriticalSection(&g_cs);
	// --------------------------------------------------------------------------------------
	float * buffersection = g_bufferpool[buffer_id];
	int b_index = 0;
	int b_p_count;
	int b_r_count;
	int b_ur_count;
	float continue_flag;

	int pushback_c = 0;
	int targetID;

	vector<pair<int, set<int>>> set_queue;

	while (true)
	{
		memset(g_pr_pool[buffer_id], 0, 2 * g_vertexnum * sizeof(float));

		// 解析数据块
		targetID = (int)round(buffersection[b_index++]);
		b_p_count = (int)round(buffersection[b_index++]);
		b_r_count = (int)round(buffersection[b_index++]);
		b_ur_count = (int)round(buffersection[b_index++]);
		
		// decode estimate value
		for (int i = 0; i < b_p_count; i++)
		{
			int p_id = (int)round(buffersection[b_index++]);
			g_pr_pool[buffer_id][p_id] = buffersection[b_index++];
		}

		// decode residual value
		for (int i = 0; i < b_r_count; i++)
		{
			int r_id = g_vertexnum + (int)round(buffersection[b_index++]);
			g_pr_pool[buffer_id][r_id] = buffersection[b_index++];
		}

		// decode update residual value
		for (int i = 0; i < b_ur_count; i++)
		{
			int ur_id = g_vertexnum + (int)round(buffersection[b_index++]);
			for (int j = ATTRIBUTE_1; j <= ATTRIBUTE_NUM; j++)
			{
				g_pr_pool[buffer_id][ur_id] += (g_edgeweight[STRUCTURE][j] - g_base_weight[j]) * buffersection[b_index++];
			}
		}

		// 更新
		set<int> pushback_queue;

		for (int i = 0; i < g_vertexnum; i++)
		{
			if (g_pr_pool[buffer_id][g_vertexnum + i] > g_epsilon)
			{
				pushback_queue.insert(i);
			}
		}

		// 根据更新结果再次进行 pushback
		while (pushback_queue.size() > 0)
		{
			//push_back u
			int uID = *pushback_queue.begin();
			pushback_queue.erase(pushback_queue.begin());

			Vertex * u = g_vertices[uID];

			g_pr_pool[buffer_id][uID] += g_alpha * g_pr_pool[buffer_id][g_vertexnum + uID];      // estimated value
			pushback_c++;

			//遍历u能够到达的点（reserve push）
			for (auto & wID : u->neighborvertex)
			{
				Vertex * w = g_vertices[wID];
				g_pr_pool[buffer_id][g_vertexnum + wID] += (1 - g_alpha) * g_pr_pool[buffer_id][g_vertexnum + uID] * getTransprob(w, u);  // residual value

				if (g_pr_pool[buffer_id][g_vertexnum + wID] > g_epsilon)
				{
					pushback_queue.insert(wID);
				}
			}
			g_pr_pool[buffer_id][g_vertexnum + uID] = 0;
		}
		pushback_queue.clear();

		// 获取邻居
		set<int> n_set;
		for (int i = g_cluster_startid; i <= g_cluster_endid; i++)
		{
			if (g_pr_pool[buffer_id][i] > g_delta)
			{
				n_set.insert(i);
			}
		}
		set_queue.push_back(pair<int, set<int>>(targetID, n_set));
		
		// 结束标志
		continue_flag = buffersection[b_index++];
		if (continue_flag < 0)
		{
			break;
		}
	}

	//互斥释放缓冲区
	EnterCriticalSection(&g_cs);

	for (auto npair : set_queue)
		g_dbscanneighborsets.insert(npair);

	g_buffer_queue.insert(buffer_id);

	g_pushbackcount += pushback_c;
	LeaveCriticalSection(&g_cs);

	// 通知可加载新的线程
	ReleaseSemaphore(g_hSemaphoreRunnerNum, 1, NULL);
	return 0;
}


void ApproximateReservePush::approximateReservePushMultiThreadUpdate()
{
	g_dbscanneighborsets.clear();   // 重置

	InitializeCriticalSection(&g_cs);        
	g_hSemaphoreRunnerNum = CreateSemaphore(NULL, 7, 7, NULL);      // 参数表示允许同时运行的线程数(与缓冲区数目对应)
	g_hThreadEvent = CreateEvent(NULL, FALSE, FALSE, NULL);         // 用来保证参数传递的同步

	// 指定线程数量(任务数量)
	HANDLE * hThread = new HANDLE[a_THREADNUM];

	// 执行线程
	for (int i = 0; i < a_THREADNUM; i++)
	{
		// 等待有空的缓冲区出现
		WaitForSingleObject(g_hSemaphoreRunnerNum, INFINITE);
		hThread[i] = (HANDLE)_beginthreadex(NULL, 0, approximateReservePushMultiThreadUpdateFunc, &i, 0, NULL);
		WaitForSingleObject(g_hThreadEvent, INFINITE);             // 等待事件被触发  
	}

	WaitForMultipleObjects(a_THREADNUM, hThread, TRUE, INFINITE);

	for (int i = 0; i < a_THREADNUM; i++)
		CloseHandle(hThread[i]);
	delete hThread;
	CloseHandle(g_hSemaphoreRunnerNum);
	DeleteCriticalSection(&g_cs);
}


void ApproximateReservePush::execute()
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

	log_ou << endl << "********************************" << endl << "Approximate Approach: " << endl;
	// ====================================== 读图 ======================================
	readGraph();
	g_vertexnum = (int)g_vertices.size();
	// ====================================== 预处理 ====================================== 
	if (g_preflag)
	{
		reservePush_Pre_Encode();   // 压缩
		cout << "PreProcess is Completed!" << endl;
		log_ou << "PreProcess is Completed!" << endl;
		return;
	}
	// 读取预存文件的数量
	ifstream pre_file_in;
	string pre_file_inputpath = g_mmfPath + to_string(g_datasetid) + "_" + to_string(g_clustertype) + "_" + to_string(g_schemeid) + ".num";
	pre_file_in.open(pre_file_inputpath, ios::in);
	pre_file_in >> a_THREADNUM;
	pre_file_in.close();
	// ====================================== 迭代计算 ======================================
	// 分配缓冲区空间
	g_buffer_queue.clear();
	for (int i = 0; i < g_buff_size; i++)
	{
		g_bufferpool[i] = new float[a_mmf_buffer_size / sizeof(float)];
		g_pr_pool[i] = new float[2 * g_vertexnum];
		g_buffer_queue.insert(i);
	}

	// 运行时间重复20次取平均值
	int runTimes = 1;
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
			if (iterTimes > 30)
			{
				log_ou << "Can not converge!" << endl;
				return;
			}

			// Compute ppr score
			g_pushbackcount = 0;
			// 1. 单线程
			//approximateReservePushSingleThreadUpdate(log_ou);
			// 2. 多线程
			approximateReservePushMultiThreadUpdate();
			total_pushback_times += g_pushbackcount;

			// Symmetrization
			symmetrization();

			// DBSCAN
			dbscan();

			// Weight update
			log_ou << "Current weight: ";
			for (int i = ATTRIBUTE_1; i <= ATTRIBUTE_NUM; i++)
			{
				log_ou << g_edgeweight[STRUCTURE][i] << "\t";
			}
			log_ou << endl;
			diff = weightUpdate_Entropy();         // 熵（只实现按主类判断）
			//diff = weightUpdate_Vote();          // 投票（只实现按主类判断）
		}
		total_end = clock();

		total_running_time += total_end - total_start;
	}

	// 回收缓冲区
	for (int i = 0; i < g_buff_size; i++)
	{
		delete g_bufferpool[i];
		delete g_pr_pool[i];
	}

	// ====================================== 统计结果 ======================================
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

	log_ou.close();

	// 存储聚类结果
	storeClusterResult(cluster_output);
}
