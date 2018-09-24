#include <fstream>
#include <time.h>

#include "baseframework.h"
#include "MyGlobalParameters.h"


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


inline float GameTheory::getResponsecost(int _vertexid, int _clusterid)
{
	// 1. assignment cost = ppr; assignment cost only
	// -----------------------------------------
	//float ac = 0.0f;
	//float sc = 1.0f;
	//if (m_GlobalTable[_vertexid][3 + 3 * _clusterid] > 0)
	//	ac = m_cn * g_gamma * m_GlobalTable[_vertexid][3 + 3 * _clusterid + 1] / m_GlobalTable[_vertexid][3 + 3 * _clusterid];
	//float cost_c = sc - ac;

	// 2. assignment cost = ppr; assignment cost & social cost
	// -----------------------------------------
	float ac = 0.0f;
	if (m_GlobalTable[_vertexid][3 + 3 * _clusterid] > 0)
		ac = m_cn * g_gamma * m_GlobalTable[_vertexid][3 + 3 * _clusterid + 1] / m_GlobalTable[_vertexid][3 + 3 * _clusterid];
	float sc = m_GlobalTable[_vertexid][3 + 3 * _clusterid + 2];
	float cost_c = sc - ac;

	// 3. assignment cost = entropy; assignment cost & social cost
	// -----------------------------------------
	//float sc = m_GlobalTable[_vertexid][3 + 2 * _clusterid + 1];
	//float ac = m_cn * g_gamma * m_GlobalTable[_vertexid][3 + 2 * _clusterid];
	//float cost_c = sc + ac;
	
	return cost_c;
}


inline float GameTheory::getClusterEntropy(set<int> & _cluster)
{
	float entropy = 0.0f;
	for (int i = ATTRIBUTE_1; i < ATTRIBUTE_1 + ATTRIBUTE_NUM; i++)
	{
		entropy += getEntropy(_cluster, STRUCTURE, i);
	}
	return (entropy / ATTRIBUTE_NUM);
}


inline float GameTheory::AET_getClusterEntropy(int _vertexid, int _clusterid)
{
	float total_entropy = 0.0f;
	for (int i = ATTRIBUTE_1; i < ATTRIBUTE_1 + ATTRIBUTE_NUM; i++)
	{
		int total = 0;    // 统计连边总数
		for (auto & pair : AET[_vertexid][_clusterid][i - 1])   // 此处还可以预存
		{
			total += pair.second;
		}
		 
		float entropy = 0.0f;
		for (auto & pair : AET[_vertexid][_clusterid][i - 1])
		{
			float prob = pair.second * 1.0f / total;
			entropy += -1 * prob * log2(prob);
		}

		total_entropy += entropy;
	}
	return (total_entropy / ATTRIBUTE_NUM);
}


void GameTheory::gameTheory_ReservePush()
{
	m_pprDistances.clear();  // 初始化

	for (int targetID = g_cluster_startid; targetID <= g_cluster_endid; targetID++)
	{
		vector<float> p(g_vertexnum, 0.0);
		vector<float> r(g_vertexnum, 0.0);

		set<int> pushback_queue;      // 存放带pushback的节点

		r[targetID] = 1.0f;   // target point
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

		// 获取距离
		map<int, float> m_map;
		for (int i = g_cluster_startid; i <= g_cluster_endid; i++)
		{
			if (p[i] > g_epsilon)
			{
				m_map.insert(pair<int, float>(i, p[i]));
			}
		}
		m_pprDistances.insert(pair<int, map<int, float>>(targetID, m_map));
	}
}


void GameTheory::buildGlobalTable()
{
	// == Initialize
	m_GlobalTable.clear();  
	swap(m_happy_queue, stack<int>());  // 清空队列  
	cost_queue.clear();   

	// == Build glable table
	for (auto v_id : m_valid_cluster_points)
	{
		m_GlobalTable.insert(pair<int, vector<float>>(v_id, vector<float>(3 + 3 * m_clusters.size(), 0)));
	}

	for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
	{
		for (auto vertexid : m_clusters[clusterid])  
		{
			// 初始化表头
			m_GlobalTable[vertexid][0] = (float)clusterid;      // 实际类别                 

			// 计算maxSC
			int valid_neighbor_num = 0;
			Vertex * v = g_vertices[vertexid];
			for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
			{
				int f = v->neighborvertex[i];     // f is a friend of v
				if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())    // 邻居点中只选择 valid point
				{
					valid_neighbor_num++;
				}
			}
			m_GlobalTable[vertexid][2] = (float)valid_neighbor_num;    // 节点的度
			float maxSC = (1 - g_gamma) * 0.5f * valid_neighbor_num;   // 连边都是同类型的主类节点，并且权重为1，否则需要另行处理
			for (int c_id = 0; c_id < m_clusters.size(); c_id++)       // 对每个簇
			{
				m_GlobalTable[vertexid][3 + 3 * c_id + 2] = maxSC;
			}
		}
	}

	// == Initialize happy queue
	for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
	{
		for (auto vertexid : m_clusters[clusterid]) 
		{
			// social cost
			Vertex * v = g_vertices[vertexid];
			for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
			{
				int f = v->neighborvertex[i];  // f is friend of v
				if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())  // 邻居点中只选择 valid point
				{
					auto iter = m_GlobalTable.find(f);
					int f_clusterid = (int)round(iter->second[0]);   // 获取f所在的类别
					m_GlobalTable[vertexid][3 + 3 * f_clusterid + 2] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];  // 修改social cost
				}
			}

			// assignment cost
			for (auto & dis_map : m_pprDistances[vertexid])
			{
				int dis_id = dis_map.first;

				if (dis_id == vertexid)  // 距离集合中包含自身，去掉
					continue;

				if (m_valid_cluster_points.find(dis_id) != m_valid_cluster_points.end())
				{
					int dis_clusterid = (int)round(m_GlobalTable[dis_id][0]);               // 查找该点所在的类别
					m_GlobalTable[vertexid][3 + 3 * dis_clusterid] += 1;                    // num
					m_GlobalTable[vertexid][3 + 3 * dis_clusterid + 1] += dis_map.second;   // sum
				}
			}

			// 存储代价与类别
			set<CostNode> ss;
			for (int c_id = 0; c_id < m_clusters.size(); c_id++)
			{
				ss.insert(CostNode(c_id, getResponsecost(vertexid, c_id)));
			}

			cost_queue.insert(pair<int, set<CostNode>>(vertexid, ss));

			m_GlobalTable[vertexid][1] = (float)(ss.begin()->s_clusterid);       // 存储最小代价类别

			// 如果当前类别和最小代价类别不同，则需要进行调整， 将需要调整的点放入happy queue
			if ((int)round(m_GlobalTable[vertexid][0]) != (int)round(m_GlobalTable[vertexid][1]))
			{
				m_happy_queue.push(vertexid);
			}
		}
	}

	// == Compute normalization factor
	if (m_cn_flag)
	{
		int sum_degree = 0;
		float sum_max_ppr_dist = 0.0f;
		int ccount = 0;

		for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
		{
			for (auto vertexid : m_clusters[clusterid])  // 遍历所有节点
			{
				float max_ac = 1e8f;    // maximum ppr distance 
				for (int c_id = 0; c_id < m_clusters.size(); c_id++)  // 对每个簇
				{
					if (m_GlobalTable[vertexid][3 + 3 * c_id] > 0)
					{
						float aacc = m_GlobalTable[vertexid][3 + 3 * c_id + 1] / m_GlobalTable[vertexid][3 + 3 * c_id];
						max_ac = min(max_ac, aacc);
					}
				}

				if ((int)round(m_GlobalTable[vertexid][2]) != 0 && max_ac > 0)
				{
					ccount++;                                                 // size of valid points
					sum_degree += (int)round(m_GlobalTable[vertexid][2]);     // degree 
					sum_max_ppr_dist += max_ac;                               // ppr
				}
			}
		}

		float average_degree = (float)sum_degree / ccount;
		float average_weight = 1.0f;          // 连边都是同类型的主类节点，并且权重为1，否则需要另行处理       
		float average_max_ppr_dist = sum_max_ppr_dist / ccount;
		m_cn = average_degree * average_weight / (2.0f * average_max_ppr_dist * (float)sqrt(m_clusters.size()));
		m_cn_flag = false;
	}
	
}


void GameTheory::bestResponseDynamics()
{
	m_updatetimes = 0;    // 初始化

	while (!m_happy_queue.empty())
	{
		// 获取需要调整的节点
		int response_vertexid = m_happy_queue.top();
		m_happy_queue.pop();

		m_updatetimes++;  // 统计更新的点的次数

		// 进行类别调整与代价更新
		int currentclusterid = (int)round(m_GlobalTable[response_vertexid][0]);
		int nextclusterid = (int)round(m_GlobalTable[response_vertexid][1]);

		if (currentclusterid == nextclusterid)      // 在本类代价进一步减小的情况
			continue;

		m_GlobalTable[response_vertexid][0] = m_GlobalTable[response_vertexid][1];  // 调整类别

		/*所有的节点的代价在且仅在currentclusterid和nextclusterid发生变化*/

		// adjust social cost
		Vertex * v = g_vertices[response_vertexid];
		for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
		{
			int f = v->neighborvertex[i];  // f a is friend of v
			if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())  // f is a valid point
			{
				// old cluster
				if ((int)round(m_GlobalTable[f][0]) == currentclusterid)  
				{
					// 处理 response_vertexid 的 social cost
					m_GlobalTable[response_vertexid][3 + 3 * currentclusterid + 2] += (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];
					// 处理邻居节点的 social cost
					m_GlobalTable[f][3 + 3 * currentclusterid + 2] += (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];
				}

				// new cluster
				if ((int)round(m_GlobalTable[f][0]) == nextclusterid)
				{
					// 处理 response_vertexid 的 social cost
					m_GlobalTable[response_vertexid][3 + 3 * nextclusterid + 2] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];
					// 处理邻居节点的 social cost
					m_GlobalTable[f][3 + 3 * nextclusterid + 2] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];
				}	
			}
		}

		// adjust assignment cost
		for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
		{
			for (auto vertexid : m_clusters[clusterid]) 
			{
				float cur_cost = 1e8;
				float next_cost = 1e8;
				
				if (vertexid == response_vertexid)     // response_vertexid 自身的 assignment cost 不变
				{
					for (auto c_n : cost_queue[vertexid]) 
					{
						if (c_n.s_clusterid == currentclusterid)
							cur_cost = c_n.s_cost;
						if (c_n.s_clusterid == nextclusterid)
							next_cost = c_n.s_cost;

						if (cur_cost < 1e8 && next_cost < 1e8)
							break;
					}

					// old cluster
					cost_queue[vertexid].erase(cost_queue[vertexid].find(CostNode(currentclusterid, cur_cost)));                 // 删除旧的代价
					cost_queue[vertexid].insert(CostNode(currentclusterid, getResponsecost(vertexid, currentclusterid)));        // 插入新的代价
					// new cluster
					cost_queue[vertexid].erase(cost_queue[vertexid].find(CostNode(nextclusterid, next_cost)));                 // 删除旧的代价
					cost_queue[vertexid].insert(CostNode(nextclusterid, getResponsecost(vertexid, nextclusterid)));            // 插入新的代价

					m_GlobalTable[vertexid][1] = (float)(cost_queue[vertexid].begin())->s_clusterid;
					if ((int)round(m_GlobalTable[vertexid][0]) != (int)round(m_GlobalTable[vertexid][1]))  // 需要进行类别更新
					{
						m_happy_queue.push(vertexid);
					}

					continue;
				}
					
				/*
				response_vertexid调整, 其他vertexid的currentclusterid处的assignment cost可能发生变化
				如果vertexid和response_vertexid之间的ppr score > 0, 则会引起assignment cost的变化, 进行调整; 否则, 不进行调整
				*/
				
				if (m_pprDistances[vertexid].find(response_vertexid) != m_pprDistances[vertexid].end())   // vertexid 和 response_vertexid 之间的 ppr > 0
				{
					for (auto c_n : cost_queue[vertexid])
					{
						if (c_n.s_clusterid == currentclusterid)
							cur_cost = c_n.s_cost;
						if (c_n.s_clusterid == nextclusterid)
							next_cost = c_n.s_cost;

						if (cur_cost < 1e8 && next_cost < 1e8)
							break;
					}

					// old cluster
					m_GlobalTable[vertexid][3 + 3 * currentclusterid] -= 1;
					m_GlobalTable[vertexid][3 + 3 * currentclusterid + 1] -= m_pprDistances[vertexid][response_vertexid];

					cost_queue[vertexid].erase(cost_queue[vertexid].find(CostNode(currentclusterid, cur_cost)));                 // 删除旧的代价
					cost_queue[vertexid].insert(CostNode(currentclusterid, getResponsecost(vertexid, currentclusterid)));        // 插入新的代价

					// new cluster
					m_GlobalTable[vertexid][3 + 3 * nextclusterid] += 1;
					m_GlobalTable[vertexid][3 + 3 * nextclusterid + 1] += m_pprDistances[vertexid][response_vertexid];

					cost_queue[vertexid].erase(cost_queue[vertexid].find(CostNode(nextclusterid, next_cost)));                 // 删除旧的代价
					cost_queue[vertexid].insert(CostNode(nextclusterid, getResponsecost(vertexid, nextclusterid)));            // 插入新的代价
				}

				m_GlobalTable[vertexid][1] = (float)(cost_queue[vertexid].begin())->s_clusterid;
				if ((int)round(m_GlobalTable[vertexid][0]) != (int)round(m_GlobalTable[vertexid][1]))  // 需要进行类别更新
				{
					m_happy_queue.push(vertexid);
				}
			}
		}
	}
}


void GameTheory::E_buildGlobalTable()
{
	// 初始化
	m_GlobalTable.clear();
	AET.clear();
	for (auto v_id : m_valid_cluster_points)
	{
		m_GlobalTable.insert(pair<int, vector<float>>(v_id, vector<float>(3 + 2 * m_clusters.size(), 0)));
		AET.insert(pair<int, vector<vector<unordered_map<int, int>>>>(v_id, vector<vector<unordered_map<int, int>>>
			(m_clusters.size(), vector<unordered_map<int, int>>(ATTRIBUTE_NUM, unordered_map<int, int>()))));
	}

	for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
	{
		for (auto vertexid : m_clusters[clusterid])
		{
			// 初始化表头
			m_GlobalTable[vertexid][0] = (float)clusterid;      // 实际类别
			m_GlobalTable[vertexid][1] = -1.0f;                 // 最小代价类别
			m_GlobalTable[vertexid][2] = 1e8;                   // minCost = 无穷大

			// 计算maxSC
			int valid_neighbor_num = 0;
			Vertex * v = g_vertices[vertexid];
			for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
			{
				int f = v->neighborvertex[i];     // f is a friend of v
				if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())    // 邻居点中只选择 valid point
				{
					valid_neighbor_num++;
				}
			}
			float maxSC = (1 - g_gamma) * 0.5f * valid_neighbor_num;   // 连边都是同类型的主类节点，并且权重为1，否则需要另行处理

			for (int c_id = 0; c_id < m_clusters.size(); c_id++)       // 对每个簇
			{
				// 保存maxSC
				m_GlobalTable[vertexid][3 + 2 * c_id + 1] = maxSC;
				
				// 计算当前节点分配到各个簇的熵
				set<int> temp_cluster(m_clusters[c_id]);   // 复制cluster
				if (c_id != clusterid)
					temp_cluster.insert(vertexid);         // 待分配节点加入各个簇中

				// 计算熵
				//m_GlobalTable[vertexid][3 + 2 * c_id] = getClusterEntropy(temp_cluster);
			}
		}
	}
}


void GameTheory::AET_buildGlobalTable()
{
	// 初始化
	m_GlobalTable.clear();
	AET.clear();
	for (auto v_id : m_valid_cluster_points)
	{
		m_GlobalTable.insert(pair<int, vector<float>>(v_id, vector<float>(3 + 2 * m_clusters.size(), 0)));
		AET.insert(pair<int, vector<vector<unordered_map<int, int>>>>(v_id, vector<vector<unordered_map<int, int>>>
			(m_clusters.size(), vector<unordered_map<int, int>>(ATTRIBUTE_NUM, unordered_map<int, int>()))));
	}

	for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
	{
		for (auto vertexid : m_clusters[clusterid])
		{
			// 初始化表头
			m_GlobalTable[vertexid][0] = (float)clusterid;      // 实际类别
			m_GlobalTable[vertexid][1] = -1.0f;                 // 最小代价类别
			m_GlobalTable[vertexid][2] = 1e8;                   // minCost = 无穷大

			// 计算maxSC
			int valid_neighbor_num = 0;
			Vertex * v = g_vertices[vertexid];
			for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
			{
				int f = v->neighborvertex[i];     // f is a friend of v
				if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())    // 邻居点中只选择 valid point
				{
					valid_neighbor_num++;
				}
			}
			float maxSC = (1 - g_gamma) * 0.5f * valid_neighbor_num;   // 连边都是同类型的主类节点，并且权重为1，否则需要另行处理

			for (int c_id = 0; c_id < m_clusters.size(); c_id++)       // 对每个簇
			{
				// 保存maxSC
				m_GlobalTable[vertexid][3 + 2 * c_id + 1] = maxSC;

				// 计算当前节点分配到各个簇的熵
				set<int> temp_cluster(m_clusters[c_id]);   // 复制cluster
				if (c_id != clusterid)
					temp_cluster.insert(vertexid);         // 待分配节点加入各个簇中

				// 构建appearances表
				// ------------------------------------------------------------
				for (int i = ATTRIBUTE_1; i < ATTRIBUTE_1 + ATTRIBUTE_NUM; i++)
				{
					int type1 = STRUCTURE;
					int type2 = i;

					//枚举簇里面的每一个点作为源点
					for (int vid : temp_cluster)
					{
						Vertex * v = g_vertices[vid];

						if (v->vertextype == type1)             //如果该点属于type1
						{
							//枚举这个点能到达的其他点
							for (auto & uid : v->neighborvertex)
							{
								Vertex * u = g_vertices[uid];

								if (u->vertextype == type2)      //如果终点类型是type2
								{
									AET[vertexid][c_id][i - 1][u->vertexid] += 1;  // 统计appearances
								}
							}
						}
					}
				}
				// ------------------------------------------------------------

				// 计算熵
				m_GlobalTable[vertexid][3 + 2 * c_id] = AET_getClusterEntropy(vertexid, c_id);
			}
		}
	}
}


void GameTheory::E_initializeHappyQueue()
{
	swap(m_happy_queue, stack<int>());  // 清空队列  

	for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
	{
		for (auto vertexid : m_clusters[clusterid])  // 遍历所有节点
		{
			// 处理每个类别
			float minCost = 1e8;
			for (int c_id = 0; c_id < m_clusters.size(); c_id++)  // 对每个簇
			{
				// ** 计算代价
				float cost_c = getResponsecost(vertexid, c_id);

				if (cost_c < minCost)
				{
					minCost = cost_c;
					m_GlobalTable[vertexid][1] = (float)c_id;
				}
			}

			// 处理有效的邻居节点
			Vertex * v = g_vertices[vertexid];
			for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
			{
				int f = v->neighborvertex[i];  // f is friend of v
				if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())  // 邻居点中只选择 valid point
				{
					auto iter = m_GlobalTable.find(f);
					int f_clusterid = (int)round(iter->second[0]);   // 获取f所在的类别
					m_GlobalTable[vertexid][3 + 2 * f_clusterid + 1] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];  // 修改social cost

					// ** 计算代价
					float cost_c = getResponsecost(vertexid, f_clusterid);
					if (cost_c < minCost)
					{
						minCost = cost_c;
						m_GlobalTable[vertexid][1] = (float)f_clusterid;
					}
				}
			}

			// 保存 minCost 
			m_GlobalTable[vertexid][2] = minCost;

			// 如果当前类别和最小代价类别不同，则需要进行调整， 统计需要调整的点
			if (m_GlobalTable[vertexid][0] != m_GlobalTable[vertexid][1])
			{
				m_happy_queue.push(vertexid);
			}
		}
	}
}


void GameTheory::E_bestResponseDynamics()
{
	m_updatetimes = 0;
	while (!m_happy_queue.empty())
	{
		// 获取需要调整的节点
		int response_vertexid = m_happy_queue.top();
		m_happy_queue.pop();

		m_updatetimes++;  // 统计更新的点的次数

		// 进行类别调整与代价更新
		int currentclusterid = (int)round(m_GlobalTable[response_vertexid][0]);
		int nextclusterid = (int)round(m_GlobalTable[response_vertexid][1]);

		if (currentclusterid == nextclusterid)      // 在本类代价进一步减小的情况
			continue;

		m_GlobalTable[response_vertexid][0] = (float)nextclusterid;  // 调整类别

		// 从旧类中移除response_vertexid并加入新类
		m_clusters[currentclusterid].erase(m_clusters[currentclusterid].find(response_vertexid));
		m_clusters[nextclusterid].insert(response_vertexid);

		// 根据调整的变化调节其他所有点的 assignment cost
		for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
		{
			for (auto vertexid : m_clusters[clusterid])  // 遍历所有节点
			{
				if (vertexid == response_vertexid)  // vertexid自身的 assignment cost 不变
					continue;

				/*
					response_vertexid进行类别调整之后，所有其他点在currentclusterid和nextclusterid的熵都需要重新计算
				*/
				// 处理旧类, 计算熵之前需要移除response_vertexid
				set<int> temp_cluster(m_clusters[currentclusterid]);   // 复制cluster
				if (clusterid != currentclusterid)
				{
					temp_cluster.insert(vertexid);  // vertexid如果不在currentclusterid中，插入
				}
				m_GlobalTable[vertexid][3 + 2 * currentclusterid] = getClusterEntropy(temp_cluster);   // 更新熵

				// 处理新类, 计算熵之前需要添加response_vertexid
				set<int> temp_cluster2(m_clusters[nextclusterid]);   // 复制cluster
				if (clusterid != nextclusterid)
				{
					temp_cluster2.insert(vertexid);  // vertexid如果不在nextclusterid中，插入
				}
				m_GlobalTable[vertexid][3 + 2 * nextclusterid] = getClusterEntropy(temp_cluster2);   // 更新熵

				// ** 计算代价
				float cost_old = getResponsecost(vertexid, currentclusterid);
				float cost_new = getResponsecost(vertexid, nextclusterid);
				float cost_c = min(cost_old, cost_new);
				if (cost_c < m_GlobalTable[vertexid][2])
				{
					m_GlobalTable[vertexid][2] = cost_c;
					m_GlobalTable[vertexid][1] = (cost_old < cost_new) ? (float)currentclusterid : (float)nextclusterid;
					m_happy_queue.push(vertexid);
				}
			}
		}

		// 根据邻居节点调整 social cost
		Vertex * v = g_vertices[response_vertexid];
		for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
		{
			int f = v->neighborvertex[i];  // f is friend of v
			if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())  // 邻居点中只选择 valid point
			{
				// 旧类, 此处代价只会增大，不会有更小的解
				m_GlobalTable[f][3 + 2 * currentclusterid + 1] += (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];

				// 新类
				m_GlobalTable[f][3 + 2 * nextclusterid + 1] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];

				// ** 计算代价
				float cost_c = getResponsecost(f, nextclusterid);
				if (cost_c < m_GlobalTable[f][2])
				{
					m_GlobalTable[f][2] = cost_c;
					m_GlobalTable[f][1] = (float)nextclusterid;
					m_happy_queue.push(f);
				}
			}
		}
	}
}


void GameTheory::AET_bestResponseDynamics()
{
	m_updatetimes = 0;
	while (!m_happy_queue.empty())
	{
		// 获取需要调整的节点
		int response_vertexid = m_happy_queue.top();
		m_happy_queue.pop();

		m_updatetimes++;  // 统计更新的点的次数

		// 进行类别调整与代价更新
		int currentclusterid = (int)round(m_GlobalTable[response_vertexid][0]);
		int nextclusterid = (int)round(m_GlobalTable[response_vertexid][1]);

		if (currentclusterid == nextclusterid)      // 在本类代价进一步减小的情况
			continue;

		m_GlobalTable[response_vertexid][0] = (float)nextclusterid;  // 调整类别

		Vertex * v = g_vertices[response_vertexid];

		// 根据调整的变化调节其他所有点的 assignment cost
		for (int clusterid = 0; clusterid < m_clusters.size(); clusterid++)
		{
			for (auto vertexid : m_clusters[clusterid])  // 遍历所有节点
			{
				if (vertexid == response_vertexid)  // vertexid自身的 assignment cost 不变
					continue;

				/*
				更新: response_vertexid进行类别调整之后，所有其他点在currentclusterid和nextclusterid的熵都需要重新计算
				*/
				for (int ii = v->edgetype_outdegree[STRUCTURE]; ii < v->neighborvertex.size(); ii++)   // 处理属性邻居点
				{
					int uid = v->neighborvertex[ii];  // 取出该属性邻居点
					Vertex * u = g_vertices[uid];

					// 处理旧类, 在vertexid对应的currentclusterid的类别中删除response_vertexid产生的appearance
					if (AET[vertexid][currentclusterid][u->vertextype - 1][uid] < 2)  // 只有一次appearance, 删除这个键值对
						AET[vertexid][currentclusterid][u->vertextype - 1].erase(AET[vertexid][currentclusterid][u->vertextype - 1].find(uid));
					else
						AET[vertexid][currentclusterid][u->vertextype - 1][uid] -= 1;

					// 处理新类, 在vertexid对应的nextclusterid的类别中添加response_vertexid产生的appearance
					AET[vertexid][nextclusterid][u->vertextype - 1][uid] += 1;
				}

				m_GlobalTable[vertexid][3 + 2 * currentclusterid] = AET_getClusterEntropy(vertexid, currentclusterid);  // 旧类
				m_GlobalTable[vertexid][3 + 2 * nextclusterid] = AET_getClusterEntropy(vertexid, nextclusterid);        // 新类

				// ** 计算代价
				float cost_old = getResponsecost(vertexid, currentclusterid);
				float cost_new = getResponsecost(vertexid, nextclusterid);
				float cost_c = min(cost_old, cost_new);
				if (cost_c < m_GlobalTable[vertexid][2])
				{
					m_GlobalTable[vertexid][2] = cost_c;
					m_GlobalTable[vertexid][1] = (cost_old < cost_new) ? (float)currentclusterid : (float)nextclusterid;
					m_happy_queue.push(vertexid);
				}
			}
		}

		// 根据邻居节点调整 social cost
		for (int i = 0; i < v->edgetype_outdegree[STRUCTURE]; i++)
		{
			int f = v->neighborvertex[i];  // f is friend of v
			if (m_valid_cluster_points.find(f) != m_valid_cluster_points.end())  // 邻居点中只选择 valid point
			{
				// 旧类, 此处代价只会增大，不会有更小的解
				m_GlobalTable[f][3 + 2 * currentclusterid + 1] += (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];

				// 新类
				m_GlobalTable[f][3 + 2 * nextclusterid + 1] -= (1 - g_gamma) * 0.5f * g_edgeweight[STRUCTURE][STRUCTURE];
				// ** 计算代价
				float cost_c = getResponsecost(f, nextclusterid);
				if (cost_c < m_GlobalTable[f][2])
				{
					m_GlobalTable[f][2] = cost_c;
					m_GlobalTable[f][1] = (float)nextclusterid;
					m_happy_queue.push(f);
				}
			}
		}
	}
}


void GameTheory::gatherClusterResult()
{
	vector<set<int>> old_clusters(m_clusters);
	vector<set<int>> new_clusters(m_clusters.size());
	m_clusters.clear();

	for (int clusterid = 0; clusterid < old_clusters.size(); clusterid++)
	{
		for (auto vertexid : old_clusters[clusterid])  // 遍历所有节点
		{
			// 重新统计各个点的聚类结果
			new_clusters[(int)round(m_GlobalTable[vertexid][0])].insert(vertexid);
		}
	}

	// 对聚类后的结果进行筛选
	for (auto iter = new_clusters.begin(); iter != new_clusters.end();)
	{
		if ((*iter).empty())
		{
			iter = new_clusters.erase(iter);
		}
		else
			iter++;
	}

	m_clusters = new_clusters;
}


void GameTheory::gameTheoryModulation()
{
	// 1. assignment cost == ppr
	// --------------------------------------------------
	buildGlobalTable();        // 构建 Global Table, 计算 cN
	bestResponseDynamics();    // best-response dynamics

	// 2. assignment cost = entropy
	// --------------------------------------------------
	//E_buildGlobalTable();        // 构建 Global Table, 计算 cN
	//E_initializeHappyQueue();    // 初始化 happy_queue
	//E_bestResponseDynamics();    // best-response dynamics

	// 3. assignment cost = entropy(使用Appearance Entropy Table)
	// --------------------------------------------------
	//AET_buildGlobalTable();        // 构建 Global Table, 计算 cN
	//E_initializeHappyQueue();      // 初始化 happy_queue(与方案2相同)
	//AET_bestResponseDynamics();    // best-response dynamics(AET优化方案)

	// 重新统计聚类结果
	gatherClusterResult();
}


void GameTheory::execute()
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
	int total_Update_times = 0;
	clock_t total_start, total_end;

	log_ou << endl << "********************************" << endl << "GameTheory Approach: " << endl;

	// ====================================== 读图 ======================================
	readGraph();
	g_vertexnum = (int)g_vertices.size();

	// ====================================== 迭代计算 ======================================
	m_cn_flag = true;   // 只计算一次cn

	// 运行时间重复20次取平均值
	int runTimes = 1;  // default = 20
	long long total_running_time = 0;
	for (int i = 0; i < runTimes; i++)
	{
		// Initialize
		diff = 1e10;
		iterTimes = 0;
		total_pushback_times = 0;
		total_Update_times = 0;

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
			gameTheory_ReservePush();
			total_pushback_times += g_pushbackcount;

			// Symmetrization
			symmetrizationWithMemory();

			// DBSCAN
			dbscan();

			// Game Theory
			gameTheoryModulation();
			total_Update_times += m_updatetimes;

			// Weight update
			diff = weightUpdate_Entropy();        // 熵（只实现按主类判断）
			//diff = weightUpdate_Vote();         // 投票（只实现按主类判断）
		}
		total_end = clock();

		total_running_time += total_end - total_start;
	}

	// ====================================== 统计结果 ======================================
	log_ou << "Total runningtime: " << total_running_time / runTimes << endl;
	log_ou << "Iteration Times: " << iterTimes << endl;
	log_ou << "Total Pushback Times: " << total_pushback_times << endl;
	log_ou << "Total Update Times: " << total_Update_times << endl;
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
