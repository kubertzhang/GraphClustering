#ifndef BASEFRAMEWORK_H
#define BASEFRAMEWORK_H

#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <queue>
#include <stack>

#include "Vertex.h"
#include "mmf.h"

using namespace std;


class BaseFramework
{
protected:
	string m_inputpath;         // 输入数据路径
	string m_stoc_file;
	int m_minPts;               // dbscan 参数 minPts

	map<int, map<int, float>> m_pprDistances;		// 存储ppr score
	unordered_set<int> m_valid_cluster_points;		// 统计有效的聚类点
	vector<set<int>> m_clusters;					// 聚类结果

	void readGraph();					// 读图
	void symmetrization();			// 不保存距离的对称化方式, approximate, partial 方案使用
	void symmetrizationWithMemory();		// 保存距离的对称化方式, game theory 方案使用
	void dbscan();						// DBSCAN
	float getEntropy(set<int> & cluster, int type1, int type2);    // 计算特定类型(type1, type2)的熵
	float weightUpdate_Entropy();                                  // 基于熵的权重更新方法
	float weightUpdate_Vote();                                     // 基于投票的权重更新算法

	void storeClusterResult(string cluster_output);                // 存储聚类结果
	void storeClusterResultForComparison(string cluster_output_path);   // 在baseline的代码里计算评价指标

	void getClusterResult(std::string _cluster_result_path);

	float clusterEvaluate_Density();                               // 聚类效果评价 - density
	float clusterEvaluate_NormalEntropy();                              // 聚类效果评价 - entropy (直接求平均)
	float clusterEvaluate_WeightedEntropy();                              // 聚类效果评价 - entropy (加权平均)
	float clusterEvaluate_WithinClusterAveDistance();              // 聚类效果评价 - 簇间平均距离
	float clusterEvaluate_NMI();                                   // 聚类效果评价 - NMI

public:
	explicit BaseFramework(const char * argv[]);
	virtual ~BaseFramework(){}
};


// 基本方案 
class BaseReservePush : public BaseFramework
{
private:
	void baseReservePush();
	void baseReservePushWithMemory();
public:
	explicit BaseReservePush(const char * argv[]) : BaseFramework(argv){}
	~BaseReservePush(){}
	void execute();
};


// 近似计算方案
class ApproximateReservePush : public BaseFramework
{
public:
	explicit ApproximateReservePush(const char * argv[]) : BaseFramework(argv){}
	~ApproximateReservePush(){}
	void execute();
private:
	void reservePush_Pre_Encode();
	void approximateReservePushMultiThreadUpdate();
};


// 部分计算方案
class PartialReservePush : public BaseFramework
{
public:
	explicit PartialReservePush(const char * argv[]) : BaseFramework(argv){}
	~PartialReservePush(){}
	void execute();
private:
	void reservePush_Partial_Encode();
	void partialReservePushMultiThreadUpdate();
};


// 博弈论方案
class GameTheory : public BaseFramework
{
public:
	explicit GameTheory(const char * argv[]) : BaseFramework(argv){};
	~GameTheory(){}
	void execute();
private:
	struct DegreenNode    // 保存节点及其对应的度 
	{
		int s_vertexid;
		int s_degree;

		DegreenNode(int _vertexid, int _degree)
		{
			this->s_vertexid = _vertexid;
			this->s_degree = _degree;
		}

		bool operator < (const DegreenNode & _dn) const  // 降序
		{
			if (this->s_degree != _dn.s_degree)
				return (this->s_degree > _dn.s_degree);   
			if (this->s_vertexid != _dn.s_vertexid)
				return (this->s_vertexid < _dn.s_vertexid);
			return false;
		}
	};

	struct CostNode     // 保存节点及其对应的调整代价
	{
		int s_clusterid;
		float s_cost;

		CostNode(int _clusterid, float _cost)
		{
			this->s_clusterid = _clusterid;
			this->s_cost = _cost;
		}

		bool operator < (const CostNode & _cn) const  // 降序
		{
			if (this->s_cost != _cn.s_cost)
				return (this->s_cost < _cn.s_cost);             // s_cost 不同，按照 s_cost 升序排序
			if (this->s_clusterid != _cn.s_clusterid)     
				return this->s_clusterid < _cn.s_clusterid;     // s_cost 相同， 按照 s_clusterid 的升序排序
			return false;                                       // 两边不能同时为true: strict weak ordering
		}
	};

	
	unordered_map<int, vector<float>> m_GlobalTable;   // Global Table, 用于best response dynamics的优化

	unordered_map<int, set<CostNode>> cost_queue;      // priority queue for cost, ascending order of cost

	//set<DegreenNode> m_happy_queue;       // store candidate points - descending order of degree
	//unordered_set<int> m_happy_queue;     // store candidate points - no order
	//queue<int> m_happy_queue;             // store candidate points - first in first out
	stack<int> m_happy_queue;             // store candidate points - first in last out

	float m_cn;             // 归一化参数(暂不考虑)
	bool m_cn_flag;         // 是否计算并使用归一化参数
	int m_updatetimes;      // update points times per iteration
	

	unordered_map<int, vector<vector<unordered_map<int, int>>>> AET;   // 存储各属性点的出现的次数, 用来优化熵的计算
	/*            点     簇    属性    点出现次数                */

	float getResponsecost(int _vertexid, int _clusterid);
	void gameTheory_ReservePush();

	void gameTheory_dbscan();

	void buildGlobalTable();
	void initializeHappyQueue();
	void bestResponseDynamics();

	// （弃用）
	float getClusterEntropy(set<int> & _cluster);   // 计算簇内的平均熵
	void E_buildGlobalTable();
	void E_initializeHappyQueue();
	void E_bestResponseDynamics();

	// （弃用）
	float AET_getClusterEntropy(int _vertexid, int _clusterid);
	void AET_buildGlobalTable();
	void AET_bestResponseDynamics();
	
	void gatherClusterResult();       // 收集每轮迭代调整后的结果
	void gameTheoryModulation();      // 博弈论优化算法
};

#endif