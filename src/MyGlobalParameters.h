#ifndef MYGLOBALPARAMETERS_H
#define MYGLOBALPARAMETERS_H
#include <map>
#include <set>
#include <iostream>
#include <fstream>

#include "Vertex.h"
#include "mmf.h"

using namespace std;

                                 // 数据集       DBLP			Flickr		Football
#define STRUCTURE			0    // 主节点		 paper			image		player
#define ATTRIBUTE_1			1    // 属性点1		 author			user		list
#define ATTRIBUTE_2			2    // 属性点2		 conference		tag			tweet
#define ATTRIBUTE_3			3    // 属性点2		 keywords
extern int ATTRIBUTE_NUM;		 // 属性数目     3				2

#define ATTRIBUTE_BUFF_SIZE 4

// global parameters
extern int g_clustertype;                            // 聚类类型，默认为STRUCTURE
extern int g_datasetid;                              // 数据集编号
extern int g_schemeid;                               // 方案编号
extern float g_epsilon;                              // reverse push 误差参数
extern float g_delta;                                // dbscan range
extern int g_minPts;                                 // dbscan minPts
extern float g_alpha;                                // reverse push 重启概率
extern float g_gamma;                                // game theory 代价函数参数
extern int g_cluster_startid;                        // 聚类类型开始节点(包含该点)
extern int g_cluster_endid;                          // 聚类类型结束节点(包含该点)
extern int g_vertexnum;                              // 总的节点数目
extern long long g_pushbackcount;                    // pushback 次数
extern vector<Vertex*> g_vertices;                   // 点结构
extern map<int, set<int>> g_dbscanneighborsets;      // 存邻居节点(已指定节点类型)

const int g_buff_size = 7;                    // 缓冲区个数
extern float * g_bufferpool[];                // 预存数据缓冲区
extern MMF * g_mmfpool[];                     // MMF对象池
extern float * g_mmfm_address_pool[];         // MMF 地址池
extern set<int> g_buffer_queue;               // buffer池分配队列
extern float * g_pr_pool[];                   // pr缓冲区

extern HANDLE  g_hSemaphoreRunnerNum;      // 信号量，设定最大的并行线程数
extern HANDLE  g_hThreadEvent;             // 信号量，用于线程分布的同步
extern CRITICAL_SECTION g_cs;              // 关键段

// 边权重
extern float g_structure_edgeweight;              // 主类节点边权重 { = 1}
extern float g_total_edgeweight;                  // 所有类型边的权重和
extern float g_attribute_total_edgeweight;        // 属性节点边权重和
extern float g_base_weight[ATTRIBUTE_BUFF_SIZE];                                // 初始权重
extern float g_former_edgeweight[ATTRIBUTE_BUFF_SIZE][ATTRIBUTE_BUFF_SIZE];     // 在BaseFramework::BaseFramework()中初始化
extern float g_edgeweight[ATTRIBUTE_BUFF_SIZE][ATTRIBUTE_BUFF_SIZE];            // 在BaseFramework::BaseFramework()中初始化

extern string g_inputpath;        // 输入数据文件绝对路径
extern string g_resultpath;       // 输出结果文件绝对路径

extern string g_mmfPath;          // 预存数据保存绝对路径
extern int g_preflag;             // 是否预处理标记

// ApproximateReversePush
extern DWORD a_mmfsizehigh;                // 高位文件大小（x4G）
extern DWORD a_mmfsizelow;                 // 低位文件大小
extern unsigned int a_mmf_buffer_size;     // 单个映射文件的大小
extern int a_THREADNUM;                    // 近似方案多线程任务数

// PartialReversePush
extern DWORD p_mmfsizehigh;                // 高位文件大小（x4G）
extern DWORD p_mmfsizelow;                 // 低位文件大小
extern unsigned int p_mmf_buffer_size;     // 单个映射文件的大小
extern int p_THREADNUM;                    // 部分计算方案多线程任务数

#endif