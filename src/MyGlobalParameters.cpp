#include "MyGlobalParameters.h"

using namespace std;

// global parameters
int g_clustertype = STRUCTURE;
float g_alpha = 0.2f;
float g_gamma = 0.0f;

int g_datasetid = 0;
int g_schemeid = 0; 

int ATTRIBUTE_NUM = 0;
float g_epsilon = 0.0f;
float g_delta = 0.0f;
int g_minPts = 0;
int g_cluster_startid = 0;
int g_cluster_endid = 0;
int g_vertexnum = 0;
long long g_pushbackcount = 0;
vector<Vertex*> g_vertices;
map<int, set<int>> g_dbscanneighborsets;

float * g_bufferpool[g_buff_size];
MMF * g_mmfpool[g_buff_size];
float * g_mmfm_address_pool[g_buff_size];
set<int> g_buffer_queue;
float * g_pr_pool[g_buff_size];

HANDLE  g_hSemaphoreRunnerNum;
HANDLE  g_hThreadEvent;
CRITICAL_SECTION g_cs;

// ±ﬂ»®÷ÿ
float g_structure_edgeweight = 1.0f;
float g_attribute_total_edgeweight = 0.0f;
float g_total_edgeweight = 0.0f;
float g_base_weight[ATTRIBUTE_BUFF_SIZE] = { 1.0f, 0.0f, 0.0f, 0.0f };
float g_former_edgeweight[ATTRIBUTE_BUFF_SIZE][ATTRIBUTE_BUFF_SIZE] =
{
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 }
};
float g_edgeweight[ATTRIBUTE_BUFF_SIZE][ATTRIBUTE_BUFF_SIZE] =
{
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0, 0.0 }
};

string g_inputpath = "F:\\WorkSpace\\GraphClustering\\GC_ApproximateReversePush\\data\\";
string g_resultpath = "F:\\WorkSpace\\GraphClustering\\GC_Result\\";

string g_mmfPath = "F:\\WorkSpace\\GraphClustering\\GC_ApproximateReversePush\\data\\PPR\\";

int g_preflag = 0;

// ApproximateReversePush
unsigned int a_mmf_buffer_size = 0;
DWORD a_mmfsizehigh = 0;
DWORD a_mmfsizelow = 0;
int a_THREADNUM = 0;

// PartialReversePush
unsigned int p_mmf_buffer_size = 0;
DWORD p_mmfsizehigh = 0;
DWORD p_mmfsizelow = 0;
int p_THREADNUM = 0;