#ifndef Clustering_Vertex_h
#define Clustering_Vertex_h

#include <vector>
#include <string>

using std::string;
using std::vector;

class Vertex
{
public:
	int vertexid;                             // 节点Id
	int vertextype;                           // 节点类型
	string vertexname;                        // 节点名称（统一32位编号）
	vector<int> edgetype_outdegree;           // 每种类型的邻边数目，下标是类型. DBLP中顺序是 p a c k; Flickr中顺序是image, user, tag
	vector<int> neighborvertex;               // 邻居节点编号（按照outDegrees类别的顺序排列）

	Vertex(int _vertexid, int _vertextype, string vertexname) :vertexid(_vertexid), vertextype(_vertextype), vertexname(vertexname){}
};

#endif
