#ifndef EDGEHEAP_H
#define EDGEHEAP_H

#include <queue>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <limits>
#include "Cartesian3.h"
#include "Matrix4.h"

struct EdgeDelta
{
    int32_t edgeId;
    int32_t from;
    int32_t to;
    double delta;
    Cartesian3 newPos;
    Matrix4 edgeQ;
    EdgeDelta(int32_t from = -999999,int32_t to= -999999):from(from),to(to){
        delta = std::numeric_limits<double>::max();
    }

    auto operator() (const EdgeDelta & X, const EdgeDelta & Y) const -> bool{
         return X.delta > Y.delta;
    }
};

class EdgeHeap
{
public:
    EdgeHeap();
    auto addEdge(EdgeDelta & edge) -> void;
    auto delEdge(const EdgeDelta & e)  -> void;
    auto getMinDelta(const std::vector<bool> & vertexDeleted) -> EdgeDelta;
    auto clear() -> void;
private:
    std::priority_queue<EdgeDelta,std::vector<EdgeDelta>,EdgeDelta> heap;
    std::map<std::pair<int, int>, int> mapEdgeToID;
    int32_t numerOfEdge = 0;
    std::unordered_map<int32_t,bool> deletedMap;


};

#endif // EDGEHEAP_H
