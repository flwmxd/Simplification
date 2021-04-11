#include "EdgeHeap.h"

EdgeHeap::EdgeHeap()
{

}

auto EdgeHeap::addEdge(EdgeDelta & e) -> void
{
    numerOfEdge++;
    e.edgeId = numerOfEdge;
    int u = std::min(e.from,e.to);
    int v = std::max(e.from,e.to);
    if(mapEdgeToID.find(std::make_pair(u,v)) != mapEdgeToID.end()){
        auto id = mapEdgeToID[std::make_pair(u,v)];
        e.edgeId = id;
        deletedMap[id] = false;
    }else{
        mapEdgeToID[std::make_pair(u,v)] = numerOfEdge;
    }
    heap.push(e);
}

auto EdgeHeap::delEdge(const EdgeDelta & e)  -> void
{
    int32_t u = std::min(e.from,e.to);
    int32_t v = std::max(e.from,e.to);
    int32_t ID = mapEdgeToID[std::make_pair(u,v)];
    deletedMap[ID] = true;
}

auto EdgeHeap::clear() -> void
{
    deletedMap.clear();

    heap = {};
    numerOfEdge = 0;
}
auto EdgeHeap::getMinDelta(const std::vector<bool> & vertexDeleted) -> EdgeDelta
{
    if(heap.size() <= 0){
        return EdgeDelta(0,0);
    }
    while(deletedMap[heap.top().edgeId] || vertexDeleted[heap.top().to] || vertexDeleted[heap.top().from]){
        heap.pop();
    }
    auto e = heap.top();
    heap.pop();
    deletedMap[e.edgeId] = true;
    return e;
}
