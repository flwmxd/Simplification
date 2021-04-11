#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

#include <unordered_set>
#include "Matrix4.h"
#include "EdgeHeap.h"


/**
 * @brief The Simplifier class
 * Surface Simplification Using Quadric Error Metrics
 */
class Simplifier
{
public:
    Simplifier();
    auto generateNeighbours(int32_t vertexId) -> void;
    auto simplification(float ratio) -> void;
    auto calculateError(EdgeDelta & delta) -> void;
    auto calVertexQ(int32_t vertexId) -> Matrix4;
    auto calVertexNewPos(const EdgeDelta & edge,const Matrix4 & matrix) -> Cartesian3;
    auto getCommonNum(int32_t v1, int32_t v2) -> int32_t;
    auto isConnected(int32_t v1, int32_t v2) -> bool;
    auto removeConnect(int32_t v1, int32_t v2) -> void;
    auto removeConnect(int32_t v) -> void;
    auto calculateCost() -> void;
    auto createNewMesh( std::vector< signed long > & outMesh) -> void;

private:
    auto start() -> void;
    friend class GeometricSurfaceDirectedEdge;
    std::unordered_map<int32_t,std::unordered_set<int32_t>> neighboursMap;
    std::vector<Matrix4> errorQ;
    float ratio = 1.f;
    EdgeHeap heap;
    std::vector<bool> deletedMarks;
    std::vector<Cartesian3> position;
    std::vector< signed long > faceVertices;
};

#endif // SIMPLIFIER_H
