#ifndef TARJAN_H
#define TARJAN_H

#include <vector>
#include <cstdint>
#include <Cartesian3.h>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <memory>


class GeometricSurfaceDirectedEdge;

struct Mesh
{
    std::vector<Cartesian3> position;
    std::vector<signed long> faceVertices;
    std::vector<signed long> firstDirectedEdge;
    std::vector<signed long> otherHalf;
    std::vector<Cartesian3> normal;
};

class Tarjan
{
public:
    Tarjan();
    //Tarjan's algorithm with none-recursive
    auto start(GeometricSurfaceDirectedEdge * edge) -> void;
private:
    auto generateMeshes(GeometricSurfaceDirectedEdge * edge,const std::vector<std::vector<int32_t>> & faces) -> void;
    auto generateDirectedEdge(Mesh & mesh) -> void;
    auto split(std::string input, const std::string& delimiter) -> std::vector<std::string>;
    auto generateNeighbours(GeometricSurfaceDirectedEdge * edge) -> void;
    auto isConnected(int32_t v1, int32_t v2) -> bool;
    std::unordered_map<int32_t,std::unordered_set<int32_t>> neighboursMap;
};

#endif // TARJAN_H
