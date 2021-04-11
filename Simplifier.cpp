#include "Simplifier.h"
#include <cassert>
#include <cmath>

Simplifier::Simplifier()
{

}



//calculate the  neighbours for every vertex.
auto Simplifier::generateNeighbours(int32_t vertexId) -> void
{
    std::unordered_set<int32_t> neighbours;

    for(auto i = 0; i < faceVertices.size(); i+=3)
    {
        if(faceVertices[i] == vertexId)
        {
            neighbours.emplace(faceVertices[i+1]);
            neighbours.emplace(faceVertices[i+2]);
        }
        else if (faceVertices[i+1] == vertexId)
        {
            neighbours.emplace(faceVertices[i]);
            neighbours.emplace(faceVertices[i+2]);

        }
        else if (faceVertices[i+2] == vertexId)
        {
            neighbours.emplace(faceVertices[i+1]);
            neighbours.emplace(faceVertices[i]);
        }
    }

    neighboursMap[vertexId]= neighbours;
}



auto Simplifier::simplification(float ratio) -> void
{
    //clear pre-result
    deletedMarks.clear();
    deletedMarks.resize(position.size(),false);
    neighboursMap.clear();

    for(auto id = 0 ;id < position.size();id++){
       generateNeighbours(id);
    }
    //calculate cost
    calculateCost();

    this->ratio = ratio;

    //start simplification
    start();
}


/**
 * 1. first step Compute the Q matrices for all the initial vertices
 * and
 * 2. select all valid pair..
  */
auto Simplifier::calculateCost() -> void
{
    heap.clear();
    errorQ.clear();
    errorQ.resize(position.size());

    for(auto id = 0; id <position.size(); id++){
        errorQ[id] = calVertexQ(id);
    }
    std::cout<<"finished errorQ "<<std::endl;


    for(size_t id = 0 ;id < position.size();id++)
    {
        if(!deletedMarks[id])
        {
            for(auto nId : neighboursMap[id])
            {
                if(id < nId){
                    break;
                }
                EdgeDelta edge;
                edge.from = nId;
                edge.to = id;
                calculateError(edge);

                if(edge.delta != 0){
                    //push it to the heap...
                    heap.addEdge(edge);
                }
            }
        }
     }
}

auto Simplifier::start() -> void
{

    auto delFace = (1.f - ratio)  * (faceVertices.size() / 3.f);

    std::cout<<"start ..total face is "<<faceVertices.size() / 3.f<<"  delFace : "<<delFace<<std::endl;

    //made a backup ,because we can't delete any elements when iterate..
    std::unordered_map<int32_t,std::unordered_set<int32_t>> backupNeighbours(neighboursMap.begin(),neighboursMap.end());

    if(faceVertices.size() / 3.f - delFace < 2){
        delFace = faceVertices.size() / 3.f  - 2;
    }

    /**
      * final step
      * Iteratively remove the pair (v1, v2 ) of least cost from the heap,
      */
    for(int32_t i = 0;i<delFace;i+= 2){
        auto e = heap.getMinDelta(deletedMarks);

        if(e.from == 0 && e.to == 0){
            break;
        }

        position[e.from] = e.newPos;
        errorQ[e.from] = e.edgeQ;

        //delete e.to
        deletedMarks[e.to] = true;

        std::unordered_set<int32_t> connected;
        connected.clear();

        //here, two vertices (e.from and e.to) will be combined,
        //the neighbours of original vertix (e.from) will keep relation with the original vertex
        //the neighbours of second vertex(e.to) will keep relation with e.from

        removeConnect(e.from,e.to);

        auto & neighbours = backupNeighbours[e.to];
        for (auto id : neighbours){
            if(id != e.from){
                heap.delEdge({id,e.to});
                removeConnect(id,e.to);
                connected.emplace(id);
            }
        }


        removeConnect(e.to);


        for(auto to : connected){
            neighboursMap[e.from].emplace(to);
            neighboursMap[to].emplace(e.from);
        }

        //update the costs of all valid pairs involving v1
        for(auto to : neighboursMap[e.from]){
            if(e.from != to){
                EdgeDelta edge(e.from,to);
                heap.delEdge(edge);

                calculateError(edge);
                if(edge.delta != 0)
                    heap.addEdge(edge);
            }
        }

        //backup new neighbours;
        backupNeighbours = neighboursMap;
    }
}

/**
 * @brief Simplifier::create New Mesh
 * @param outMesh
 */
auto Simplifier::createNewMesh(std::vector< signed long > & outMesh) -> void
{
    outMesh.clear();

    for(int32_t i = 0;i<position.size();i++)
    {
        //if current position is deleted,continue the loop
        if(deletedMarks[i]){
            continue;
        }

        auto & nei = neighboursMap[i];
        for (auto id : nei)
        {
            if( i >= id){
                continue;
            }

            for (auto id2 : nei)
            {
                //check the three dots are connected with each other.
                //so add into faceVertices.
                if(id2 > id && isConnected(id2,id)){
                    outMesh.emplace_back(i);
                    outMesh.emplace_back(id);
                    outMesh.emplace_back(id2);
                }
            }
        }
    }

    std::cout<<"the number of current faces : " <<outMesh.size() / 3.f<<std::endl;
}


//Q = SUM of Kp(P*Pt) , P = [a,b,c,d]T
auto Simplifier::calVertexQ(int32_t vertexId) -> Matrix4
{

    Matrix4 ans;
    ans.SetZero();

    auto & connected = neighboursMap[vertexId];

    const auto & p = position[vertexId];

    //iterate every related faces
    for(auto t1 : connected){
        for(auto t2 : connected){
            if(t1 < t2 && isConnected(t1,t2)){

                const auto & v1 = position[t1];
                const auto & v2 = position[t2];

                auto normal =  (v1 - p).cross(v2 - p).normalise();

                //P = [a,b,c,d]T
                float vec4[4] = {normal.x,normal.y,normal.z,-(normal.dot(p))};
                /* Ax + Bx + Cx + D = 0 */
                /* A,B,C,D should be normalized */


                /*
                 * p = [a b c d]
                 * p * pT = Kp
                 *
                 *  [a2 ab ac ad]
                 *  [ab b2 bc bd]
                 *  [ac bc c2 cd]
                 *  [ad bd cd d2]
                 */

                for(int i = 0;i < 4;i++){
                    for(int j = 0;j < 4;j++){
                        ans[i][j] += vec4[i] * vec4[j];
                    }
                }
            }
        }
    }
   return ans;
}

auto Simplifier::calVertexNewPos(const EdgeDelta & edge,const Matrix4 & matrix) -> Cartesian3
{

    auto m = matrix;
    m[3][0] = 0;
    m[3][1] = 0;
    m[3][2] = 0;
    m[3][3] = 1;

    Homogeneous4 y = {0,0,0,1};

    assert(!deletedMarks[edge.from] && !deletedMarks[edge.to]);
    //the matrix(m) can be inverted...
    if(std::fabs(m.det()) < 0.000001){
        //some optimization ..
        auto mid = (position[edge.from] + position[edge.to]) / 2;
        auto cost1 = (matrix * position[edge.to]).dot(position[edge.to]);
        auto cost2 = (matrix * position[edge.from]).dot(position[edge.from]);
        auto costMid = (matrix * mid).dot(mid);

        if(cost1 < cost2){
            if(cost1 < costMid){
                return position[edge.to];
            }
            return mid;
        }
        if(cost2 < costMid){
            return position[edge.from];
        }
        /*
         * If the matrix is not invertible
         * we use the mid-point
         */
        return  mid;//center point ..
    }

    //V = A^-1 * Y
    return (m.inverse() * y).Point();
}

/**
 * 3. step 3
 * Compute the optimal contraction target v¯ for each valid pair
 * (v1, v2 ). The error v¯T(Q1 +Q2 )v¯ of this target vertex becomes
 * the cost of contracting that pair.
 */
auto Simplifier::calculateError(EdgeDelta & e) -> void
{
    //compute the error and calculate the sum;
    Matrix4 mat = errorQ[e.to] + errorQ[e.from];
    e.newPos = calVertexNewPos(e,mat);
    e.edgeQ = mat;

    /*if (getCommonNum(e.from,e.to) != 2) {
        e.delta = 0;
        return;
    }*/
    Homogeneous4 x = e.newPos;
    e.delta = (mat * x).dot(x);
}

//remove the connectivity between v1 and v2..
auto Simplifier::removeConnect(int32_t v1, int32_t v2) -> void
{
    neighboursMap[v1].erase(v2);
    neighboursMap[v2].erase(v1);
}

//remove vertex in neighboursMaps..
auto Simplifier::removeConnect(int32_t v) -> void
{
    neighboursMap.erase(v);
}

//judge whether v1 and v2 can be constructed as an edge..
auto Simplifier::isConnected(int32_t v1, int32_t v2) -> bool
{
   return neighboursMap[v1].count(v2) > 0 || neighboursMap[v2].count(v1) > 0;
}

auto Simplifier::getCommonNum(int32_t v1, int32_t v2) -> int32_t
{
    auto & conn1 = neighboursMap[v1];
    auto & conn2 = neighboursMap[v2];
    auto count = 0;
    for(auto id : conn1){
        for(auto id2 : conn2){
            if(id == id2){
                count++;
            }
        }
    }
    return count;
}




















