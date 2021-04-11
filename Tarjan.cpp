#include "Tarjan.h"
#include "GeometricSurfaceDirectedEdge.h"
#include <stack>
#include <algorithm>

Tarjan::Tarjan()
{

}

//Tarjan's algorithm with none-recursive
auto Tarjan::start(GeometricSurfaceDirectedEdge * edge) -> void
{
    //the result. vertices collection for every mesh.
     std::vector<std::vector<int32_t>> components ;


    //generate the neighbours for mesh because I need the relation among all edges;
    generateNeighbours(edge);


    //the stack for visited vertex
    std::stack<int32_t> sccStack;
    //indicite for visited vertex.
    std::unordered_map<int32_t,bool> sccFound;
    //indicite for visited time.
    std::unordered_map<int32_t,int32_t> preorder;

    //marks for the min value for vertex in a graph
    std::vector<int32_t> lowValue;
    //init it as -1.
    lowValue.resize(edge->position.size(),-1);

    int32_t time = 0;

//iterate every vertex;
    for(size_t i = 0;i<edge->position.size();i++){
//to avoid running out of the stack memeory,
//using stack to emulate the recursive process..
        std::stack<int32_t> stack;
        if(!sccFound.count(i))
        {
            stack.push(i);
            while(!stack.empty()){
                auto v = stack.top();
                if(!preorder.count(v)){//calculate the time for vertex;
                    preorder[v] = ++time;
                }
                bool done = true;
                //iterate neighbours..
                for(auto w : neighboursMap[v]){
                    if(!preorder.count(w)){
                        stack.push(w);
                        done = false;
                        break;
// after break the iteration, the program will execute the next vertex.
                    }
                }

                //visited all vertices
                if(done)
                {

                    lowValue[v] = preorder[v];
//set all lower vale for every vertex in current mesh or graph
                    for(auto w : neighboursMap[v]){
                        if(!sccFound.count(w)){
                            if (preorder[w] > preorder[v]) {
                                lowValue[v] =std::min(lowValue[v], lowValue[w]);
                            } else {
                                lowValue[v] =std::min(lowValue[v], preorder[w]);
                            }
                        }
                    }

                    stack.pop();

                    if(lowValue[v] == preorder[v])
                    {
                        sccFound[v] = true;
                        std::vector<int32_t> mesh;
                        mesh.emplace_back(v);
//if all vertex lower value are the same. it belong to a common mesh or graph.
//craete it..
                        while (!sccStack.empty() && preorder[sccStack.top()] > preorder[v])
                        {
                             auto k = sccStack.top();
                            sccStack.pop();
                            sccFound[k] = true;
                            mesh.emplace_back(k);
                        }
                        components.emplace_back(mesh);
                    }
                    else
                    {
                         sccStack.push(v);
                    }
                }
            }
        }
    }


    std::cout<<"there are "<<components.size()<<" components"<<std::endl;

    generateMeshes(edge,components);

}


auto Tarjan::generateMeshes(GeometricSurfaceDirectedEdge * edge,const std::vector<std::vector<int32_t>> & faces) -> void
{

    std::unordered_map<int32_t,int32_t> vertexIdMap;
    std::vector<Mesh> meshes;

    for(auto & vertices : faces)
    {
        //allocate a mesh
        meshes.emplace_back(Mesh{});
        auto & m = meshes.back();
        for(auto id : vertices)
        {
            //re-mapping the vertex id
            if(vertexIdMap.count(id) == 0)
            {
                auto & pos = edge->position[id];
                vertexIdMap[id] = m.position.size();
                m.position.emplace_back(pos);
                m.normal.emplace_back(edge->normal[id]);
            }
        }

        //create new faces with new position.
        for(int32_t i = 0;i<vertices.size();i++){

            auto original = vertices[i];
            auto & nei = neighboursMap[original];
            for (auto id : nei)
            {
                if( original >= id){
                    continue;
                }

                for (auto id2 : nei)
                {
                    //check the three dots are connected with each other.
                    //so add into faceVertices.
                    if(id2 > id && isConnected(id2,id)){
                        m.faceVertices.emplace_back(vertexIdMap[id2]);
                        m.faceVertices.emplace_back(vertexIdMap[id]);
                        m.faceVertices.emplace_back(vertexIdMap[original]);
                    }
                }
            }
        }

        //finish a mesh
        vertexIdMap.clear();
        //crate direct edge..
        generateDirectedEdge(m);
    }


    //create new filename ..
    std::string name = edge->originalName;
    auto pos = name.find_last_of('/');
    if(pos != std::string::npos)
    {
        name = name.substr(pos + 1,name.size() - 1);
    }
    auto ret = split(name,".");



    //write out
    for(int32_t i = 0;i<meshes.size();i++)
    {
        char fileName[256] = {0};

        sprintf(fileName,"connected-component-%s-%d.%s",ret[0].c_str(),i,ret[1].c_str());
        edge->WriteFileDirEdge(fileName,
                               meshes[i].position,
                               meshes[i].firstDirectedEdge,
                               meshes[i].faceVertices,
                               meshes[i].otherHalf,
                               meshes[i].normal
                               );
        std::cout<<"write out "<<fileName<<" success"<<std::endl;
    }

}


//split string by delimiter.... just a tools function..
auto Tarjan::split(std::string input, const std::string& delimiter) -> std::vector<std::string>
{
    std::vector<std::string> ret;
    size_t pos = 0;
    std::string token;
    while ((pos = input.find(delimiter)) != std::string::npos)
    {
        token = input.substr(0, pos);
        ret.push_back(token);
        input.erase(0, pos + delimiter.length());
    }
    ret.push_back(input);
    return ret;
}


//generate DirectedEdge.. comes from prev assignments...
auto Tarjan::generateDirectedEdge(Mesh & mesh) -> void
{

    auto getPrevEdge = [](int32_t index) -> int32_t
    {
        return (index % 3) != 0 ? index - 1 : index + 2;
    };

    constexpr auto INIT_VALUE = std::numeric_limits<int32_t>::max();

    mesh.firstDirectedEdge.clear();
    mesh.otherHalf.clear();
    //init wit a unreachable value to mark it is not used
    mesh.otherHalf.resize(mesh.faceVertices.size(),INIT_VALUE);
    mesh.firstDirectedEdge.resize(mesh.faceVertices.size(),INIT_VALUE);

    //generate the half data structure;
    for(auto edgeId = 0;edgeId<mesh.faceVertices.size();edgeId++)
    {
       auto faceId = edgeId / 3;
       //to -> from

       // Hint 1 :
       // if a face with vertices was interpreted as from
       // we can do like that
       // auto endVert = faceVertices[getNextEdge(edgeId)];
       auto startVert = mesh.faceVertices[edgeId];
       auto endVert = mesh.faceVertices[getPrevEdge(edgeId)];

       if(startVert != endVert)
       {
           if(mesh.firstDirectedEdge[startVert] == INIT_VALUE){
               mesh.firstDirectedEdge[startVert] = edgeId;
           }

           if(mesh.otherHalf[edgeId] == INIT_VALUE){
               //foreach the remain faces;
               for(auto nextEdgeId = 3 * (faceId + 1); nextEdgeId < mesh.faceVertices.size();nextEdgeId++)
               {
                   auto nextStart = mesh.faceVertices[nextEdgeId];
                   //same with Hint1
                   auto nextEnd= mesh.faceVertices[getPrevEdge(nextEdgeId)];

                   if((startVert == nextEnd && nextStart == endVert)
                           ||
                       (startVert == nextStart && nextEnd == endVert)
                           )
                   {
                       mesh.otherHalf[edgeId] = nextEdgeId;
                       mesh.otherHalf[nextEdgeId] = edgeId;
                   }
               }
           }
       }
    }
}


//same with simplifier...
auto Tarjan::isConnected(int32_t v1, int32_t v2) -> bool
{
    return neighboursMap[v1].count(v2) > 0 || neighboursMap[v2].count(v1)  > 0 ;
}



auto Tarjan::generateNeighbours(GeometricSurfaceDirectedEdge * edge) -> void
{
    for(int32_t vertexId = 0; vertexId <edge->position.size(); vertexId++){

        std::unordered_set<int32_t> neighbours;

        for(int32_t i = 0; i < edge->faceVertices.size(); i+=3)
        {
            if(edge->faceVertices[i] == vertexId)
            {
                neighbours.emplace(edge->faceVertices[i+1]);
                neighbours.emplace(edge->faceVertices[i+2]);
            }
            else if (edge->faceVertices[i+1] == vertexId)
            {
                neighbours.emplace(edge->faceVertices[i]);
                neighbours.emplace(edge->faceVertices[i+2]);

            }
            else if (edge->faceVertices[i+2] == vertexId)
            {
                neighbours.emplace(edge->faceVertices[i+1]);
                neighbours.emplace(edge->faceVertices[i]);
            }
        }
        neighboursMap[vertexId]= neighbours;
    }
}
