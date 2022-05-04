#ifndef GRAPHNO_HPP
#define GRAPHNO_HPP

#include <vector>
using namespace std;

typedef unsigned int  Vertex;

class GraphNO {

private:
   vector<bool> matrix;
   Vertex nbVertex;
   void get_params(char * Preamble);
   void removeEdge(Vertex v1, Vertex v2);

public:
   void initEmptyGraph(Vertex nbVertices);
   bool isEdge(Vertex v1, Vertex v2);
   bool hasNoNeighbor(Vertex v1, const std::vector<Vertex> & v2);
   void markingNeighbors(Vertex v1, std::vector<bool> & v);
   void addEdge(Vertex v1, Vertex v2);
   Vertex getNbVertices(){return nbVertex;};
   void deleteGraph(){matrix.clear();nbVertex=0;};
   void importGraphDIMACS( char *file);
   void display();
   Vertex getDegree(Vertex vertex);
   Vertex getDegreeWithUnmarked(Vertex vertex,  const vector<bool> & marked);
   void toComplementGraph();

   vector<Vertex> getNeighbors(Vertex vertex);
   void intersectNeighbors(Vertex vertex ,  vector<Vertex> & remaining);
};
#endif
