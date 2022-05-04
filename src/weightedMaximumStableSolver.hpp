#ifndef WEIGHTEDMAXIMUMSTABLESOLVER_HPP
#define WEIGHTEDMAXIMUMSTABLESOLVER_HPP

#include "graphNO.hpp"
using namespace std;

#define POIDS_CROISSANTS '0'
#define POIDS_DECROISSANTS '1'
#define DEG_CROISSANTS '2'
#define DEG_DECROISSANTS '3'

#define EGA_RANDOM 'A'
#define EGA_BEST_WEIGHT 'B'
#define EGA_WORST_WEIGHT 'C'
#define EGA_BEST_DEGREE 'D'
#define EGA_WORST_DEGREE 'E'

#define PICK_BEST 'A'
#define PICK_BEST5_RANDOM 'B'
#define PICK_BEST5_RANDOM_WEIGHT 'C'
#define PICK_BEST5_RANDOM_DEGREE 'D'


class WeightedMaximumStableSolver {

private:
   GraphNO graph;
   Vertex nbVertex;

   float currentCost;
   float currentNB;
   float bestCost;

   std::vector<float> weights;
   std::vector<bool> currentSolution;
   std::vector<bool> bestSolution;

   void initSolution();

public:

   void importGraphDIMACS( char *file);
   void importWeights( char *file);
   void displayBestSolution();
   bool checkSolution();
   void updateBestSolution();

   int getBestSolutionVertexCount(){ int count = 0; for(bool b : bestSolution) if(b) count++; return count;}
   float getBestSolutionWeight() { return bestCost; }

   void solver();
   void solverGreedyWeight(char triMode);
   void solverDegNeighbors();

   void decroissant(vector<std::pair<Vertex,Vertex>> v);
   std::pair<Vertex,Vertex> getSaturation(Vertex v);
   vector<std::pair<Vertex,Vertex>> getSaturations();
   std::pair<Vertex,Vertex> getDegree(Vertex v);
   vector<std::pair<Vertex,Vertex>> getDegrees();
   vector<std::pair<Vertex,Vertex>> getNotColored(vector<std::pair<Vertex,Vertex>> v);
   vector<Vertex> getColorNeighbors(Vertex v);
   void croissant(vector<Vertex> v);
   vector<std::pair<Vertex,Vertex>> equalCandidat(vector<std::pair<Vertex,Vertex>> v, Vertex level);
   int getNBColors();
   void croissantFirst(vector<std::pair<Vertex,Vertex>> v);
   float getCurrentCost() { return currentCost; };
   float getCurrentNB() { return currentNB; };

   void solverWeightDegree(char mode);
   void egalite(std::pair<Vertex, Vertex> &i, std::pair<Vertex, Vertex> &j, char egaMode);
   void pick(vector<std::pair<Vertex,float>> sortedVertex, char pickMode);
   void solverGreedy(char triMode);
   void solverGreedy(char triMode, char egaMode, char pickMode);
   void addVertexToSolution(Vertex v);



   void decroissant(vector<std::pair<Vertex,float>>* vector);

};
#endif
