#include "../src/graphNO.hpp"
#include "../src/weightedMaximumStableSolver.hpp"
#include "../src/solver_to_csv.hpp"
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <dirent.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <string>

std::chrono::duration<double> print_solve_stats( WeightedMaximumStableSolver *solver, string display, char tMode, char eMode, char pMode)
{
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 
//   std::cout  << display << std::endl;
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 

  auto start = std::chrono::steady_clock::now();
  solver->solverGreedy(tMode, eMode, pMode);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
//   std::cout << "stable max: " << solver->getCurrentNB() << endl;
//   std::cout << "poids : " << solver->getCurrentCost() << endl;
//   std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
  solver->updateBestSolution();
//   std::cout << std::endl;

  return elapsed_seconds;
}

void print_and_save_solve_stats(string inst_name, WeightedMaximumStableSolver *solver, string display, char tMode, char eMode, char pMode)
{
   std::chrono::duration<double> elapsed_seconds = print_solve_stats(solver, display, tMode, eMode, pMode);
   tocsv_solve_stats(inst_name, display, solver, elapsed_seconds);
}

std::chrono::duration<double> print_solveweight_stats( WeightedMaximumStableSolver *solver, string display, char tMode)
{
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 
//   std::cout  << display << std::endl;
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 

   auto start = std::chrono::steady_clock::now();
   solver->solverGreedyWeight(tMode);
   auto end = std::chrono::steady_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
//   std::cout << "stable max: " << solver->getCurrentNB() << endl;
//   std::cout << "poids : " << solver->getCurrentCost() << endl;
//   std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
   solver->updateBestSolution();
//   std::cout << std::endl;

  return elapsed_seconds;
}

void print_and_save_solveweight_stats(string inst_name, WeightedMaximumStableSolver *solver, string display, char tMode)
{
   std::chrono::duration<double> elapsed_seconds = print_solveweight_stats(solver, display, tMode);
   tocsv_solve_stats(inst_name, display, solver, elapsed_seconds);
}

std::chrono::duration<double> print_solvedegree_stats( WeightedMaximumStableSolver *solver, string display, char tMode)
{
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 
//   std::cout  << display << std::endl;
//   std::cout << "___________________________________________________________________________________________________________ \n" << std::endl; 

   auto start = std::chrono::steady_clock::now();
   solver->solverWeightDegree(tMode);
   auto end = std::chrono::steady_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
//   std::cout << "stable max: " << solver->getCurrentNB() << endl;
//   std::cout << "poids : " << solver->getCurrentCost() << endl;
//   std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
   solver->updateBestSolution();
//   std::cout << std::endl;

  return elapsed_seconds;
}

void print_and_save_solvedegree_stats(string inst_name, WeightedMaximumStableSolver *solver, string display, char tMode)
{
   std::chrono::duration<double> elapsed_seconds = print_solvedegree_stats(solver, display, tMode);
   tocsv_solve_stats(inst_name, display, solver, elapsed_seconds);
}

void solver_weight(WeightedMaximumStableSolver &solver, string inst_name)
{
   {
      auto start = std::chrono::steady_clock::now();
      solver.solverDegNeighbors();
      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = end-start;
    
      // std::cout  << "SolverWeightNeighbors" << std::endl; 
      // std::cout << "stable max: " << solver.getCurrentNB() << endl;
      // std::cout << "poids : " << solver.getCurrentCost() << endl;
      // std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
      // std::cout << std::endl;
    
      solver.updateBestSolution();

      tocsv_solve_stats(inst_name, "SolverWeightNeighbors", &solver, elapsed_seconds);
   }
  
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyWeightsDecDegC", '0');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyWeightsDecDegDec", '1');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyWeightsCDegC", '2');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyWeightsCDegDec", '3');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyDegCWDec", '4');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyDegCWC", '5');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyDegDecWDec", '6');
   print_and_save_solveweight_stats(inst_name, &solver, "SolverGreedyDegDecWC", '7');
}

void solver_degree(WeightedMaximumStableSolver &solver, string inst_name)
{
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Degree/Weights Max Id Max", '0');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Degree/Weights Max Id Min", '1');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Degree/Weights Min Id Max", '2');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Degree/Weights Min Id Min", '3');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Weights/Degree Max Id Max", '4');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Weights/Degree Max Id Min", '5');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Weights/Degree Min Id Max", '6');
   print_and_save_solvedegree_stats(inst_name, &solver, "SolverRatio Weights/Degree Min Id Min", '7');
}

void solver_random(WeightedMaximumStableSolver &solver, string inst_name)
{
   // Simple sort -> randoms -> best
  {
    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsCr " 
    "EgaRand "
    "Best", 
    POIDS_CROISSANTS, EGA_RANDOM, PICK_BEST);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsDecr " 
    "EgaRand "
    "Best",  
    POIDS_DECROISSANTS, EGA_RANDOM, PICK_BEST);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegCr " 
    "EgaRand "
    "Best", 
    DEG_CROISSANTS, EGA_RANDOM, PICK_BEST);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegDecr " 
    "EgaRand "
    "Best", 
    DEG_DECROISSANTS, EGA_RANDOM, PICK_BEST);
  }
  
  // Simple sort -> randoms -> best 5 -> random
  {
    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsCr " 
    "EgaRand "
    "Best5 rand",   
    POIDS_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsDecr " 
    "EgaRand "
    "Best5 rand",  
    POIDS_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegCr " 
    "EgaRand "
    "Best5 rand",   
    DEG_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegDecr " 
    "EgaRand "
    "Best5 rand",  
    DEG_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM);
  }
  
  // Simple sort -> randoms -> best 5 -> random (weight)
  {
    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsCr " 
    "EgaRand "
    "Best5 rand (poids)",   
    POIDS_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_WEIGHT);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsDecr " 
    "EgaRand "
    "Best5 rand (poids)",  
    POIDS_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_WEIGHT);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegCr " 
    "EgaRand "
    "Best5 rand (poids)",  
    DEG_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_WEIGHT);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegDecr " 
    "EgaRand "
    "Best5 rand (poids)", 
    DEG_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_WEIGHT);
  }
  
  // Simple sort -> randoms -> best 5 -> random (degree)
  {
    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsCr " 
    "EgaRand "
    "Best5 rand (degree)",   
    POIDS_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_DEGREE);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "PoidsDecr " 
    "EgaRand "
    "Best5 rand (degree)",  
    POIDS_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_DEGREE);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegCr " 
    "EgaRand "
    "Best5 rand (degree)",  
    DEG_CROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_DEGREE);

    print_and_save_solve_stats(inst_name, &solver, 
    "SolverRandom " 
    "DegDecr " 
    "EgaRand "
    "Best5 rand (degree)", 
    DEG_DECROISSANTS, EGA_RANDOM, PICK_BEST5_RANDOM_DEGREE);
  }
}

int main(int argc, char** argv) 
{
   vector<char*> filesNames;
   struct dirent *dir;
   DIR *d = opendir("../weights"); 
   if (d) {
      while ((dir = readdir(d)) != NULL) {
         //if ((dir->d_name).compare(".") && (dir->d_name).compare("..")){
         //if (!strcmp(dir->d_name, ".") && !strcmp(dir->d_name, "..")){
            filesNames.push_back(strdup(dir->d_name));
         //}
      }
      closedir(d);
   }
   std::cout << filesNames.size() << std::endl;
   
   WeightedMaximumStableSolver solver;
   char* file;
   if(argc == 1) file = filesNames[9];
   //else file = filesNames[stoi(string(argv[1]))];
   else file = argv[1];

   std::cout << "file name : " << file << std::endl;
   char buffer[1000]={};
   solver.importWeights(file);
   char* path = strcat(strcat(buffer,strdup("../instances/")),file);
   std::cout << "instances path : " << path << std::endl;
   solver.importGraphDIMACS(path);

   //On coupe juste le nom de l'instance pour le fichier csv
   string inst_name = "";
   string tmp = "";
   for(char c = *file; c; c=*(++file)) 
   {
      if(c != '.') tmp += c;
      else { inst_name+= tmp; tmp = ""; }
   }

   //On clear le fichier csv de l'instance
   clear_instance_file(inst_name);

   //On sauvegarde les stats des solvers dans le nouveau fichier csv
   solver_random(solver, inst_name);
   solver_degree(solver, inst_name);
   solver_weight(solver, inst_name);

   tocsv_custom_stats(inst_name, "Best Solution", solver.getBestSolutionVertexCount(), solver.getBestSolutionWeight(), 0);

   return 0;
}