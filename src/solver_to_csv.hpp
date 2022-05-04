#ifndef SOLVER_TO_CSV
#define SOLVER_TO_CSV

using namespace std;

#include "weightedMaximumStableSolver.hpp"
#include <fstream>
#include <chrono>
#include <iostream>

void tocsv_solve_stats( string instance_name, string display, WeightedMaximumStableSolver *solver, std::chrono::duration<double> time);
void tocsv_custom_stats(string instance_name, string display, int nb_vertex, float weight, double time);
void clear_instance_file(string instance_name);
#endif // !SOLVER_TO_CSV

