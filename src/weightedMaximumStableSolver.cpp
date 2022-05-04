#include "weightedMaximumStableSolver.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <random>
#include <time.h>
#include <chrono>

#include <time.h>


void WeightedMaximumStableSolver::initSolution(){
	bestCost=0;
	currentCost=0;
	currentSolution.clear();
	currentSolution.resize(nbVertex);
	currentSolution.shrink_to_fit();
	bestSolution.clear();
	bestSolution.resize(nbVertex);
	bestSolution.shrink_to_fit();
	for(Vertex row =0; row < nbVertex; row++) currentSolution[row]=0;
	updateBestSolution();
}

void WeightedMaximumStableSolver::updateBestSolution(){
	if (!checkSolution()) 
	{
		std::cout << "Erreur : stable non valide" << std::endl;
		return;
	}
	if (currentCost>bestCost){
		bestCost = currentCost;
		bestSolution = currentSolution;
	}
}

bool WeightedMaximumStableSolver::checkSolution() {
 	for(Vertex row1 =0; row1 < nbVertex; row1++){
 		if (!currentSolution[row1]) continue;
 	 	for(Vertex row2 =row1+1; row2 < nbVertex; row2++){
 	 		if (currentSolution[row2] && graph.isEdge(row1,row2)) return false;
 	 	}	
 	}
 	float temp=0;
 	for(Vertex r =0; r < nbVertex; r++){
 		if (currentSolution[r]) 
		 {
			 temp += weights[r];
		 }
 	}
 	if (temp!= currentCost) {
 	std::cout << "WARNING: different value, objective is " <<  currentCost << " different from the computed sum " << temp << std::endl;
 	currentCost=temp;
 	}
	return true;
}

void WeightedMaximumStableSolver::displayBestSolution() {
      std::cout << "Stable of weighted cost " << bestCost << " with  vertices : ";
      for(Vertex nod =0; nod < nbVertex; nod++){
        if (bestSolution[nod]) std::cout <<  nod << " ";
      }
      std::cout << std::endl;
}

void WeightedMaximumStableSolver::importGraphDIMACS( char *file){
	graph.importGraphDIMACS(file);
	nbVertex = graph.getNbVertices();
	bestCost=0;
	initSolution();	
}

void WeightedMaximumStableSolver::importWeights(char *file){
	char buffer[5000] = {};
	char* path = strcat(strcat(buffer,strdup("../weights/")),file);
	std::cout << path << std::endl;
	std::ifstream source(path);
	if(source.is_open()) {
		int nbWeigth;
		source >> nbWeigth;
		//nbVertex = nbWeigth;
		for (int i = 0; i < nbWeigth; i++) {
			float x; 
			source >> x; 
			weights.push_back(x);
		}
		source.close();
	} else {
		std::cout << "Impossible d'ouvrir le fichier !" << std::endl;
	}
	//std::cout << weights.size();
	//for (Vertex i = 0; i < weights.size(); i++) {
	//	std::cout << weights[i];
	//}
}

void WeightedMaximumStableSolver::solverGreedy(char triMode)
{
	solverGreedy(triMode, 'A', 'A');
}

/**
 * Adds the vertex to the solution, adding it's cost to total cost 
 * and incrementing the count of vertices in solution
 * 
 * @param v : vertex to add 
*/
void WeightedMaximumStableSolver::addVertexToSolution(Vertex v)
{
	currentSolution[v] = true;
	// currentCost += weights[v];
	// currentNB++;
}

/** 
 * Methode appellee pour obtenir la solution d'un maximum stable solver
 * 
 * @param triMode : comment trier la solution
 * 
 * 	0 : poids croissant 
 *  1 : poids decroissant
 *  2 : degre croissant
 *  3 : degre decroissant
 * 
 * @param egaMode : voir egalite
 * @param pickMode : voir pick
*/
void WeightedMaximumStableSolver::solverGreedy(char triMode, char egaMode, char pickMode)
{
	currentSolution.clear();
	currentSolution.resize(nbVertex);
	currentSolution.shrink_to_fit();
	currentCost = 0;
	currentNB = 0;

	vector<std::pair<Vertex,float>> sortedVertex;
	sortedVertex.resize(nbVertex);
	sortedVertex.shrink_to_fit();
	vector<std::pair<Vertex,Vertex>> tmp;
	tmp.resize(nbVertex);
  
	//on trie 
	auto sort_croissant = [](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) { return i.second < j.second ; };
	auto sort_decroissant = [](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second > j.second ;};
  
	std::srand(time(0));
	switch (triMode)
	{
		case '0': //poids croissant
		{
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, weights[i]);
			std::sort(tmp.begin(), tmp.end(), sort_croissant);
		} break;
    
		case '1': //poids decroissant
		{
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, weights[i]);
			std::sort(tmp.begin(), tmp.end(), sort_decroissant);
		} break;

		case '2': //deg croissant
		{
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(), sort_croissant);
		} break;
    
		case '3': //deg decroissant
		{
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(), sort_decroissant);
		} break;
	}

	for (int i = 0; i < nbVertex - 1; i++)
	{
		if (tmp[i].second == tmp[i+1].second) egalite(tmp[i], tmp[i+1], egaMode);
	}

	for (int i = 0; i < nbVertex; i++) sortedVertex[i] = make_pair(tmp[i].first, weights[tmp[i].first]);

	pick(sortedVertex, pickMode);
}

/** 
 * Appel lorsque deux elements sont egaux et a separer
 * 
 * @param i : premier element
 * @param j : deuxieme elemen
 * 
 * @param egaMode : comment gerer une egalite
 *  A : random
 *  B : Poids le plus grand
 *  C : Poids le plus petit
 *  D : Degre le plus grand
 *  E : Degre le plus petit
 */
void WeightedMaximumStableSolver::egalite(std::pair<Vertex, Vertex> &i, std::pair<Vertex, Vertex> &j, char egaMode)
{
	switch (egaMode)
	{
		case 'A': //Random
		{
			if(rand()%2) std::swap(i, j);
		} break;

		case 'B': // Poids le plus grand
		{
			if(weights[i.first] < weights[j.first]) std::swap(i, j);
		} break;

		case 'C': // Poids le plus petit
		{
			if(weights[i.first] > weights[j.first]) std::swap(i, j);
		} break;

		case 'D': // Degre le plus grand
		{
			if(graph.getDegree(i.first) < graph.getDegree(j.first)) std::swap(i, j);
		} break;

		case 'E': // Degre le plus petit
		{
			if(graph.getDegree(i.first) > graph.getDegree(j.first)) std::swap(i, j);
		} break;
		
		default: break;
	}
}

/**
 * Methode appellee pour obtenir la solution finale apres avoir trie nos sommets
 * 
 * @param sortedVertex 
 * 
 * @param pickMode
 * A : Prendre les 5 meilleurs jusqu'a ne plus pouvoir
 * B : Prendre les 5 meilleurs puis random
 * C : Prendre les 5 meilleurs puis random base sur le poids
 * D : Prendre les 5 meilleurs puis random base sur les degres
*/
void WeightedMaximumStableSolver::pick(vector<std::pair<Vertex,float>> sortedVertex, char pickMode)
{
	//inits
	bool consider = true;
	vector<Vertex> currentStable; currentStable.reserve(nbVertex);
	vector<Vertex> voisins; voisins.reserve(nbVertex);

	switch (pickMode)
	{
		case 'A': // Best first
		{
			for(int i = 0; i < nbVertex; i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}
				if (consider)
				{
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}
		} break;

		case 'B': // Best first then random (after 5 picks)
		{
			int picked = 0;
			//Pick the 5 bests
			for(int i = 0; i < nbVertex && picked < 5; i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
				}
				if (consider)
				{
					picked++;
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}
			////int cpt = 0;

			// Remove them from sorted thing
			std::remove_if(sortedVertex.begin(), sortedVertex.end(), [this](std::pair<Vertex, float> k) {return this->currentSolution[k.first];});

			//Shuffle the remaining vertices
			std::random_shuffle(sortedVertex.begin(), sortedVertex.end());
			
			//Pick the rest randomly
			for(int i = 0; i < sortedVertex.size(); i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}
				if (consider)
				{
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}
		} break;

		case 'C' : // Best first then random based on weights
		{
			int picked = 0;
			std::vector<std::pair<Vertex, float>> to_remove;
			//Pick the 5 bests
			for(int i = 0; i < nbVertex && picked < 5; i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}
				if (consider)
				{
					picked++;
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}

			//Remove them from sorted thing
			std::remove_if(sortedVertex.begin(), sortedVertex.end(), [this](std::pair<Vertex, float> k) {return this->currentSolution[k.first];});

			//Now we shuffle but we give a better chance to higher weights
			{
				unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
				unsigned int n = sortedVertex.end() - sortedVertex.begin();
				
				for (auto i= n - 1; i>0; --i) 
				{
					std::uniform_int_distribution<decltype(i)> d(0,i);
					auto g = std::default_random_engine(seed);
					int j = d(g);
					// On a 1 chance sur 2 de skip si les weights sont decroissantes
					if(weights[sortedVertex[i].first] < weights[sortedVertex[j].first] && i > n/2 && rand()%2 == 0) continue;
					if(weights[sortedVertex[i].first] > weights[sortedVertex[j].first] && i <= n/2 && rand()%2 == 0) continue;
					swap (sortedVertex[i], sortedVertex[j]);
				}
			}

			//Pick the rest randomly
			for(int i = 0; i < sortedVertex.size(); i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);
				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}

				if (consider)
				{
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}
		} break;

		case 'D' : // Best first then random based on degree
		{
			int picked = 0;
			std::vector<std::pair<Vertex, float>> to_remove;
			//Pick the 5 bests
			for(int i = 0; i < nbVertex && picked < 5; i++)
			{
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}
				if (consider)
				{
					picked++;
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}

			//Remove them from sorted thing
			std::remove_if(sortedVertex.begin(), sortedVertex.end(), [this](std::pair<Vertex, float> k) {return this->currentSolution[k.first];});

			//Now we shuffle but we give a better chance to higher weights
			{
				unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
				unsigned int n = sortedVertex.end() - sortedVertex.begin();
				
				for (unsigned int i = n - 1; i > 0; --i) 
				{
					std::uniform_int_distribution<unsigned int> d(0,i);
					auto g = std::default_random_engine(seed);
					int j = d(g);
					// On a 1 chance sur 2 de skip si les degree sont decroissants
					if(graph.getDegree(sortedVertex[i].first) < graph.getDegree(sortedVertex[j].first) && i > n/2 && rand()%2 == 0) continue;
					if(graph.getDegree(sortedVertex[i].first) > graph.getDegree(sortedVertex[j].first) && i <= n/2 && rand()%2 == 0) continue;
					swap (sortedVertex[i], sortedVertex[j]);
				}
			}

			//Pick the rest randomly
			for(int i = 0; i < sortedVertex.size(); i++)
			{
				
				
				consider = true;
				voisins = graph.getNeighbors(sortedVertex[i].first);

				if(voisins.size() != 0)
				{
					for (Vertex v : voisins)
					{
						if(currentSolution[v]) 
						{
							consider = false;
							break;
						}
					}
				}
				if (consider)
				{
					currentStable.push_back(sortedVertex[i].first);
					addVertexToSolution(sortedVertex[i].first);
				}
			}
		} break;
		
		default: break;
	}

	for(int i = 0; i < nbVertex; ++i)
		if(currentSolution[i])
		{
			currentCost += weights[i];
			currentNB++;
		}
}


void WeightedMaximumStableSolver::solverGreedyWeight(char triMode){
    currentSolution.clear();
	currentSolution.resize(nbVertex);
    currentCost = 0;
	currentNB = 0;
    vector<std::pair<Vertex,float>> sortedVertex;
    sortedVertex.resize(nbVertex);

	vector<std::pair<Vertex,Vertex>> tmp;
	tmp.resize(nbVertex);
    //on trie 
	switch (triMode){
		case '0': //par poids décroissant si égal par degré croissant
			//std::sort(sortedVertex.begin(), sortedVertex.end(), sortWeights);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertex[i] = make_pair (i, weights[i]);
			std::sort(sortedVertex.begin(), sortedVertex.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				//std::cout << i << std::endl;
				if (sortedVertex[i].second == sortedVertex[i+1].second){
					if (graph.getDegree(sortedVertex[i].first) > graph.getDegree(sortedVertex[i +1].first)){
						std::pair<Vertex, Vertex> tmp = sortedVertex[i];
						sortedVertex[i] = sortedVertex[i+1];
						sortedVertex[i+1] = tmp;
					}
				}
			}
		break;
		case '1': //par poids décroissant si égal par degré décroissant
			//std::sort(sortedVertex.begin(), sortedVertex.end(), sortWeights);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertex[i] = make_pair (i, weights[i]);
			std::sort(sortedVertex.begin(), sortedVertex.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				//std::cout << i << std::endl;
				if (sortedVertex[i].second == sortedVertex[i+1].second){
					if (graph.getDegree(sortedVertex[i].first) < graph.getDegree(sortedVertex[i +1].first)){
						std::pair<Vertex, Vertex> tmp = sortedVertex[i];
						sortedVertex[i] = sortedVertex[i+1];
						sortedVertex[i+1] = tmp;
					}
				}
			}
		break;
		case '2': //par poids croissant si égal par degré croissant
			//std::sort(sortedVertex.begin(), sortedVertex.end(), sortWeights);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertex[i] = make_pair (i, weights[i]);
			std::sort(sortedVertex.begin(), sortedVertex.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				//std::cout << i << std::endl;
				if (sortedVertex[i].second == sortedVertex[i+1].second){
					if (graph.getDegree(sortedVertex[i].first) > graph.getDegree(sortedVertex[i +1].first)){
						std::pair<Vertex, Vertex> tmp = sortedVertex[i];
						sortedVertex[i] = sortedVertex[i+1];
						sortedVertex[i+1] = tmp;
					}
				}
			}
		break;
		case '3': //par poids croissant si égal par degré décroissant
			//std::sort(sortedVertex.begin(), sortedVertex.end(), sortWeights);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertex[i] = make_pair (i, weights[i]);
			std::sort(sortedVertex.begin(), sortedVertex.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				//std::cout << i << std::endl;
				if (sortedVertex[i].second == sortedVertex[i+1].second){
					if (graph.getDegree(sortedVertex[i].first) < graph.getDegree(sortedVertex[i +1].first)){
						std::pair<Vertex, Vertex> tmp = sortedVertex[i];
						sortedVertex[i] = sortedVertex[i+1];
						sortedVertex[i+1] = tmp;
					}
				}
			}
		break;
		case '4': //par degré croissant si égal par poids décroissant
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(),
					[](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second < j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				if (tmp[i].second == tmp[i+1].second){
					if (weights[tmp[i].first] > weights[tmp[i +1].first]){
						std::pair<Vertex, Vertex> pair = tmp[i];
						tmp[i] = tmp[i+1];
						tmp[i+1] = pair;
					}
				}
			}
			for (int i = 0; i < nbVertex; i++){
				sortedVertex[i] = make_pair(tmp[i].first, weights[tmp[i].first]);
			}
		break;
		case '5': //par degré croissant si égal par poids croissant
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(),
					[](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second < j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				if (tmp[i].second == tmp[i+1].second){
					if (weights[tmp[i].first] < weights[tmp[i +1].first]){
						std::pair<Vertex, Vertex> pair = tmp[i];
						tmp[i] = tmp[i+1];
						tmp[i+1] = pair;
					}
				}
			}
			for (int i = 0; i < nbVertex; i++){
				sortedVertex[i] = make_pair(tmp[i].first, weights[tmp[i].first]);
			}
		break;
		case '6': //par degré décroissant si égal par poids décroissant
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(),
					[](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second > j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				if (tmp[i].second == tmp[i+1].second){
					if (weights[tmp[i].first] < weights[tmp[i +1].first]){
						std::pair<Vertex, Vertex> pair = tmp[i];
						tmp[i] = tmp[i+1];
						tmp[i+1] = pair;
					}
				}
			}
			for (int i = 0; i < nbVertex; i++){
				sortedVertex[i] = make_pair(tmp[i].first, weights[tmp[i].first]);
			}
		break;
		case '7': //par degré décroissant si égal par poids croissant
			for (int i = 0; i < nbVertex; i++) tmp[i] = make_pair(i, graph.getDegree(i));
			std::sort(tmp.begin(), tmp.end(),
					[](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second > j.second ;}
			);
			for (int i = 0; i < nbVertex - 1; i++){
				if (tmp[i].second == tmp[i+1].second){
					if (weights[tmp[i].first] > weights[tmp[i +1].first]){
						std::pair<Vertex, Vertex> pair = tmp[i];
						tmp[i] = tmp[i+1];
						tmp[i+1] = pair;
					}
				}
			}
			for (int i = 0; i < nbVertex; i++){
				sortedVertex[i] = make_pair(tmp[i].first, weights[tmp[i].first]);
			}
		break;

	}
	
	float totWeight = 0;
	bool consider = true;
	Vertex nbVertexInStable = 0;

	Vertex currentVertex = 0;
	Vertex firstpossible = sortedVertex[0].first;
	vector<Vertex> currentStable;
	currentStable.reserve(nbVertex);
	vector<Vertex> voisins;
	voisins.reserve(nbVertex);
	Vertex nbParcoured = 0;
	currentStable.push_back(firstpossible);
	totWeight += weights[firstpossible];
	currentSolution[firstpossible] = 1;
	nbVertexInStable++;
	nbParcoured = 1;
	voisins = graph.getNeighbors(firstpossible);
	/*std::cout << "voisins de " << firstpossible << " : [ ";
	for (int j = 0; j < voisins.size(); j++){
		std::cout << voisins[j] << " ; ";
	}
	std::cout << " ]" << std::endl;*/

	
	for(int i = 1; i < nbVertex; i++){
		voisins.clear();
		currentVertex = sortedVertex[i].first;
		consider = true;
		voisins = graph.getNeighbors(currentVertex);
		for (int j = 0; j < voisins.size(); j++){
			if (currentSolution[voisins[j]] == 1){
				consider = false; 
			}
		}
		if (consider){
			currentStable.push_back(currentVertex);
			currentSolution[currentVertex] = 1;
			totWeight += sortedVertex[i].second;
			nbVertexInStable++;
			/*std::cout << "voisins de " << currentVertex << " : [ ";
			for (int j = 0; j < voisins.size(); j++){
				std::cout << voisins[j] << " ; ";
			}
			std::cout << " ]" << std::endl;*/
		}
	}
	/*std::cout << " stable de fin: [ ";
	for (int i = 0; i < currentStable.size(); i++){
		std::cout << currentStable[i] << " ; ";
	}
	std::cout << " ]" << std::endl;*/
	currentCost = totWeight;
	currentNB = nbVertexInStable;
	//std::cout << currentCost << std::endl;
}

void WeightedMaximumStableSolver::solverDegNeighbors(){
	//tri par degré minimal
	currentSolution.clear();
    currentCost = 0;
	currentNB = 0;
    vector<std::pair<Vertex,Vertex>> sortedVertex;
    sortedVertex.resize(nbVertex);

	//paires + tri par degré minimal (donc tri croissant)
	for (int i = 0; i < nbVertex; i++) sortedVertex[i] = make_pair(i, graph.getDegree(i));
	std::sort(sortedVertex.begin(), sortedVertex.end(),
			[](std::pair<Vertex,Vertex> i, std::pair<Vertex,Vertex> j) {return i.second < j.second ;}
	);


	//toutes ces def à trier
	float totWeight = 0;
	Vertex currentVertex = 0;
	Vertex firstpossible = 0;
	vector<bool> currentStable;
	currentStable.reserve(nbVertex);
	vector<Vertex> voisins;
	voisins.reserve(nbVertex);
	Vertex nbParcoured = 0;

	vector<bool> parcouru;
	parcouru.resize(nbVertex);
	vector<Vertex> vNeighbors;
	vNeighbors.reserve(nbVertex);

	float totVoisins = 0;
	float perteMin = 0;
	Vertex sommetAPrendre = 0;

	while (nbParcoured < nbVertex){
		//le premier pas parcouru
		while (parcouru[(sortedVertex[firstpossible]).first]) firstpossible++;
		currentVertex = (sortedVertex[firstpossible]).first;

		voisins.clear();
		voisins = graph.getNeighbors(currentVertex);
		for (int i = 0; i < voisins.size(); i++){
			if (!(parcouru[voisins[i]]) ){
				totVoisins += weights[voisins[i]];
			}
		}
		perteMin = weights[currentVertex] - totVoisins;
		sommetAPrendre = currentVertex;

		for (int i = 0; i < voisins.size(); i++){
			if (!(parcouru[voisins[i]])){
				totVoisins = 0;
				vNeighbors.clear();
				vNeighbors = graph.getNeighbors(voisins[i]);
				for (int j = 0; j < vNeighbors.size(); j++){
					totVoisins += weights[vNeighbors[j]];	
				}
				if (weights[voisins[i]] - totVoisins > perteMin){
					perteMin = weights[voisins[i]] - totVoisins;
					sommetAPrendre = voisins[i];
				}
			}
		}
		voisins = graph.getNeighbors(sommetAPrendre);
		
		for (int j = 0; j < voisins.size(); j++){
			if (!(parcouru[voisins[j]])){
				parcouru[voisins[j]] = true;
				nbParcoured++;
			}
		}
		parcouru[sommetAPrendre] = true;
		nbParcoured++;
		currentSolution[sommetAPrendre] = 1;
		currentNB++;
		totWeight += weights[sommetAPrendre];
	}
	currentCost = totWeight;
/* 
	@Function : solverWeightDegree
	@Parameters : (int)mode
		- 0 : tri degree/weights max + id max 
		- 1 : tri degree/weights max + id min 
		- 2 : tri degree/weights min + id max 
		- 3 : tri degree/weights min + id min 
		- 4 : tri weights/degree max + id max 
		- 5 : tri weights/degree max + id min 
		- 6 : tri weights/degree min + id max 
		- 7 : tri weights/degree min + id min 

*/
}

void WeightedMaximumStableSolver::solverWeightDegree(char mode) {
    currentSolution.clear();
	currentSolution.resize(nbVertex);

    vector<std::pair<Vertex,float>> sortedVertex;
	sortedVertex.resize(nbVertex);
	sortedVertex.clear();

	switch(mode) {
		case '0': {
			
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i, float(graph.getDegree(i)/weights[i]));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);

			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;

			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first > j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '1': {

			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i, float(graph.getDegree(i)/weights[i]));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);


			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;


			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first < j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '2': {
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i, float(graph.getDegree(i)/weights[i]));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);



			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;

			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first > j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '3': {
			
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i, float(graph.getDegree(i)/weights[i]));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);

			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;

			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first < j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '4': {
			
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i,weights[i]/float(graph.getDegree(i)));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);



			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;
	

			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first > j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '5': {
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i,weights[i]/float(graph.getDegree(i)));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second > j.second ;}
			);


			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;


			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first < j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '6': {
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i,weights[i]/float(graph.getDegree(i)));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);


			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;

			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first > j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}
		case '7': {
			vector<std::pair<Vertex,float>> sortedVertexRatio;
			sortedVertexRatio.resize(nbVertex);
			for (Vertex i = 0; i < nbVertex; i++)  sortedVertexRatio[i] = make_pair (i,weights[i]/float(graph.getDegree(i)));

			
			std::sort(sortedVertexRatio.begin(),sortedVertexRatio.end(),
					[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.second < j.second ;}
			);



			
			int i = 0;
			int ratio = 0;
			int cpt = 0;
			vector<std::pair<Vertex,float>> equalCandidats;

			
			while(i < nbVertex) {
				std::pair<Vertex,float> p = sortedVertexRatio[i];
				if (equalCandidats.size() == 0) {
					ratio = p.second;
					equalCandidats.push_back(p);
					i++;
				} else {
					if (p.second == ratio) {
						equalCandidats.push_back(p);
						i++;
					} else {
						std::sort(equalCandidats.begin(),equalCandidats.end(),
						[](std::pair<Vertex,float> i, std::pair<Vertex,float> j) {return i.first < j.first ;}
						);

						for (int j = 0; j < equalCandidats.size(); j++) {
							sortedVertex.push_back(equalCandidats[j]);
						}

						equalCandidats.clear();
					}
				}
			}
			break;
		}

	}

	currentNB = 0;
	currentCost = 0;
	
	for (int i = 0; i < nbVertex; i++) {

		std::pair<Vertex,float> candidat = sortedVertex[i];

		vector<Vertex> neighbors = graph.getNeighbors(candidat.first);


		/*
		std::cout << "-----------------------------------------NEIGHBORS-----------------------------------------"  << std::endl;
		for (int i = 0; i < neighbors.size(); i++) {
			std::cout << "(" << candidat.first << "," << neighbors[i] << ")" << ";";
		}
		std::cout << std::endl;
		std::cout << "-----------------------------------------NEIGHBORS-----------------------------------------"  << std::endl;
		*/

		bool isSolution = true; 

		if (neighbors.size() != 0) {
			for (Vertex i = 0; i < neighbors.size(); i++) {
				if(currentSolution[neighbors[i]] == 1) {
					isSolution = false;
				}
			}
		} 

		if (isSolution) {
			currentSolution[candidat.first] = 1;
			currentCost += weights[candidat.first];
			currentNB += 1;
		}
	}

	/*
	std::cout << "-----------------------------------------CURRENT SOLUTION-----------------------------------------"  << std::endl;
	for (int i = 0; i < currentSolution.size(); i++) {
		std::cout << "(" << i << "," << currentSolution[i] << ")" << ";";
	}
	std::cout << std::endl;
	std::cout << "-----------------------------------------CURRENT SOLUTION-----------------------------------------"  << std::endl;

	
	std::cout << "-------------------------------------------SOMMETS TRIES-------------------------------------------"  << std::endl;
	for (int i = 0; i < sortedVertex.size(); i++) {
		std::cout << "(" << sortedVertex[i].first << "," << sortedVertex[i].second << ")" << ";";
	std::cout << std::endl;
	std::cout << "-----------------------------------------SOMMETS TRIES-----------------------------------------"  << std::endl;
	*/

}