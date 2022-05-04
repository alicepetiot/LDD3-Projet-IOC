#include "../src/solver_to_csv.hpp"

/**
 * Methode servant a génerer un csv des solutions de chaque solver
 * Ajoute une ligne a un fichier "instance_XXXXXXXX.csv" dans le dossier csv
 * (Si celui-ci n'existe pas, il est créé)
 * 
 * @param intance_name : nom de l'instance
 * @param display : texte explicant comment est obtenue la solution
 * @param solver : le solver qui donne la solution (nombre de vertex + poids)
 * @param time : temps (en secondes) pour obtenir la solution
 */
void tocsv_solve_stats(string instance_name, string display, WeightedMaximumStableSolver *solver, std::chrono::duration<double> time)
{
  ofstream file;

  //Append the line at the end
  file.open("../csv/instance_" + instance_name + ".csv", std::ios_base::app);

  if(file.bad()) cout << "Erreur d'ouverture du fichier csv" << endl;
  else
  {
    //Infos :
    file 
      << display << " , " // Texte explicatif de la solution
      << solver->getCurrentNB() << " , " // Nombre de sommet 
      << solver->getCurrentCost() << " , " // Poids
      << time.count()  // Temps d'execution
      << endl;
  }

  file.close();
}

void tocsv_custom_stats(string instance_name, string display, int nb_vertex, float weight, double time)
{
  ofstream file;

  //Append the line at the end
  file.open("../csv/instance_" + instance_name + ".csv", std::ios_base::app);

  if(file.bad()) cout << "Erreur d'ouverture du fichier csv" << endl;
  else
  {
    //Infos :
    file 
      << display << " , " // Texte explicatif de la solution
      << nb_vertex << " , " // Nombre de sommet 
      << weight << " , " // Poids
      << time  // Temps d'execution
      << endl;
  }

  file.close();
}

void clear_instance_file(string instance_name)
{
  ofstream file;
  file.open("../csv/instance_" + instance_name + ".csv", std::ios_base::trunc);
  if(file.bad()) cout << "Erreur d'ouverture du fichier csv" << endl;
  //NB : COMMENTER LA LIGNE SUIVANTE SUPPRIME LA PREMIERE LIGNE DES FICHIER 
  else file << "Algorithme , Nombre de sommets , Poids de la solution , Temps d'execution" << endl; 
  file.close();
}