# Importation des librairies nécéssaires à la bonne exécution du code
import os # Permet d'accéder aux chemins des fichiers de l'ordinateur qui exécute le programme
import random # Bibliothèque permettant de générer des nombres flottants aléatoires (utile).
              # Attention, le nombre n'est pas purement aléatoire, il dépend d'une seed générée au moment de
              # l'exécution du code
import re

output_folder = "../weights"
if not (os.path.isdir(f'./{output_folder}')):
    os.mkdir(output_folder)
list_of_files = os.listdir('../instances')


for file in list_of_files:
    if(file[-1] is not 'b'):
        with open("../instances/" + file, 'r') as fp:
            if (file[-1] != "b"):
                for line in fp:
                    if line[0] == 'p':
                        words = re.split(r"\s{1,}", line)
                        content = int(words[2])
                        break
                        # La variable content contient désormais le nombre de noeuds du graphe

                with open(f"./{output_folder}/" + file, 'w') as output:
                    # On écrit au début du fichier le nombre de noeuds pour l'utilisateur
                    output.write(f"{content} ")

                    # Pour chaque noeud, on lui attribue un poids aléatoire entre 0 et 1
                    for i in range(content):
                        output.write(f"{random.random():.4f} ") # .4f permet de ne garder que 4 chiffres significatifs

                    # Nous n'oublions pas de fermer le fichier
                    output.close()
                    fp.close()
                    print(f"Generation completed for file {file}") # Message de terminaison pour l'utilisateur :)
            

# --- FIN DU SCRIPT ---