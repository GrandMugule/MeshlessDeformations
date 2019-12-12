# MeshlessDeformations
Deformation algorithms for animations based on shape matching

Le projet contient quatre fichiers executables correspondant chacun à une animation

##Structure
Le projet contient un dossier animations avec quatre fichiers executables correspondant chacun à une animation : stretching, pulling, falling et falling_with_rebound.
Les fonctions utilisées par ces animations sont dans les fichiers du dossier utils

##Ligne de commande
Pour lancer une animation, il faut mettre en premier argument le fichier de la forme chargée (en extension .off)
L'utilisateur peut aussi préciser les paramètres alpha, beta, step, en écrivant par exemple --alpha 0.1. Il peut enfin rajouter l'extension de clustering avec la commande --clusters nombre_de_clusters.
Exemple de ligne de commande : 

'''
./stretching ../data/bunny.off --alpha 0.1 -- beta 0.5 --step 0.01 --clusters 5
'''


