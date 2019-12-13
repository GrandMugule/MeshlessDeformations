# MeshlessDeformations
Deformation algorithms for animations based on shape matching


## Structure
Le projet contient un dossier animations avec quatre fichiers executables correspondant chacun à une animation : stretching, pulling, falling et falling_with_rebound.
Les fonctions utilisées par ces animations sont dans les fichiers du dossier utils

## Ligne de commande
Pour lancer une animation, il faut mettre en premier argument le fichier de la forme chargée (en extension .off)
L'utilisateur peut aussi préciser les paramètres alpha, beta, step (pour h) en écrivant par exemple --alpha 0.1.
Il peut aussi rajouter l'extension de clustering avec la commande --clusters nombre_de_clusters,
et décider de la méthode de shape matching avec --deformation suivi de r pour rigide, l pour linéaire et q pour quadratique.
Enfin, il peut ajuster les paramètres c_yield, c_creep et c_max pour les déformations plastiques.
Exemple de ligne de commande : 

    ./stretching ../data/bunny.off --alpha 0.1 --beta 0.5 --step 0.01 --clusters 5 --deformation l


## Animations

### Stretching
* Cliquer sur un point de la forme
* Appuyer sur x, y ou z (w sur un clavier azerty) pour sélectionner l'axe dans lequel on veut déplacer le point
* Appuyer sur les touches flèches de gauche ou flèche de droite pour déplacer le point
* Appuyer sur espace pour mettre à jour la forme
* Appuyer sur d pour lancer l'animation
* Appuyer sur s pour arrêter l'animation

### Pulling
* Cliquer sur un point de la forme
* Appuyer sur x, y ou z (w sur un clavier azerty) pour sélectionner l'axe dans lequel on veut déplacer le point
* Appuyer sur les touches flèches de gauche ou flèche de droite pour déplacer le point
* Appuyer sur d pour lancer l'animation
* Appuyer sur s pour arrêter l'animation

### Falling et falling with rebound
* Appuyer sur x, y ou z (w sur un clavier azerty) pour choisir l'axe de la chute
* Appuyer sur m (virgule pour les claviers azerty) pour faire monter la forme
* Appuyer sur c pour déclencher la chute


