# INF5171-243-TP1 - Potentiel électrique

On demande votre aide pour paralléliser un programme de simulation de particules chargées.

Une particule chargée crée un champs électrique autour d'elle. La particule attire les charges de signe opposé et éloigne les charges du même signe. Dans la simulation, on assume que les particules sont dans une solution ou le vide et peuvent se déplacer. On fait bouger les particules dans le temps en fonction des forces électromotrices. La méthode des particules s'applique à un grand éventail de situations, par exemple la simulation des fluides, des planètes et des atomes.

Voici la structure de l'algorithme:

* Création de particules avec leur état initial (position, vitesse et charge)
* Pour chaque pas de temps de la simulation
  * Calculer les forces et déplacer les particules
  * Calculer le champs électrique résultant
  * Sauvegarder le champs électrique sous forme d'une image

Le déplacement des particules se fait en intégrant l'accélération deux fois pendant une très courte durée `dt`.

* force = somme des forces sur la particule
* accélération = force/masse
* vitesse = vitesse + accélération * dt
* position = position + vitesse * dt

Comme la force dépend de la position, on doit calculer toutes les forces, puis appliquer le déplacement des particules. Pour améliorer la précision de la simulation quand les particules sont très proche et que les forces sont élevées, on applique le déplacement en divisant le pas de temps en `substeps`. On doit faire les pas de temps l'un à la suite de l'autre.

Votre mission est de paralléliser le code série fourni. Les étapes à réaliser sont les suivantes.

1. Profiler le code pour déterminer le temps passé dans chaque fonction
1. Paralléliser la simulation avec Threading Building Bloc
1. Réaliser un test pour vérifier le bon fonctionnement de l'algorithme parallèle
1. Réaliser un banc d'essai pour mesurer la mise à l'échelle de votre solution
1. Écrire un court rapport sous forme de présentation PowerPoint (4-5 diapos max) pour rassembler votre démarche et vos résultats.

Le programme principal `potential` prend en paramètre plusieurs options. Voici un résumé des principales options.

* -c : permet de modifier l'échelle de couleurs utilisée (voir fichiers dans le répertoire data)
* -i : nombre d'itérations de calcul (nombre de pas de temps de la simulation)
* -dt : durée d'un pas de temps
* -r : résolution de l'image de sortie (ex: 100 équivaut à une image 100x100)
* -e : numéro d'expérience (différents agencements de particules)
* -n : nombre de particules à générer pour la simulation (selon l'expérience)
* -o : gabarit des fichiers d'image en sortie, incluant un répertoire de sortie et le numéro d'itération dans le nom
* -p : engin parallèle (0: série, 1: tbb)
* -v : mode verbose

Voici quelques indications sur les fichiers source.

* `src/main.cpp` : point d'entrée du programme principal `potential`
* `src/experiments.cpp` : fonctions de génération de particules
* `src/particle.cpp` : modèle de particule
* `src/potential.cpp` : contient l'interface `IPotential` et l'implémentation série `PotentialSerial`
* `src/potentialparallel.cpp` : implémentation parallèle (à réaliser)
* `test/test_parallel.cpp` : test de la classe parallèle (à réaliser)
* `src/bench_potential.cpp` : banc d'essai de performance (à réaliser)

Vous pouvez générer par exemple une animation 4K@60fps de l'interaction de particules dans un arrangement périodique. Voici les étapes, incluant l'installation des dépendances (Ubuntu 24.04) et l'exécution du programme principal.

```
# Installation des dépendances et outils nécessaires
sudo apt install valgrind mplayer ffmpeg libeigen3-dev libpng-dev libtbb-dev catch2 kcachegrind

# Compilation en mode Release
cmake -G Ninja -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release

# Astuce: ajouter les programmes dans le $PATH
source build-release/env.sh

# Exécution de la simulation
# Vous avez le temps d'aller prendre un café ;-)
potential -i 2000 -dt 2e-9 -n 25 -r 2000 -e 1

# Compilation du vidéo
ffmpeg -framerate 30 -i 'results/potential-%06d.png' -c:v libx264 -pix_fmt yuv420p -y potential.mp4

# Affichage
mplayer potential.mp4
```

> Note: la complexité de l'algorithme de la simulation est élevée. Par exemple, le calcul des forces entre les particules est O(n^2), car pour chaque particule, on doit faire la somme des forces de toutes les autres particules sur celle-ci. En réalité, les particules éloignées peuvent être négligées. À l'aide d'un index de recherche spatial, on peut déterminer les voisins proches plus rapidement et réduire la complexité algorithmique et ainsi simuler des millions de particules. Pour l'exercice, nous utilisons la version force brute simple avec un nombre limité de particules (n < 1000).

### Étape 1 : Profilage

Quelle fonction utilise le plus de temps de calcul? Vous allez utiliser un profileur pour déterminer la proportion du temps des principales fonction.

Déterminez une commande `potential` qui effectue quelques secondes de calcul (environ 5 secondes), en variant au besoin

* le nombre de particules
* le nombre de pas de temps
* la résolution

Profilez la commande avec l'outil Callgrind, faisant parti de Valgrind. Le programme doit être compilé avec le mode `-DCMAKE_BUILD_TYPE=RelWithDebInfo`, pour compiler avec les optimisations et les symboles de débogage (sinon, le profileur ne peut pas obtenir le nom des fonctions correctement). Voici un exemple:

```
# Profiler le programme
# crée un fichier callgrind.out.PID
valgrind --tool=callgrind ma-commande ...

# Visualiser le profile (le plus récent par défaut)
kcachegrind
```

Ensuite, remplissez le tableau suivant concernant les principales fonctions de la simulation (classe Potential). Vous allez paralléliser ces fonctions à l'étape 2.

| Fonction                        | % de temps | Accélération attendue |
|-------------------------------  |------------|-----------------------|
| PotentialSerial::compute_field  |          % |                       |
| PotentialSerial::move_particles |          % |                       |
| PotentialSerial::save_solution  |          % |                       |
| Total                           |          % |                       |

Utilisez le nombre de processeur de l'ordinateur sur lequel vous travaillez pour prédire l'accélération attendue. Il s'agit d'appliquer la loi d'Amdahl pour la proportion du temps d'exécution de chaque fonction comme si elle était parfaitement parallélisable. Pour la ligne total, on effectue la somme de toutes les proportions des fonctions listées.

Autrement dit, si on arrive à parfaitement paralléliser la liste des fonctions identifiées et que l'on néglige tout surcoût (démarrage du calcul, synchronisation, débalancement, etc.), alors quelle est l'accélération que vous obtienderiez? Vous allez comparer à l'étape 4 l'accélération maximale théorique calculée ici et l'accélération réelle obtenue.

Selon ce tableau, déduisez la proportion du programme qui restera en série.

Éléments de votre rapport: commande utilisée, tableau de résultats, proportion série mesurée du programme.

### Étape 2 : Parallélisation

Le code série est réalisé dans la classe `PotentialSerial`.

Implémentez la parallélisation des trois fonctions principales de la classe `PotentialParallel` avec TBB dans le fichier `src/potentialparallel.cpp`.

Vous pouvez créer des classes, des fonctions ou des champs comme cela vous avez besoin. Vous pouvez utiliser les fonctions avec une classe, ou lambda. Par contre, il ne faut pas modifier les arguments des méthodes surchargées de l'interface `IPotential`. La fonctionnalité de simulation doit demeurer la même, mais s'exécuter plus rapidement en utilisant tous les processeurs disponibles de votre ordinateur.

Éléments de votre rapport: résumé du choix des fonctions TBB et de l'implémentation

### Étape 3 : Test fonctionnel

Dans le fichier `test/test_parallel.cpp` réalisez trois tests fonctionnels (un test pour chaque méthode de `PotentialParallel`) pour vérifier que la version parallèle de cette fonction produit le même résultat que la version série.

Le test doit êre réalisé avec [Catch2](https://catch2-temp.readthedocs.io/).

Élément de votre rapport: résumé de l'implémentation des tests

### Étape 4 : Banc d'essai

Dans le fichier `src/bench_potential.cpp` réalisez un banc d'essai. Voici la structure: 

* Mesurer le temps d'exécution série
* Mesurer le temps d'exécution parallèle de 1 à n processeurs (inclusivement)

Le nombre de processeurs doit être détecté à l'exécution (et non pas une constante codée en dur en fonction de votre ordinateur!).

Chaque mesure doit être répétée 10 fois pour obtenir une moyenne. Reprenez les paramètres de votre mesure à l'étape 1.

Votre programme doit produire un fichier de résultat au format texte `bench_potential.dat` contenant sur chaque ligne, le nombre de processeur utilisé, le temps d'exécution moyen (en nanoseconde) et l'accélération par rapport à la version série. On parle ici de la mise à l'échelle forte (taille de problème fixe). Voici un exemple du format, incluant une ligne d'entête.

```
# ncpu temps acceleration
1 1234 1.0
2 1234 1.0
```

N'utilisez pas de sous-processus. Faites plutôt appel directement aux méthodes de simulation.

Pour obtenir vos résultats finaux, assurez-vous de compiler en mode `Release`. Aussi, il ne faut pas que le programme produise des affichages fréquents dans la boucle de mesure, car ceux-ci peuvent causer un ralentissement important. Finalement, portez attention à fermer les applications inutiles pendant la prise de mesure, comme le navigateur, pour réduire les sources de perturbation de mesure.

Vous devez réaliser deux graphiques, soit le temps d'exécution en fonction du nombre de processeurs, et l'accélération en fonction du nombre de processeurs. Ajoutez vos graphiques à votre rapport. Indiquez un titre et identifiez les axes. Vous pouvez utiliser l'outil de votre choix pour réaliser vos graphiques (gnuplot, matplotlib, Microsoft Excel, OpenOffice Calc, etc).

Élément de votre rapport: graphiques et analyse commaire des résultats.

### Remise

 * Faire une archive du code de votre solution avec la commande `make remise` ou `ninja remise`
 * Remise sur TEAMS (code et rapport), une seule remise par équipe
 * Identifiez le ou les codes permanents dans le fichier source `potentialparallel.cpp` et le rapport
 * Remettez votre rapport au format PowerPoint ou PDF
 * Code: 70%
 * Rapport: 30%
 * Respect du style: pénalité max 10% (guide de style dans le fichier `.clang-format`)
 * Qualité (fuite mémoire, gestion d'erreur, avertissements, etc): pénalité max 10%
 * Total sur 100

Bon travail !

# Note sur les logiciels externes

Le code intègre les librairies et les programmes suivants.

* https://github.com/catchorg/Catch2
* https://www.nongnu.org/pngpp
* https://eigen.tuxfamily.org
