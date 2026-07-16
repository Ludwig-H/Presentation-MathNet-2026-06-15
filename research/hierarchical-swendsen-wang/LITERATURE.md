# État de l'art ciblé

Recherche ciblée effectuée jusqu'au 16 juillet 2026. Elle ne constitue pas une preuve exhaustive de nouveauté.

## Positionnement prudent

La combinaison exacte suivante n'a pas été retrouvée dans les sources consultées :

1. horloges exponentielles dépendant de \(|W_e|\) ;
2. filtration complète des composantes calculée par Kruskal/minimum spanning forest ;
3. loi jointe exacte de type Edwards–Sokal sur \((\sigma,D)\) ;
4. heat baths aux nœuds internes utilisant tous les liens entre les deux fils ;
5. emploi de cette hiérarchie pour borner ou caractériser la weak recovery.

Des morceaux importants existent séparément. La nouveauté défendable n'est donc ni « utiliser des horloges » ni « construire un MST », mais la **factorisation conditionnelle sur le dendrogramme et son usage comme géométrie de l'information postérieure**.

## Fondations : FK, Swendsen–Wang et Edwards–Sokal

1. [Fortuin–Kasteleyn, *On the random-cluster model I*](https://doi.org/10.1016/0031-8914(72)90045-6). Représentation random-cluster, fondation géométrique des algorithmes de clusters.
2. [Swendsen–Wang, *Nonuniversal critical dynamics in Monte Carlo simulations*](https://doi.org/10.1103/PhysRevLett.58.86). Algorithme de clusters classique.
3. [Edwards–Sokal, *Generalization of the Fortuin–Kasteleyn–Swendsen–Wang representation and Monte Carlo algorithm*](https://doi.org/10.1103/PhysRevD.38.2009). Mesure jointe spins–liens ; modèle méthodologique pour prouver la marginale de Gibbs.
4. [Kandel–Ben-Av–Domany, *Cluster dynamics for fully frustrated systems*](https://doi.org/10.1103/PhysRevB.45.4700). Dynamiques par cellules pour interactions frustrées ; antécédent direct des règles triangulaires d'ordre supérieur.

## Révélation progressive et algorithmes voisins

5. [Machta et al., *Invaded cluster algorithm for equilibrium critical points*](https://arxiv.org/abs/cond-mat/9507094). Révélation progressive de liens jusqu'à invasion/percolation. Très proche de l'idée d'horloges ordonnées, mais la cible est l'auto-ajustement au point critique, pas une conditionnelle de Gibbs sur un dendrogramme fixé.
6. [Hauseux–Soprano-Loto–Avrachenkov, *Higher-order Monte Carlo cluster dynamics for community detection in Euclidean graphs*](https://inria.hal.science/hal-05267074). Travail directement antérieur du projet : dynamique d'ordre supérieur, frustration et bornes de percolation.

## GSBM et synchronisation

7. [Abbe–Baccelli–Sankararaman, *Community Detection on Euclidean Random Graphs*](https://arxiv.org/abs/1706.09942). Modèle GSBM de référence et bornes de weak recovery.
8. [Saade–Krzakala–Lelarge–Zdeborová, *Spectral Detection in the Censored Block Model*](https://arxiv.org/abs/1502.00163). Seuil de reconstruction dans un modèle signé localement arborescent ; comparaison naturelle avec Kesten–Stigum.
9. [Abbe–Massoulié–Montanari–Sly–Srivastava, *Group Synchronization on Grids*](https://arxiv.org/abs/1706.08561). Récupération multiscale sur des grilles et lien avec la ligne de Nishimori.

## Information-percolation et SDPI

10. [Polyanskiy–Wu, *Application of Information-Percolation Method to Reconstruction Problems on Graphs*](https://arxiv.org/abs/1806.04195). Domination de l'information par une percolation dont les probabilités sont des contractions de canaux.
11. [Abbe–Boix, *An Information-Percolation Bound for Spin Synchronization on General Graphs*](https://arxiv.org/abs/1806.03227), [DOI](https://doi.org/10.1214/19-AAP1523). Formulation \(\chi^2\) directement adaptée aux canaux binaires du GSBM ; donne la baseline \(p<0.794659\ldots\) sur la grille triangulaire homogène.
12. [Gu–Polyanskiy, *Weak Recovery Threshold for the Hypergraph Stochastic Block Model*](https://arxiv.org/abs/2303.14689). SDPI multi-terminales : outil probable pour traiter simultanément tous les liens ou triangles traversant une fusion.

## Reconstruction et capacités sur arbres

13. [Pemantle–Peres, *The Critical Ising Model on Trees, Concave Recursions and Nonlinear Capacity*](https://arxiv.org/abs/math/0503137). Critères de capacité pour des arbres non homogènes ; modèle conceptuel d'une capacité portée par les nœuds du dendrogramme.
14. [Evans–Kenyon–Peres–Schulman, *Broadcasting on Trees and the Ising Model*](https://doi.org/10.1214/aop/1019160259). Reconstruction, flux et capacité électrique ; base des critères de transmission le long d'un arbre.

## Triangles, frustration et ligne de Nishimori

15. [Chayes–Lei, *Random Cluster Models on the Triangular Lattice*](https://arxiv.org/abs/cond-mat/0508254), [DOI](https://doi.org/10.1007/s10955-005-8078-7). Seuil autodual des états triangulaires corrélés utilisé au chapitre 11.
16. [Nishimori–Ohzeki, *Location of the Multicritical Point for the Ising Spin Glass on the Triangular and Hexagonal Lattices*](https://arxiv.org/abs/cond-mat/0601356), [DOI](https://doi.org/10.1143/JPSJ.75.034004). Valeur conjecturée \(p_{\mathrm N}=0.8358058\ldots\).
17. [de Queiroz, *Multicritical Point of Ising Spin Glasses on Triangular and Honeycomb Lattices*](https://arxiv.org/abs/cond-mat/0510816), [DOI](https://doi.org/10.1103/PhysRevB.73.064410). Étude numérique par matrices de transfert ; classe critique distincte de la percolation ordinaire.
18. [Yamaguchi, *Conjectured Exact Percolation Thresholds of the FK Cluster for the ±J Ising Spin Glass*](https://arxiv.org/abs/1004.0654). Retrouve le seuil FK triangulaire correspondant à \(p_c^{\mathrm{edge}}\).
19. [Fajen–Hartmann–Young, *Percolation of Fortuin–Kasteleyn Clusters for the Random-Bond Ising Model*](https://arxiv.org/abs/1905.04220), [DOI](https://doi.org/10.1103/PhysRevE.102.012131). Montre numériquement que, avec frustration, la percolation FK peut précéder l'ordre magnétique : raison précise pour laquelle la percolation seule n'est pas suffisante.

## Deux répliques et représentations graphiques

20. [Chayes–Machta–Redner, *Graphical Representations for Ising Systems in External Fields*](https://arxiv.org/abs/cond-mat/9806312), [DOI](https://doi.org/10.1023/B:JOSS.0000026726.43558.80). Représentation graphique à deux répliques dans laquelle l'ordre est caractérisé par une percolation, pour le cadre ferromagnétique traité.
21. [Machta–Newman–Stein, *The Percolation Signature of the Spin Glass Transition*](https://arxiv.org/abs/0707.0073), [DOI](https://doi.org/10.1007/s10955-007-9446-2). Relie overlap de répliques et géométrie de percolation dans les verres de spins.

## Lecture stratégique

Ordre recommandé :

1. Polyanskiy–Wu et Abbe–Boix pour fixer la meilleure borne d'impossibilité existante ;
2. Gu–Polyanskiy pour définir une contraction de fusion multi-terminale ;
3. Pemantle–Peres et Evans et al. pour convertir ces contractions en capacité ;
4. Chayes–Machta–Redner et Machta–Newman–Stein pour relier cette capacité au recouvrement de deux répliques ;
5. Chayes–Lei puis Nishimori–Ohzeki pour le cas triangulaire.

## Critère de nouveauté à viser

Un résultat qui ne ferait que retrouver la percolation de \(\Pi_1\) serait déjà couvert par les représentations FK. Un résultat qui ne ferait que retrouver la contraction \((2p-1)^2\) serait couvert par information-percolation.

La contribution spécifique doit exploiter au moins un des éléments suivants :

- la valeur de \(\beta_u\) ;
- l'ensemble des liens entre \(C_1\) et \(C_2\) ;
- une SDPI de bloc plus fine que le produit arête par arête ;
- le spectre ou la capacité de la matrice de persistance hiérarchique ;
- une caractérisation à deux répliques qui devienne nécessaire et suffisante dans un cas non trivial.
