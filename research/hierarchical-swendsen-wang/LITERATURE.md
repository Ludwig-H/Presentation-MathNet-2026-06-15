# État de l'art ciblé

Recherche ciblée effectuée jusqu'au 18 juillet 2026. Elle ne constitue pas une preuve exhaustive de nouveauté.

## Positionnement prudent

La combinaison exacte suivante n'a pas été retrouvée dans les sources consultées :

1. horloges exponentielles dépendant de $`|W_e|`$ ;
2. filtration complète des composantes calculée par Kruskal/minimum spanning forest ;
3. loi jointe exacte de type Edwards–Sokal sur $(\sigma,D)$ ;
4. heat baths aux nœuds internes utilisant tous les liens entre les deux fils ;
5. emploi de cette hiérarchie pour borner ou caractériser la weak recovery.

Des morceaux importants existent séparément. La nouveauté défendable n'est donc ni « utiliser des horloges » ni « construire un MST », mais la **factorisation conditionnelle sur le dendrogramme et son usage comme géométrie de l'information postérieure**.

## Fondations : FK, Swendsen–Wang et Edwards–Sokal

1. [Fortuin–Kasteleyn, *On the random-cluster model I*](https://doi.org/10.1016/0031-8914(72)90045-6). Représentation random-cluster, fondation géométrique des algorithmes de clusters.
2. [Swendsen–Wang, *Nonuniversal critical dynamics in Monte Carlo simulations*](https://doi.org/10.1103/PhysRevLett.58.86). Algorithme de clusters classique.
3. [Edwards–Sokal, *Generalization of the Fortuin–Kasteleyn–Swendsen–Wang representation and Monte Carlo algorithm*](https://doi.org/10.1103/PhysRevD.38.2009). Mesure jointe spins–liens ; modèle méthodologique pour prouver la marginale de Gibbs.
4. [Kandel–Ben-Av–Domany, *Cluster dynamics for fully frustrated systems*](https://doi.org/10.1103/PhysRevLett.65.941). Dynamiques par cellules pour interactions frustrées ; antécédent direct des règles triangulaires d'ordre supérieur.

## Révélation progressive et algorithmes voisins

5. [Machta et al., *Invaded cluster algorithm for equilibrium critical points*](https://arxiv.org/abs/cond-mat/9507094). Révélation progressive de liens jusqu'à invasion/percolation. Très proche de l'idée d'horloges ordonnées, mais la cible est l'auto-ajustement au point critique, pas une conditionnelle de Gibbs sur un dendrogramme fixé.
6. [Hauseux–Soprano-Loto–Avrachenkov, *Higher-order Monte Carlo cluster dynamics for community detection in Euclidean graphs*](https://inria.hal.science/hal-05267074). Travail directement antérieur du projet : dynamique d'ordre supérieur, frustration et bornes de percolation.

- [Barbu–Zhu, *Generalizing Swendsen–Wang to sampling arbitrary posterior probabilities*](https://doi.org/10.1109/TPAMI.2005.161). Chaîne réversible de split/merge/relabel sur des partitions de graphe pour des postérieures générales. C'est un voisin important pour la généralisation de Swendsen--Wang, mais il ne construit ni filtration Kruskal/MSF ni heat bath indexé par les nœuds d'un dendrogramme.
- [Houdayer, *A cluster Monte Carlo algorithm for 2-dimensional spin glasses*](https://arxiv.org/abs/cond-mat/0101116), [DOI](https://doi.org/10.1007/PL00011151). Mouvement de clusters à deux répliques adapté aux verres de spins ; voisin naturel du point de vue de l'overlap, sans hiérarchie de coalescence.

## Dendrogrammes, minimum spanning forests et percolation proche-critique

7. [Gower–Ross, *Minimum Spanning Trees and Single Linkage Cluster Analysis*](https://doi.org/10.2307/2346439). Référence classique pour l'équivalence entre la hiérarchie single-linkage et l'ordre des arêtes d'un minimum spanning tree.
8. [Lyons–Peres–Schramm, *Minimal Spanning Forests*](https://arxiv.org/abs/math/0412263), [DOI](https://doi.org/10.1214/009117906000000269). Géométrie des MSF en volume infini et lien avec l'invasion percolation.
9. [Garban–Pete–Schramm, *The scaling limits of the Minimal Spanning Tree and Invasion Percolation in the plane*](https://arxiv.org/abs/1309.0269). Relie une construction planaire du MST à la percolation proche-critique ; c'est un guide naturel pour étudier la mesure des temps $`\beta_{ij}`$. Leur couplage triangulaire particulier ne s'applique toutefois pas directement aux horloges d'arêtes, possiblement non i.i.d., du présent modèle.

### Bande critique, pivots et sprinkling

- [Garban–Pete–Schramm, *The scaling limits of near-critical and dynamical percolation*](https://arxiv.org/abs/1305.5526), [DOI](https://doi.org/10.4171/JEMS/786). Construit la fenêtre proche-critique à partir des mesures pivotales ; c'est le cadre naturel pour la géométrie macroscopique des fusions autour de $`\beta_c`$.
- [Garban–Pete–Schramm, *Pivotal, cluster, and interface measures for critical planar percolation*](https://arxiv.org/abs/1008.1378), [DOI](https://doi.org/10.1090/S0894-0347-2013-00772-9). Fournit la mesure limite des pivots qui motive le flux pivotal pondéré du nouveau critère.
- [Smirnov–Werner, *Critical exponents for two-dimensional percolation*](https://arxiv.org/abs/math/0109120). Calcule les exposants d'arms pour la percolation critique de sites sur le réseau triangulaire. Ces exposants ne sont pas transférés automatiquement à la percolation par arêtes triangulaire du présent dossier ; le programme annulaire du fichier 23 vise d'abord RSW et quasi-multiplicativité sans invoquer cette universalité.
- [Nolin, *Near-critical percolation in two dimensions*](https://arxiv.org/abs/0711.4948), [DOI](https://doi.org/10.1214/EJP.v13-565). Présente les estimées de longueur de corrélation et les relations d'échelle nécessaires pour choisir une fenêtre $`\delta_n`$ autour du seuil.
- [Damron–Sapozhnikov, *Outlets of 2D invasion percolation and multiple-armed incipient infinite clusters*](https://arxiv.org/abs/0903.4496), [DOI](https://doi.org/10.1007/s00440-010-0274-y). Montre que les goulots d'invasion approchent le seuil tout en conservant une structure d'outlets ; cette distinction avertit que $`\beta_{ij}`$ ne doit pas être supposé concentré exactement en $`\beta_c`$ pour une paire ponctuelle.
- [Duminil-Copin–Raoufi–Tassion, *Sharp phase transition for the random-cluster and Potts models via decision trees*](https://arxiv.org/abs/1705.03104), [DOI](https://doi.org/10.4007/annals.2019.189.1.2). Donne décroissance exponentielle sous le seuil et densité supercritique sur les graphes transitifs pour le random-cluster $`q\ge1`$. Ces résultats justifient les dichotomies géométriques de référence, mais ne s'appliquent pas directement à la postérieure frustrée conditionnelle.
- [Aizenman–Barsky, *Sharpness of the Phase Transition in Percolation Models*](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-108/issue-3/Sharpness-of-the-phase-transition-in-percolation-models/cmp/1104116538.full), [DOI](https://doi.org/10.1007/BF01210689). Fournit la décroissance exponentielle sous-critique utilisée pour localiser le LCA d'une paire conditionnée à être connectée au seuil.
- [Köhler-Schindler–Tassion, *Crossing probabilities for planar percolation*](https://arxiv.org/abs/2011.04618). Forme robuste de RSW et extensions de volume fini ; combinée à la décroissance sous-critique, elle donne la séparation exponentielle/polynomiale du lemme de localisation critique.
- [Járai, *Incipient infinite percolation clusters in 2D*](https://projecteuclid.org/journals/annals-of-probability/volume-31/issue-1/Incipient-infinite-percolation-clusters-in-2D/10.1214/aop/1046294317.full), [DOI](https://doi.org/10.1214/aop/1046294317). Identifie plusieurs conditionnements critiques lointains à l'IIC en limite locale. Cette loi est un proxy pour le voisinage d'un endpoint ; elle ne fournit pas la loi Palm à deux points du backbone entier ni les buckets du dendrogramme de Kruskal.

**Contre-audit pour le chemin hiérarchique.** La mesure pivotale, la limite
d'échelle du MST et les outlets d'invasion ne donnent pas directement la loi
des tailles de buckets $`m_w`$ entre les deux fils successifs du dendrogramme.
L'exposant pivotal $`3/4`$ ne doit notamment pas être identifié au coefficient
$`\alpha`$ d'un éventuel régime $`m_w\sim\alpha\log H_L`$. Le nombre d'outlets
ne contrôle pas non plus $`N_{L,M}`$ : une arête record isolée peut produire
un bucket $`m_w=1`$, qui est un canal parfait dans l'oracle descendant. Le
verrou propre au projet est la loi jointe des variables
$`(H_L,(m_w,t_w)_w)`$ sous le conditionnement de paire lointaine dans la
géante et de LCA critique ; aucun des articles cités ci-dessus n'énonce ce
résultat.

## GSBM et synchronisation

10. [Abbe–Baccelli–Sankararaman, *Community Detection on Euclidean Random Graphs*](https://arxiv.org/abs/1706.09942). Modèle GSBM de référence et bornes de weak recovery.
11. [Saade–Krzakala–Lelarge–Zdeborová, *Spectral Detection in the Censored Block Model*](https://arxiv.org/abs/1502.00163). Seuil de reconstruction dans un modèle signé localement arborescent ; comparaison naturelle avec Kesten–Stigum.
12. [Abbe–Massoulié–Montanari–Sly–Srivastava, *Group Synchronization on Grids*](https://arxiv.org/abs/1706.08561). Récupération multiscale sur des grilles et lien avec la ligne de Nishimori.

## Information-percolation et SDPI

13. [Polyanskiy–Wu, *Application of Information-Percolation Method to Reconstruction Problems on Graphs*](https://arxiv.org/abs/1806.04195). Domination de l'information par une percolation dont les probabilités sont des contractions de canaux.
14. [Abbe–Boix, *An Information-Percolation Bound for Spin Synchronization on General Graphs*](https://arxiv.org/abs/1806.03227), [DOI](https://doi.org/10.1214/19-AAP1523). Formulation $\chi^2$ directement adaptée aux canaux binaires du GSBM ; combinée au seuil triangulaire exact de Wierman, elle donne la baseline $p<0.794659\ldots$.

14 bis. [Makur–Polyanskiy, *Comparison of Channels: Criteria for Domination by a Symmetric Channel*](https://arxiv.org/abs/1609.06877), [DOI](https://doi.org/10.1109/TIT.2018.2839743). Donne le critère less-noisy équivalent en termes de $`\chi^2`$. C'est le critère à fermer pour comparer le canal physique d'un triangle au canal d'effacement multi-état du fichier 11 ; une vérification sous le seul a priori uniforme ne suffit pas.

Dans l'échelle des horloges, cette baseline se réécrit $`t_\chi>\beta_c`$, où $`q_p(t_\chi)=(2p-1)^2`$. Après contraction de $`\Pi_{\beta_c}`$, elle utilise la sous-bande $`(\beta_c,t_\chi]`$, et non le graphe pur formé par toutes les arêtes de $`(\beta_c,1]`$.
15. [Gu–Polyanskiy, *Weak Recovery Threshold for the Hypergraph Stochastic Block Model*](https://arxiv.org/abs/2303.14689). SDPI multi-terminales : outil probable pour traiter simultanément tous les liens ou triangles traversant une fusion.
16. [Gu, *Exact reconstruction thresholds on hypertrees over a symmetric binary alphabet*](https://arxiv.org/abs/2606.21699). Prépublication du 19 juin 2026 : canaux BMS, comparaison de canaux et population dynamics rigoureuse. C'est un outil particulièrement proche du programme « loi complète du message $`B_u`$ » sur cactus de triangles ; il ne traite pas directement le dendrogramme de Kruskal.

## Reconstruction et capacités sur arbres

17. [Pemantle–Peres, *The Critical Ising Model on Trees, Concave Recursions and Nonlinear Capacity*](https://arxiv.org/abs/math/0503137). Critères de capacité pour des arbres non homogènes ; modèle conceptuel d'une capacité portée par les nœuds du dendrogramme.
18. [Evans–Kenyon–Peres–Schulman, *Broadcasting on Trees and the Ising Model*](https://doi.org/10.1214/aoap/1019487349). Reconstruction, flux et capacité électrique ; base des critères de transmission le long d'un arbre.
19. [Peres–Roch, *Reconstruction on Trees: Exponential Moment Bounds for Linear Estimators*](https://arxiv.org/abs/0908.2056). Contrôle quantitatif des estimateurs au-dessus de Kesten–Stigum ; utile pour comparer la chaîne LCA aux récursions de broadcast.

## Triangles, frustration et ligne de Nishimori

- [Wierman, *Bond percolation on honeycomb and triangular lattices*](https://doi.org/10.2307/1426685). Établit le seuil exact $`q_c=2\sin(\pi/18)`$ de la percolation par arêtes sur la grille triangulaire.
20. [Chayes–Lei, *Random Cluster Models on the Triangular Lattice*](https://arxiv.org/abs/cond-mat/0508254), [DOI](https://doi.org/10.1007/s10955-005-8078-7). Seuil autodual des états triangulaires corrélés utilisé au chapitre 11.
21. [Nishimori–Ohzeki, *Location of the Multicritical Point for the Ising Spin Glass on the Triangular and Hexagonal Lattices*](https://arxiv.org/abs/cond-mat/0601356), [DOI](https://doi.org/10.1143/JPSJ.75.034004). Condition principale face-vers-face et valeur conjecturée $`p_{\mathrm N}^{(0)}=0.8358058\ldots`$. Le fichier 13 réduit exactement leur équation (28) à $`H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1`$ bit ; cette simplification ne transforme pas l'ansatz principal en théorème de seuil.
22. [Ohzeki, *Locations of Multicritical Points for Spin Glasses on Regular Lattices*](https://arxiv.org/abs/0811.0464), [DOI](https://doi.org/10.1103/PhysRevE.79.021129). Remplace la condition principale de niveau zéro par des sommes partielles sur des amas suivant davantage le flot de renormalisation. La première approximation triangulaire donne $`0.835985`$ et montre pourquoi une hiérarchie ancestrale peut corriger, plutôt que reproduire à chaque profondeur, la valeur $`0.8358058`$.
23. [de Queiroz, *Multicritical Point of Ising Spin Glasses on Triangular and Honeycomb Lattices*](https://arxiv.org/abs/cond-mat/0510816), [DOI](https://doi.org/10.1103/PhysRevB.73.064410). Étude numérique par matrices de transfert ; classe critique distincte de la percolation ordinaire.
24. [Yamaguchi, *Conjectured Exact Percolation Thresholds of the FK Cluster for the ±J Ising Spin Glass*](https://arxiv.org/abs/1004.0654). Retrouve le seuil FK triangulaire correspondant à $`p_c^{\mathrm{edge}}`$.
25. [Fajen–Hartmann–Young, *Percolation of Fortuin–Kasteleyn Clusters for the Random-Bond Ising Model*](https://arxiv.org/abs/1905.04220), [DOI](https://doi.org/10.1103/PhysRevE.102.012131). Montre numériquement que, avec frustration, la percolation FK peut précéder l'ordre magnétique : raison précise pour laquelle la percolation seule n'est pas suffisante.

## Deux répliques et représentations graphiques

26. [Contucci–Giardinà–Nishimori, *Spin Glass Identities and the Nishimori Line*](https://arxiv.org/abs/0805.0754). Benchmark d'identités overlap–magnétisation sur la ligne de Nishimori. L'identité bayésienne plantée utilisée ici est redérivée dans le dossier ; cette référence ne traite pas directement le canal binaire discret ni le couplage MCMC présent.
27. [Aizenman, *Geometric Analysis of φ⁴ Fields and Ising Models. Parts I and II*](https://doi.org/10.1007/BF01205659). À champ nul dans l'Ising ferromagnétique, la représentation par courants aléatoires relie le carré d'une corrélation à la connexion dans l'union de deux courants indépendants sans sources. C'est un benchmark conceptuel pour le score LCA, pas une formule directement transférable au modèle frustré.
28. [Chayes–Machta–Redner, *Graphical Representations for Ising Systems in External Fields*](https://arxiv.org/abs/cond-mat/9806312), [DOI](https://doi.org/10.1023/B:JOSS.0000026726.43558.80). Représentation graphique à deux répliques dans laquelle l'ordre est caractérisé par une percolation, pour le cadre ferromagnétique traité.
29. [Machta–Newman–Stein, *The Percolation Signature of the Spin Glass Transition*](https://arxiv.org/abs/0707.0073), [DOI](https://doi.org/10.1007/s10955-007-9446-2). Relie overlap de répliques et géométrie de percolation dans les verres de spins.

## Heat baths, projections et censoring

- [Blackwell, *Equivalent Comparisons of Experiments*](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-24/issue-2/Equivalent-Comparisons-of-Experiments/10.1214/aoms/1177729032.full), [DOI](https://doi.org/10.1214/aoms/1177729032). Pour deux états, la domination des courbes ROC équivaut à une dégradation par noyau stochastique. Le fichier 19 l'applique à un bucket ; le fichier 20 tensorise les noyaux de dégradation sur un corridor fixé avec un prior latent arbitrairement corrélé.
- [Diaconis–Khare–Saloff-Coste, *Stochastic alternating projections*](https://projecteuclid.org/journals/illinois-journal-of-mathematics/volume-54/issue-3/Stochastic-alternating-projections/10.1215/ijm/1336568522.full), [DOI](https://doi.org/10.1215/ijm/1336568522). Relie explicitement Gibbs sampler, espérances conditionnelles et produits de projections ; c'est le cadre opératoriel de l'identité palindromique du fichier 19 et de l'enveloppe collapsed du fichier 20.
- [Dyer–Greenhill–Ullrich, *Structure and eigenvalues of heat-bath Markov chains*](https://arxiv.org/abs/1301.4055), [DOI](https://doi.org/10.1016/j.laa.2014.04.018). Caractérise les opérateurs de heat bath et leur positivité ; le produit systématique de plusieurs projections n'est toutefois pas lui-même une projection.
- [Peres–Winkler, *Can extra updates delay mixing?*](https://arxiv.org/abs/1112.0603), [DOI](https://doi.org/10.1007/s00220-013-1776-0). Inégalité de censoring pour systèmes de spins monotones démarrés d'un état extrémal. Le modèle frustré conditionnel du présent dossier ne satisfait pas automatiquement ces hypothèses.
- [Holroyd, *Some circumstances where extra updates can delay mixing*](https://arxiv.org/abs/1101.4690). Contre-exemples hors du cadre monotone ; ils imposent de certifier chaque programme hiérarchique au lieu de supposer qu'ajouter des heat baths améliore toujours la contraction.

## Déplacement strict par enhancement

- [Aizenman–Grimmett, *Strict monotonicity for critical points in percolation and ferromagnetic models*](https://doi.org/10.1007/BF01029985). Un enhancement local essentiel déplace strictement le seuil d'une percolation indépendante. C'est une piste pour transformer une vraie contraction multi-arêtes en gain strict, mais pas encore un théorème applicable aux signes frustrés ni aux coupes adaptatives de Kruskal.

## Lecture stratégique

Ordre recommandé :

1. Polyanskiy–Wu et Abbe–Boix pour fixer la borne d'information-percolation de référence dans ce dossier ;
2. Gu (2026) et Gu–Polyanskiy pour les canaux BMS et les contractions multi-terminales ;
3. Pemantle–Peres et Evans et al. pour convertir ces contractions en capacité ;
4. Gower–Ross, Nolin puis les trois articles de Garban–Pete–Schramm pour la géométrie LCA, les pivots et la MSF proche-critique ;
5. Aizenman, Chayes–Machta–Redner et Machta–Newman–Stein pour le lien deux répliques–connexion ;
6. Chayes–Lei, Nishimori–Ohzeki puis Ohzeki pour le cas triangulaire et les corrections de blocs.

La synthèse mathématique et les implications exactes pour la weak recovery sont centralisées dans [Bande critique et transmission de l'information](07_CRITICAL_BAND_CRITERION.md).
La stratégie de clôture actuellement privilégiée est détaillée dans
[Stratégie optimale pour une obstruction de weak recovery](23_OPTIMAL_WEAK_RECOVERY_OBSTRUCTION.md).

## Critère de nouveauté à viser

Un résultat qui ne ferait que retrouver la percolation de $`\Pi_1`$ serait déjà couvert par les représentations FK. Un résultat qui ne ferait que retrouver la contraction $(2p-1)^2$ serait couvert par information-percolation.

La contribution spécifique doit exploiter au moins un des éléments suivants :

- la valeur de $`\beta_u`$ ;
- l'ensemble des liens entre $`C_1`$ et $`C_2`$ ;
- une SDPI de bloc plus fine que le produit arête par arête ;
- le spectre ou la capacité de la matrice de persistance hiérarchique ;
- une caractérisation à deux répliques qui devienne nécessaire et suffisante dans un cas non trivial.
