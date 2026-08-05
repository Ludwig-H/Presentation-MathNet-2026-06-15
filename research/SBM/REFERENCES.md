# Références et correspondance avec le dépôt

## 1. Manuscrit

- [Manuscrit de thèse](../../Manuscrit_de_thèse.pdf#page=133), chapitre 11,
  pages imprimées 123–149 :
  - §§11.1.1–11.1.4 : modèle bayésien, algorithmes, random guess et trois
    régimes de recouvrement ;
  - §11.1.5 : rééchantillonnage postérieur et couplages invariants ;
  - §11.2 : mesures de Gibbs signées et dynamiques de clusters ;
  - §11.3.1 : percolation nécessaire et borne sur la fraction recouvrable.
- [Source LaTeX du chapitre 11](../../ChapII.tex), utile pour retrouver les
  labels exacts des définitions, théorèmes et équations cités dans la note 01.
- Le dendrogramme n'est pas introduit au chapitre 11. Son langage vient des
  définitions de Single-Linkage et de dendrogramme des chapitres 2–3,
  notamment les pages imprimées 20–21 et 33–36.

Le présent dossier est donc une extension fidèle à la méthode du chapitre
11, pas la transcription d'un résultat déjà contenu dans le manuscrit.

## 2. Seuil de weak recovery du SBM

- Elchanan Mossel, Joe Neeman et Allan Sly,
  [*Stochastic Block Models and Reconstruction*](https://arxiv.org/abs/1202.1499),
  [DOI](https://doi.org/10.1007/s00440-014-0576-6). Impossibilité sous le
  seuil et lien avec la reconstruction sur arbre.
- Elchanan Mossel, Joe Neeman et Allan Sly,
  [*A Proof of the Block Model Threshold Conjecture*](https://arxiv.org/abs/1311.4115),
  [DOI](https://doi.org/10.1007/s00493-016-3238-8). Achievability efficace
  au-dessus de Kesten--Stigum.
- Laurent Massoulié,
  [*Community Detection Thresholds and the Weak Ramanujan Property*](https://arxiv.org/abs/1311.3085),
  [DOI](https://doi.org/10.1145/2591796.2591857). Preuve spectrale
  indépendante de l'achievability.

## 3. Reconstruction sur arbre et percolation d'information

- William Evans, Claire Kenyon, Yuval Peres et Leonard Schulman,
  [*Broadcasting on Trees and the Ising Model*](https://doi.org/10.1214/aoap/1019487349).
  Seuil de reconstruction du canal binaire et méthode de second moment.
- Emmanuel Abbe et Enric Boix-Adserà,
  [*An Information-Percolation Bound for Spin Synchronization on General Graphs*](https://arxiv.org/abs/1806.03227),
  [DOI](https://doi.org/10.1214/19-AAP1523). Contraction $\chi^2$ et
  majoration par percolation d'information.

## 4. Almost exact et exact recovery

- Yury Polyanskiy et Yihong Wu,
  [*Application of the Information-Percolation Method to Reconstruction Problems on Graphs*](https://arxiv.org/abs/1806.04195),
  [DOI](https://doi.org/10.4171/MSL/10). La méthode d'information-percolation
  en version KL, de la même famille que la borne $`\chi^2`$ d'Abbe--Boix ;
  citée par la comparaison de la [note 08](08_PREUVES_COMPLETES_SEUILS.md).
- Elchanan Mossel, Joe Neeman et Allan Sly,
  [*Consistency Thresholds for the Planted Bisection Model*](https://arxiv.org/abs/1407.1591),
  [DOI](https://doi.org/10.1214/16-EJP4185). Seuils d'almost exact
  ($`\lambda_n\to\infty`$) et d'exact recovery pour la bisection
  plantée — l'emprunt E3b de la partie II de la
  [note 08](08_PREUVES_COMPLETES_SEUILS.md).
- Anderson Y. Zhang et Harrison H. Zhou,
  [*Minimax Rates of Community Detection in Stochastic Block Models*](https://arxiv.org/abs/1507.05313),
  [DOI](https://doi.org/10.1214/15-AOS1428). Exposant minimax de la
  proportion d'erreurs.
- Emmanuel Abbe, Afonso Bandeira et Georgina Hall,
  [*Exact Recovery in the Stochastic Block Model*](https://arxiv.org/abs/1405.3267),
  [DOI](https://doi.org/10.1109/TIT.2015.2490670). Frontière
  $(\sqrt A-\sqrt B)^2=2$ dans le SBM symétrique logarithmique.
- Emmanuel Abbe et Colin Sandon,
  [*Community Detection in General Stochastic Block Models: Fundamental Limits and Efficient Recovery Algorithms*](https://arxiv.org/abs/1503.00609),
  [DOI](https://doi.org/10.1109/FOCS.2015.47). Divergence
  Chernoff--Hellinger et seuil général d'exact recovery.
- Emmanuel Abbe et Colin Sandon,
  [*Recovering Communities in the General Stochastic Block Model without Knowing the Parameters*](https://arxiv.org/abs/1506.03729).
  Algorithmes almost exact et exact dans les régimes de degré croissant.

## 5. Dynamiques et ligne de Nishimori

- Robert Edwards et Alan Sokal,
  [*Generalization of the Fortuin--Kasteleyn--Swendsen--Wang Representation and Monte Carlo Algorithm*](https://doi.org/10.1103/PhysRevD.38.2009).
- Robert Swendsen et Jian-Sheng Wang,
  [*Nonuniversal Critical Dynamics in Monte Carlo Simulations*](https://doi.org/10.1103/PhysRevLett.58.86).
- Yukito Iba,
  [*The Nishimori Line and Bayesian Statistics*](https://arxiv.org/abs/cond-mat/9809190),
  [DOI](https://doi.org/10.1088/0305-4470/32/21/302).
- Roy Glauber,
  [*Time-Dependent Statistics of the Ising Model*](https://doi.org/10.1063/1.1703954).

## 6. Notes techniques déjà présentes dans ce dépôt

Les notes suivantes contiennent les dérivations longues dont ce dossier
donne la version pédagogique :

- [baseline du chapitre 11](../hierarchical-swendsen-wang/foundations/02_CHAPTER_11_BASELINE.md) ;
- [critère hiérarchique à deux répliques](../hierarchical-swendsen-wang/foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) ;
- [cadre exact des horloges et heat baths](../hierarchical-swendsen-wang/foundations/01_MATHEMATICAL_FRAMEWORK.md) ;
- [pilote SBM et calcul de Kesten--Stigum](../hierarchical-swendsen-wang/active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) ;
- [double Gibbs répliqué](../hierarchical-swendsen-wang/active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) ;
- [port global et trois régimes de recovery](../hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md).

En cas de divergence entre une ancienne feuille de route et une identité
présentée ici, le
[statut canonique du projet](../hierarchical-swendsen-wang/CURRENT_STATUS.md)
reste la référence pour le GSBM triangulaire. Le présent dossier ne revendique
aucune amélioration de son seuil.
