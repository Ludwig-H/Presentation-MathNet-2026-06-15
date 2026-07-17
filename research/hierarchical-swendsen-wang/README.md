# Swendsen–Wang hiérarchique — dossier de recherche

Ce dossier rassemble le programme théorique ouvert par la dynamique fondée sur les horloges exponentielles et le dendrogramme de Kruskal. L'objectif prioritaire est de remplacer la seule obstruction de percolation du chapitre 11 par une quantité qui mesure la transmission d'information à toutes les échelles, puis d'obtenir des conditions nécessaires et suffisantes de weak recovery.

## Question centrale

Pour $K=2$, quand l'observation $`(X_n,W_n)`$ contient-elle une information macroscopique sur $`\Sigma_n`$, c'est-à-dire quand existe-t-il un algorithme $`\tau_n`$ tel que
```math
\mathbb P\left[
\left|\frac1n\sum_{i=1}^n \Sigma_{n,i}\tau_{n,i}\right|\ge \varepsilon
\right]
```
reste strictement positive pour un $\varepsilon>0$ ? Avec l'a priori i.i.d. uniforme, cela équivaut à battre le random guess dans la définition du manuscrit.

La dynamique hiérarchique doit servir de **couplage invariant de deux répliques postérieures**. La taille des composantes à la coupe $t=1$ n'est alors que le premier cas d'une observable plus riche : la persistance de l'information sous des heat baths effectués à différents nœuds du dendrogramme.

La voie actuellement prioritaire suit chaque paire jusqu'à son nœud de coalescence
```math
u_{ij}=\mathrm{LCA}_D(i,j).
```
À ce nœud, la parité des quatre flips contrôle exactement la survie de $`\sigma_i\sigma_j`$. Cela donne un score de fusion $`\eta_u=\tanh^2(L_u/2)`$, puis une borne sommable en temps linéaire dans le nombre de nœuds :
```math
Q_n
\le
H_n^{\mathrm{LCA}}
=
\frac1{n^2}\mathbb E\left[
n+2\sum_u|C_{u,1}||C_{u,2}|\eta_u
\right]
\le
\frac1{n^2}\mathbb E\sum_{R\text{ racine}}|R|^2.
```

La dernière quantité est la borne percolative de Swendsen--Wang : le nouveau score en est donc un raffinement exact, fusion par fusion.

La nouvelle réduction à la bande critique fixe maintenant le rôle du seuil $\beta$. Si

```math
S_n(\beta)
=
\frac1{n^2}\mathbb E\sum_{C\in\Pi_\beta}|C|^2
\longrightarrow0,
```

alors

```math
Q_n
\le
S_n(\beta)
+
\mathcal M_n((\beta,1]),
```

où $`\mathcal M_n`$ somme les fusions au-dessus de $\beta$, pondérées par $`\eta_u`$. La weak recovery exige donc une masse macroscopique de **connexions informatives** nées dans la bande. La connectivité seule n'est pas suffisante : le critère exact se factorise en connexion quotient, fiabilité locale et cohérence signée après marginalisation de $D$.

La réduction favorable hiérarchique fixe désormais le bon oracle pour une
borne d'impossibilité : $`i,j`$ sont lointains, appartiennent au même arbre et
se séparent à la descente au niveau $`\beta_c`$, c'est-à-dire que leur LCA
fusionne au seuil de percolation. Les paires proches sont négligeables et les
racines distinctes ont un score LCA nul. Si l'on prouve en plus que cet oracle
critique domine les paires qui fusionnent plus tard, alors

```math
Q_L
\le
b_L+S_L(\beta_c-\varepsilon)
+\Gamma_{L,\varepsilon}^{\mathrm{fav}}
+\text{erreur de domination}.
```

Cette réduction ne suppose pas que les temps LCA réels se concentrent au
seuil. Elle remplace leur expérience quatre états par une expérience critique
plus informative. Le lemme de domination correspondant, nommé HF, est formulé
dans le fichier 12.

Le calcul d'une fusion exactement critique n'est résolu au seul niveau local
que lorsque $`B_u=0`$. Si le bucket critique contient $`m`$ arêtes, sa
fiabilité $`\Gamma_m^c`$ vérifie

```math
\Gamma_m^c(p_{\mathrm{SW}})=\frac1m,
\qquad
\Gamma_m^c(p)\longrightarrow1
\quad\text{pour tout }p>p_{\mathrm{SW}}.
```

Sa fenêtre exacte est $`p-p_{\mathrm{SW}}\asymp m^{-1/2}`$. Ce calcul local
est seulement un contre-audit. Les slides 31--33 montrent que le vrai heat
bath contient le produit de tous les facteurs ancestraux. La priorité est donc
d'estimer les quatre $`\Lambda_v(\sigma^{ab})`$ pour chaque $`v\succ u`$ sous
la loi du squelette vu depuis une paire lointaine critique, puis de prouver la
domination HF. Le fichier 10 donne le noyau exact des marques et un certificat
de queue ; le fichier 12 transporte cette erreur jusqu'à la fiabilité
$`\tanh^2(L_u/2)`$.

Le canal d'un triangle physique étudié dans le fichier 11 reste un calcul
auxiliaire. Il ne remplace ni la chaîne hiérarchique des $`\Lambda_v`$, ni le
biais de la paire critique, ni HF. Sa constante conditionnelle ne constitue
donc pas l'objectif prioritaire de ce dossier.

Le fichier 13 donne en revanche une calibration exacte directement liée aux
horloges. L'équation triangulaire de Nishimori--Ohzeki est identiquement

```math
H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1\ \text{bit},
```

et sa racine supérieure unique est
$`p_{\mathrm N}^{(0)}=0.835805792367\ldots`$. Cette entropie conditionnelle à
quatre états est l'espérance de la surprise du gagnant d'une course
exponentielle. Au LCA, la course analogue a pour taux les quatre poids
$`q_u^{ab}`$ et dépend de tous les $`\Lambda_v^{ab}`$ ancestraux. L'identité
retrouve donc exactement la constante de la conjecture au niveau d'une face,
mais ne prouve pas qu'elle est le seuil de weak recovery : il manque le pont
entre autodualité de blocs, fiabilité collapsed de la paire critique, fuite
d'information par $D$ et domination HF.

## Socle de départ

Les points 1, 2 et 5 ci-dessous sont établis sous les hypothèses indiquées. Les points 3, 4 et 6 à 12 rassemblent des résultats finis ou des audits conditionnels dont l'algèbre a été vérifiée et dont le statut précis est donné dans le dossier ; les résultats qui utilisent la mesure jointe restent à intégrer dans une rédaction formelle complète avec A1.

1. La coupe $t=1$ des horloges redonne exactement les liens de Swendsen–Wang.
2. Les heat baths des orientations globales des arbres redonnent la recoloration de Swendsen–Wang lorsque l'a priori est uniforme. Aux feuilles, on obtient le heat bath mono-site de Glauber ; un noyau de Metropolis–Hastings mono-site ciblant la même conditionnelle est une variante valide.
3. Pour deux répliques postérieures indépendantes $\sigma^{(1)},\sigma^{(2)}$, la non-disparition de
```math
   Q_n=\mathbb E\left\langle
   \left(\frac1n\sum_i\sigma_i^{(1)}\sigma_i^{(2)}\right)^2
   \right\rangle
```
   caractérise exactement la weak recovery au sens « avantage avec probabilité positive » dans le cas binaire symétrique.
4. Tout parcours hiérarchique invariant fournit une matrice de persistance $`H_S`$. Si $`\mathbb E[\lambda_{\max}(H_S)/n]\to0`$ pour un parcours $S$, la weak recovery est impossible. Pour Swendsen–Wang aux racines, $`H_S(i,j)=\mathbf 1_{\{i,j\text{ dans la même composante}\}}`$ : on retrouve l'obstruction du chapitre 11.
5. Sur la grille triangulaire homogène, la borne d'information-percolation déjà connue donne l'impossibilité pour
```math
   p<\frac{1+\sqrt{2\sin(\pi/18)}}2=0.794659\ldots,
```

   ce qui est plus fort que les bornes Swendsen–Wang $0.673648\ldots$ et triangulaire d'ordre supérieur $0.719224\ldots$. Toute nouvelle borne hiérarchique doit donc être comparée à $0.794659\ldots$, pas seulement à la borne du chapitre 11.
6. Pour une paire fixée, le noyau qui rafraîchit $D$ puis met à jour son LCA est positif et réversible. Ses autocorrélations $`A_{ij}^{(m)}`$ décroissent vers $`c_{ij}(O)^2`$ sous ergodicité. Le score à un pas est ainsi le premier terme d'une suite allant vers le critère exact à deux répliques.
7. À toute coupe $\beta$ telle que $`S_n(\beta)\to0`$, la weak recovery se réduit au score signé des paires dont le LCA naît dans $`(\beta,1]`$. Sur la grille triangulaire, l'information-percolation se réécrit $`t_\chi(p)>\beta_c(p)`$, avec $`q_p(t_\chi)=(2p-1)^2`$.
8. Pour un nœud de fusion $u$, tous les taux $`\Lambda_v(\sigma^{ab})`$ de ses ancêtres se calculent exactement à partir de trois groupes par bucket. Conditionnellement au squelette de Kruskal non marqué, leur loi pondérée, leurs moments et leur covariance sont explicites. Le verrou asymptotique restant est la loi géométrique de ces groupes le long de la chaîne ancestrale biaisée par la paire critique.
9. Pour une fusion locale au temps critique et sans message ancestral, les paramètres se simplifient en
```math
   h_c(p)=\frac{2(p-p_{\mathrm{SW}})}{1-q_c},
   \qquad
   a_c(p)=2\,\mathrm{artanh}(h_c(p)).
```
   Le bord de la grande coupe informative est exactement $`p_{\mathrm{SW}}`$,
   avec une limite gaussienne explicite dans la fenêtre $`m^{-1/2}`$. Ce bord
   est oracle et ne remplace ni la masse des paires ni la contraction non
   oracle.
10. Pour un triangle observé, le profil de contraction sous un a priori de
    masse $t$ est la fonction rationnelle exacte $`c_q(t)`$ du fichier 11. La
    SDPI globale est $`2q^2/(1+q^2)`$, ce qui réfute la borne scalaire naïve
    $`0.829491\ldots`$. Le candidat multi-état
    $`0.809909\ldots`$ est explicitement étiqueté conditionnel : il dépend
    encore de la positivité d'une matrice rationnelle $`3\times3`$ dans le
    secteur polarisé.
11. Pour une paire lointaine, la réduction favorable du fichier 12 donne une
    implication globale exacte sous HF. Le message tronqué aux $`K`$ premiers
    ancêtres approche la fiabilité complète avec l'erreur certifiée

    ```math
    \min\left(1,\frac{2\mathcal R_u^{(>K)}}{3\sqrt3}\right).
    ```

    Une nouvelle borne triangulaire exige donc la convergence du squelette
    critique, la sommabilité de cette queue et la domination HF.
12. L'équation (28) de Nishimori--Ohzeki se réduit exactement à
    $`3h_2(p)-h_2((1+(2p-1)^3)/2)=1`$. Une course conditionnelle de quatre
    horloges redonne cette entropie sans répliques. La hiérarchie autoduale
    proposée dans le fichier 13 conserve les $`K`$ premiers ancêtres du LCA :
    son niveau face $`K=0`$ donne $`0.835805792367\ldots`$, tandis que les
    niveaux $`K\ge1`$ sont des problèmes de blocs à construire et certifier,
    non des seuils déjà démontrés.

## Carte du dossier

- [01_MATHEMATICAL_FRAMEWORK.md](01_MATHEMATICAL_FRAMEWORK.md) : définition exacte de $D$, loi jointe de type Edwards–Sokal, règles de mise à jour et hypothèses.
- [02_CHAPTER_11_BASELINE.md](02_CHAPTER_11_BASELINE.md) : théorème $\theta^{\max}$ corrigé, portée réelle et baseline d'information-percolation.
- [03_HIERARCHICAL_WEAK_RECOVERY.md](03_HIERARCHICAL_WEAK_RECOVERY.md) : critère à deux répliques, obstruction $`H_S`$, capacité hiérarchique et conjectures.
- [04_TRIANGULAR_GSBM.md](04_TRIANGULAR_GSBM.md) : calculs explicites sur la grille triangulaire et objectifs numériques/théoriques.
- [05_PROOF_ROADMAP.md](05_PROOF_ROADMAP.md) : lemmes à démontrer, dépendances, cas tests et critères de succès.
- [06_LCA_SPIN_CORRELATION.md](06_LCA_SPIN_CORRELATION.md) : quatre événements de flip, formule exacte faisant intervenir $`\beta_u=\xi_{e_u}`$, borne LCA, chaîne pair-spécifique et programme triangulaire.
- [07_CRITICAL_BAND_CRITERION.md](07_CRITICAL_BAND_CRITERION.md) : réduction à la bande critique, distinction bande pure/quotient, temps informationnel $`t_\chi`$, flux pivotal et capacité de blocs.
- [08_ANCESTRAL_LAMBDA_CHAIN.md](08_ANCESTRAL_LAMBDA_CHAIN.md) : formule exacte des quatre $`\Lambda_v`$ au-dessus du LCA, réduction à $`(h_1,h_2,J)`$, loi conditionnelle de Kruskal et méthodes de calcul certifiées.
- [09_CRITICAL_MERGER_ORACLE.md](09_CRITICAL_MERGER_ORACLE.md) : résolution exacte de la fusion critique locale, fenêtre $`m^{-1/2}`$, sandwich des taux ancestraux et contre-audit de la masse des paires.
- [10_ANCESTRAL_LAMBDA_ESTIMATION.md](10_ANCESTRAL_LAMBDA_ESTIMATION.md) : problème central des slides 31--33, course pondérée exacte, moments des quatre taux, concentration, certificat de queue et formulation du verrou géométrique sous le biais d'une paire critique.
- [11_TRIANGLE_BLOCK_SDPI.md](11_TRIANGLE_BLOCK_SDPI.md) : profil SDPI exact du canal de triangle, échec du regroupement scalaire, canal d'effacement multi-état et candidat conditionnel $`0.809909\ldots`$ avec son lemme manquant explicite.
- [12_FAVORABLE_HIERARCHICAL_REDUCTION.md](12_FAVORABLE_HIERARCHICAL_REDUCTION.md) : réduction exacte aux paires lointaines du même arbre, oracle de séparation au seuil, lemme de domination HF et transport certifié de la queue ancestrale vers la weak recovery.
- [13_NISHIMORI_HIERARCHICAL_CLOCKS.md](13_NISHIMORI_HIERARCHICAL_CLOCKS.md) : réduction exacte de la conjecture triangulaire à une entropie de face, représentation par course exponentielle, entropie collapsed avec tous les ancêtres et hiérarchie autoduale à contrôler.
- [LITERATURE.md](LITERATURE.md) : état de l'art primaire, voisins conceptuels et positionnement prudent de la nouveauté.
- [references.bib](references.bib) : bibliographie ciblée et autonome.
- [computations/README.md](computations/README.md) : cahier des charges des calculs exacts et simulations à ajouter.

## Convention de statut

Chaque affirmation nouvelle doit porter l'un des statuts suivants :

- **Établi** : preuve complète disponible ou résultat primaire cité avec ses hypothèses.
- **Immédiat à formaliser** : conséquence courte des identités du dossier, mais rédaction détaillée encore absente.
- **À prouver** : énoncé précis et plausible avec une stratégie identifiée.
- **Conjecture** : cible de recherche ; aucune utilisation en aval comme si elle était démontrée.

## Premier objectif publiable

Rédiger d'abord le théorème fini LCA et le théorème de réduction favorable du
fichier 12. Le premier calcul nouveau doit porter sur la loi de
$`(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succ u}`$ vue depuis une paire
lointaine du même arbre dont le LCA est critique : exactement sur cactus, puis
par matrices de transfert certifiées sur bandes triangulaires. Cette loi
alimente le noyau conditionnel exact du fichier 10. Il faut alors fermer, dans
cet ordre, la convergence des premiers ancêtres, la sommabilité du certificat
$`\mathcal R_u`$, le contrôle des coins nuls et la domination HF entre les
expériences postcritique et critique. Sur la grille homogène, l'objectif reste
une zone rigoureuse de non-recouvrement dépassant $p=0.794659\ldots$ ; aucune
nouvelle constante n'est annoncée avant ces quatre preuves. Le point
multicritique de Nishimori conjecturé $p\simeq0.8358058$ est maintenant
retrouvé exactement comme zéro de l'entropie autoduale d'une face. Le premier
test nouveau doit calculer sur un cactus, avec les mêmes taux ancestraux, le
défaut autodual $`\Psi_1`$ et la fiabilité favorable
$`\Gamma_1^{\mathrm{fav}}`$ : leur comparaison teste le lemme de pont NH3 sans
confondre une calibration locale avec un seuil. Le fichier 11 est conservé
comme audit auxiliaire, sans priorité sur la chaîne hiérarchique.

## Sources internes

- [Chapitre 11 canonique](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/e5a2f06b77a6f3ac5f2865b41ea65a3d0f7834f0/Manuscrit_de_these/Manuscrit%20these%20Louis%20Hauseux/PartIII/ChapII.tex).
- [Audit mathématique canonique](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/e5a2f06b77a6f3ac5f2865b41ea65a3d0f7834f0/AUDIT_MATHEMATIQUE.md).
- [Présentation du 16 juillet 2026](../../beamer-presentation-reunion-2026-07-16/).

Ce dossier ne modifie ni les slides ni le manuscrit. Il isole les calculs et les conjectures avant toute réintégration dans un texte principal.
