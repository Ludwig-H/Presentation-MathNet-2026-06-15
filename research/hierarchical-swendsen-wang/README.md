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

Le cas favorable d'une fusion exactement critique est résolu au seul niveau
local. Si $`B_u=0`$ et si le bucket critique contient $`m`$ arêtes, sa
fiabilité $`\Gamma_m^c`$ vérifie

```math
\Gamma_m^c(p_{\mathrm{SW}})=\frac1m,
\qquad
\Gamma_m^c(p)\longrightarrow1
\quad\text{pour tout }p>p_{\mathrm{SW}}.
```

Sa fenêtre exacte est $`p-p_{\mathrm{SW}}\asymp m^{-1/2}`$. Ce calcul ne
donne toutefois pas le seuil de weak recovery : la contribution d'une fenêtre
critique est au plus sa masse
$`S_n(\beta_c+\delta)-S_n(\beta_c)`$ parmi les paires. Une fiabilité
conditionnelle proche de $`1`$ sur une masse $`o(1)`$ reste invisible dans
l'overlap global. Surtout, les slides 31--33 montrent que le vrai heat bath
contient le produit de tous les facteurs ancestraux. La priorité est donc
d'estimer les quatre $`\Lambda_v(\sigma^{ab})`$ pour chaque $`v\succ u`$ sous
la loi du squelette vu depuis une paire lointaine critique.

## Socle de départ

Les points 1, 2 et 5 ci-dessous sont établis sous les hypothèses indiquées. Les points 3, 4 et 6 à 9 sont des résultats finis dont l'algèbre a été vérifiée et dont la preuve courte est donnée dans ce dossier ; les résultats qui utilisent la mesure jointe restent à intégrer dans une rédaction formelle complète avec A1.

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

Rédiger d'abord le théorème fini LCA et le théorème de réduction à la bande
critique. Le premier calcul nouveau doit ensuite porter sur la loi de
$`(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succ u}`$ vue depuis une paire
lointaine dont le LCA est critique : exactement sur cactus, puis par matrices
de transfert certifiées sur bandes triangulaires. Cette loi alimente le noyau
conditionnel exact du fichier 10. La contraction non oracle doit toujours être
calibrée à une arête par $`\gamma_1=(2p-1)^2`$ et la masse des fusions doit être
conservée séparément. Sur la grille homogène, l'objectif reste une zone
rigoureuse de non-recouvrement dépassant $p=0.794659\ldots$ ; le point
multicritique de Nishimori conjecturé $p\simeq0.8358058$ demeure beaucoup plus
ambitieux.

## Sources internes

- [Chapitre 11 canonique](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/e5a2f06b77a6f3ac5f2865b41ea65a3d0f7834f0/Manuscrit_de_these/Manuscrit%20these%20Louis%20Hauseux/PartIII/ChapII.tex).
- [Audit mathématique canonique](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/e5a2f06b77a6f3ac5f2865b41ea65a3d0f7834f0/AUDIT_MATHEMATIQUE.md).
- [Présentation du 16 juillet 2026](../../beamer-presentation-reunion-2026-07-16/).

Ce dossier ne modifie ni les slides ni le manuscrit. Il isole les calculs et les conjectures avant toute réintégration dans un texte principal.
