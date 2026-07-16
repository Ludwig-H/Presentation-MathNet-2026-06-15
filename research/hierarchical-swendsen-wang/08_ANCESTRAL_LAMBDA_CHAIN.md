# Chaîne ancestrale des taux $`\Lambda_v`$

Ce fichier ferme le problème algébrique laissé dans le message extérieur $`B_u`$. Pour une paire $i,j$ et son nœud de fusion $`u=u_{ij}`$, il donne :

1. les quatre valeurs exactes de $`\Lambda_v(\sigma^{ab})`$ pour chaque ancêtre $`v\succ u`$ ;
2. une réduction de toute la chaîne ancestrale à quatre poids, puis à trois scalaires lorsqu'aucun poids ne s'annule ;
3. la loi annealed exacte de ces taux dans le GSBM homogène, conditionnellement au squelette de Kruskal non marqué ;
4. des estimateurs certifiés et un programme précis pour les cactus, les bandes et la grille triangulaire.

Le résultat important est une désintégration finie exacte : la difficulté se sépare en une géométrie de Kruskal et un canal de marques indépendant conditionnellement à cette géométrie. Cela ne donne pas encore le seuil exact de weak recovery sur toute la grille, mais supprime le verrou « comment calculer tous les $`\Lambda_v`$ au-dessus de $u$ ? ».

## 1. Le problème à deux points

Soit $`(G_L)`$ une exhaustion du graphe, et soient $`i_L,j_L`$ tels que

```math
d(i_L,j_L)\longrightarrow\infty.
```

Pour le dendrogramme $D_L$, posons

```math
u_L=\mathrm{LCA}_{D_L}(i_L,j_L),
\qquad
\beta_{i_Lj_L}=\beta_{u_L}.
```

Cette notation vaut sur l'événement où le LCA existe. Sinon, on pose $`\beta_{i_Lj_L}=+\infty`$ et $`\eta_{u_L}=0`$ par la convention inter-racines du score étendu.

Dans le GSBM homogène triangulaire,

```math
\mathbb P(\beta_{i_Lj_L}\le t)
=
\tau_{i_Lj_L}(q_p(t)),
\qquad
q_p(t)=p(1-e^{-u_pt}),
\qquad
u_p=\log\frac p{1-p}.
```

Le niveau géométrique critique est

```math
\beta_c(p)=q_p^{-1}(q_c),
\qquad
q_c=2\sin(\pi/18).
```

### Une restriction proche-critique n'est pas automatique

Le LCA ponctuel de deux sommets lointains n'est pas en général concentré en $`\beta_c`$. Sous la convergence supercritique usuelle,

```math
\tau_{i_Lj_L}(q_p(t))
\longrightarrow
\theta(q_p(t))^2,
\qquad t>\beta_c.
```

Conditionnellement à une connexion avant $1$, la fonction de répartition limite est donc

```math
\frac{\theta(q_p(t))^2}{\theta(q_p(1))^2},
```

qui n'est généralement pas un saut en $`\beta_c`$. La fenêtre proche-critique décrit les pivots qui créent la géométrie macroscopique ; elle ne contient pas nécessairement les attaches tardives de deux sommets ponctuels.

Une réduction à $`\beta_u-\beta_c=o(1)`$ exige donc un lemme de localisation supplémentaire, par exemple

```math
\lim_{\delta\downarrow0}
\limsup_{L\to\infty}
\mathbb E\left[
\eta_{u_L}
\mathbf1_{\{\beta_c+\delta<\beta_{u_L}\le1\}}
\right]
=0.
```

Cette propriété est fausse si $`\eta_{u_L}`$ est remplacé par $1$ et n'est pas démontrée pour le score informationnel. La formulation sûre conserve donc toute la mesure sur $`(\beta_c,1]`$, ou contracte d'abord les attaches locales avant de définir un temps de fusion macroscopique.

Même lorsque $`\beta_u`$ appartient à une fenêtre proche-critique, les temps de ses ancêtres satisfont

```math
\beta_u<\beta_{v_1}<\cdots<\beta_{v_h}\le1.
```

Il faut donc bien contrôler toute la chaîne au-dessus de $u$.

## 2. Convention indispensable sur le dendrogramme

Les formules de ce dossier emploient le **dendrogramme de partitions non marqué** : $D$ contient les deux composantes fusionnées et le temps $`\beta_v`$, mais l'identité de l'arête qui atteint ce minimum est marginalisée.

Pour une coupe $`E_v`$, sa densité contient alors

```math
\Lambda_v(\sigma)e^{-\beta_v\Lambda_v(\sigma)}.
```

Après compensation avec l'énergie, le facteur conditionnel est

```math
F_v(\Lambda_v)
=
\Lambda_v e^{(1-\beta_v)\Lambda_v}.
```

Si l'on enrichit au contraire $D$ par l'identité $`e_v`$ de l'arête gagnante, la densité contient

```math
|W_{e_v}|\,
\mathbf1_{\{e_v\text{ satisfaite}\}}
e^{-\beta_v\Lambda_v(\sigma)},
```

et non $`\Lambda_v e^{-\beta_v\Lambda_v}`$. Les deux états auxiliaires définissent deux heat baths différents. Ils ne doivent pas être mélangés. L'arête gagnante peut être utilisée pendant Kruskal pour calculer $`\beta_v`$, puis son identité doit être oubliée si l'on veut employer les formules présentes ici.

## 3. Décomposition exacte des coupes ancestrales

Fixons un nœud

```math
u:C_u=C_1\mathbin{\dot\cup}C_2
```

et sa chaîne d'ancêtres stricts

```math
u=v_0\prec v_1\prec\cdots\prec v_h.
```

Pour $`v=v_\ell`$, notons $`P_v`$ le fils de $v$ qui contient $C_u$, et $`S_v`$ son autre fils. Toute arête de $`E_v`$ possède une extrémité dans $`P_v`$ et l'autre dans $`S_v`$. Décomposons

```math
E_v
=
E_v^{(0)}\mathbin{\dot\cup}
E_v^{(1)}\mathbin{\dot\cup}
E_v^{(2)},
```

où :

- $`E_v^{(1)}`$ contient les arêtes dont l'extrémité dans $`P_v`$ appartient à $`C_1`$ ;
- $`E_v^{(2)}`$ contient celles dont cette extrémité appartient à $`C_2`$ ;
- $`E_v^{(0)}`$ contient celles dont cette extrémité appartient à $`P_v\setminus C_u`$.

Posons, pour $`r\in\{0,1,2\}`$,

```math
T_{v,r}
=
\sum_{e\in E_v^{(r)}}|W_e|,
\qquad
\lambda_{v,r}
=
\sum_{e\in E_v^{(r)}}
|W_e|\mathbf1_{\{e\text{ satisfaite par }\sigma\}}.
```

### Proposition — quatre taux ancestraux exacts

Pour $`a,b\in\{0,1\}`$,

```math
\boxed{
\Lambda_v(\sigma^{ab})
=
\lambda_{v,0}
+
\begin{cases}
\lambda_{v,1},&a=0,\\
T_{v,1}-\lambda_{v,1},&a=1,
\end{cases}
+
\begin{cases}
\lambda_{v,2},&b=0,\\
T_{v,2}-\lambda_{v,2},&b=1.
\end{cases}
}
```

Équivalemment, avec

```math
X_{v,r}=2\lambda_{v,r}-T_{v,r},
```

on a

```math
\Lambda_v^{ab}
=
\Lambda_v^{00}-aX_{v,1}-bX_{v,2}.
```

### Preuve

Un flip de $`C_r`$ inverse la satisfaction exactement pour les arêtes ayant une seule extrémité dans $`C_r`$. Dans le bucket $`E_v`$, ce sont précisément les arêtes de $`E_v^{(r)}`$. Une arête satisfaite devient insatisfaite et réciproquement ; son poids satisfait total passe donc de $`\lambda_{v,r}`$ à $`T_{v,r}-\lambda_{v,r}`$. Les groupes sont disjoints, ce qui donne la formule.

Au nœud $u$ lui-même, la règle est différente : chaque arête joint $`C_1`$ à $`C_2`$, donc

```math
\Lambda_u^{ab}
=
\begin{cases}
\Lambda_u,&a=b,\\
T_u-\Lambda_u,&a\ne b.
\end{cases}
```

Ces identités sont déterministes et valables pour des poids non homogènes.

## 4. Calcul exact du message $`B_u`$

Posons

```math
\phi_v(x)
=
\begin{cases}
\log x+(1-\beta_v)x,&x>0,\\
-\infty,&x=0.
\end{cases}
```

Pour les quatre états, définissons

```math
\Phi_u^{ab}
=
\log\mu_0(\sigma^{ab})
+
\sum_{v\succ u}\phi_v(\Lambda_v^{ab}).
```

Alors

```math
\boxed{
B_u
=
\mathrm{LSE}(\Phi_u^{00},\Phi_u^{11})
-
\mathrm{LSE}(\Phi_u^{01},\Phi_u^{10}),
}
```

où $`\mathrm{LSE}(x,y)=\log(e^x+e^y)`$. Le log-odds complet est

```math
\boxed{
L_u
=
B_u
+
\phi_u(\Lambda_u)
-
\phi_u(T_u-\Lambda_u).
}
```

Cette écriture en quatre états est la version sûre lorsque certains taux sont nuls. Elle se calcule en un passage sur la chaîne ancestrale, soit $`O(h)`$ opérations une fois les six statistiques $`(T_{v,r},\lambda_{v,r})`$ disponibles.

## 5. Réduction exacte à trois scalaires

Supposons dans cette section que les quatre taux $`\Lambda_v^{ab}`$ soient strictement positifs pour chaque ancêtre. Toute fonction des deux signes

```math
x=(-1)^a,
\qquad
y=(-1)^b
```

s'écrit de manière unique

```math
\Phi_u^{ab}
=
C_u+h_{u,1}x+h_{u,2}y+J_uxy.
```

Pour un ancêtre $v$, ses trois coefficients utiles sont

```math
h_{v,1}
=
\frac14
\log\frac{\Lambda_v^{00}\Lambda_v^{01}}
{\Lambda_v^{10}\Lambda_v^{11}}
+
\frac{1-\beta_v}{4}
(\Lambda_v^{00}+\Lambda_v^{01}
-\Lambda_v^{10}-\Lambda_v^{11}),
```

```math
h_{v,2}
=
\frac14
\log\frac{\Lambda_v^{00}\Lambda_v^{10}}
{\Lambda_v^{01}\Lambda_v^{11}}
+
\frac{1-\beta_v}{4}
(\Lambda_v^{00}+\Lambda_v^{10}
-\Lambda_v^{01}-\Lambda_v^{11}),
```

et

```math
\boxed{
J_v
=
\frac14
\log\frac{\Lambda_v^{00}\Lambda_v^{11}}
{\Lambda_v^{01}\Lambda_v^{10}}.
}
```

Le terme linéaire $`(1-\beta_v)\Lambda_v^{ab}`$ ne contribue pas à $`J_v`$ parce que $`\Lambda_v^{ab}`$ est affine en $(a,b)$. Le couplage ancestral direct des deux flips vient donc uniquement du préfacteur de course $`\log\Lambda_v^{ab}`$.

Les coefficients de $`\log\mu_0(\sigma^{ab})`$ s'ajoutent aux sommes

```math
h_1=\sum_{v\succ u}h_{v,1},
\qquad
h_2=\sum_{v\succ u}h_{v,2},
\qquad
J=\sum_{v\succ u}J_v
```

dans le cas uniforme. La sommation des deux orientations absolues donne alors

```math
\boxed{
B_u
=
2J
+
\log\cosh(h_1+h_2)
-
\log\cosh(h_1-h_2).
}
```

Cette formule montre deux mécanismes distincts :

1. un couplage direct $J$ créé par la marginalisation de l'arête gagnante ;
2. deux champs ancestraux $`h_1,h_2`$ qui produisent un message de parité après sommation des orientations absolues.

### Bornes déterministes et troncature

Comme $`\log\cosh`$ est $1$-Lipschitz,

```math
|B_u|
\le
2|J|+2\min(|h_1|,|h_2|).
```

Si la chaîne est tronquée et si $`\Delta h_1,\Delta h_2,\Delta J`$ désignent les coefficients omis, alors

```math
|B_u-B_u^{\mathrm{tronc}}|
\le
2(|\Delta h_1|+|\Delta h_2|+|\Delta J|).
```

Enfin, puisque $`\eta(L)=\tanh^2(L/2)`$,

```math
|\eta(L)-\eta(\widetilde L)|
\le
\frac{2}{3\sqrt3}|L-\widetilde L|.
```

On obtient donc un certificat a posteriori sur l'erreur de fiabilité induite par la troncature.

Le couplage d'un ancêtre admet aussi une borne locale. Si $`\lambda_{v,\min}`$ est le minimum des quatre taux et s'il est positif, alors, avec

```math
\Delta_{v,r}=T_{v,r}-2\lambda_{v,r},
```

la formule intégrale du second incrément de $`\log`$ donne

```math
|J_v|
\le
\frac{|\Delta_{v,1}\Delta_{v,2}|}
{4\lambda_{v,\min}^2}.
```

Un ancêtre dont au moins un des deux fils de $u$ possède une influence de bord négligeable crée donc peu de couplage direct.

## 6. Loi exacte des taux dans le GSBM homogène

Cette section est annealed et utilise l'identité de Nishimori. Conditionnellement à la réplique courante, chaque arête porte une marque indépendante

```math
Y_e
=
\mathbf1_{\{e\text{ satisfaite}\}}
\sim
\mathrm{Bernoulli}(p),
```

puis un temps d'activation

```math
T_e
=
\begin{cases}
\mathrm{Exp}(u_p),&Y_e=1,\\
+\infty,&Y_e=0.
\end{cases}
```

Notons $`\mathscr D`$ le squelette non marqué : arbre de partitions, buckets $`E_v`$ et temps $`\beta_v`$, mais ni les marques $`Y_e`$ des arêtes fermées ni l'identité de l'arête gagnante.

### Théorème de course conditionnelle — statut : établi

Conditionnellement à $`\mathscr D`$, les marques portées par des buckets distincts sont indépendantes. Pour un bucket $v$ de taille

```math
m_v=m_{v,0}+m_{v,1}+m_{v,2},
\qquad
m_{v,r}=|E_v^{(r)}|,
```

posons

```math
s_v
=
s_p(\beta_v)
=
\frac{pe^{-u_p\beta_v}}
{1-p+pe^{-u_p\beta_v}}
=
\mathrm{logistic}(u_p(1-\beta_v)).
```

Il existe une catégorie gagnante latente $`G_v\in\{0,1,2\}`$ telle que

```math
\mathbb P(G_v=r\mid\mathscr D)
=
\frac{m_{v,r}}{m_v}.
```

Conditionnellement à $`G_v`$, les trois comptes sont indépendants et

```math
\boxed{
K_{v,r}
:=
\frac{\lambda_{v,r}}{u_p}
\stackrel{d}{=}
\mathbf1_{\{G_v=r\}}
+
\mathrm{Bin}\left(
m_{v,r}-\mathbf1_{\{G_v=r\}},
s_v
\right).
}
```

Au nœud de fusion lui-même,

```math
\boxed{
K_u
:=
\frac{\Lambda_u}{u_p}
\stackrel{d}{=}
1+\mathrm{Bin}(m_u-1,s_p(\beta_u)).
}
```

### Preuve

Pour une arête,

```math
\mathbb P(T_e>t)
=
1-p+pe^{-u_pt}.
```

Conditionnellement à $`T_e>t`$, sa probabilité d'être satisfaite vaut $`s_p(t)`$. Dans un bucket de taille $m$, la densité du minimum non marqué en $t$ est

```math
m\,p\,u_pe^{-u_pt}
(1-p+pe^{-u_pt})^{m-1}.
```

L'arête qui réalise ce minimum est uniforme. Elle est satisfaite ; les autres arêtes sont conditionnées seulement par $`T_e>t`$ et gardent donc des marques Bernoulli indépendantes de paramètre $`s_p(t)`$.

Les buckets $`E_v`$ sont disjoints. La préimage d'un squelette fixé impose, dans chaque bucket, un minimum égal à $`\beta_v`$, et impose $`T_e>1`$ aux arêtes entre racines distinctes. Ces contraintes portent sur des ensembles d'arêtes disjoints, d'où l'indépendance conditionnelle entre buckets. La partition en trois groupes donne la formule annoncée.

### Ce que devient le biais de Kruskal

Le choix de Kruskal biaise la géométrie

```math
(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succeq u}
```

et corrèle les buckets géométriques le long de la chaîne. En revanche, conditionnellement à ce squelette complet, il n'existe plus de correction mystérieuse à apporter aux marques : leur noyau est exactement le produit précédent. C'est la distinction qui manquait dans les formulations « coupe fixée contre coupe sélectionnée ».

### Moyennes et concentration simultanée

On a exactement

```math
\mathbb E[K_{v,r}\mid\mathscr D]
=
m_{v,r}s_v
+
\frac{m_{v,r}}{m_v}(1-s_v).
```

Conditionnellement aux catégories gagnantes, une inégalité de Hoeffding et une union bound donnent, avec probabilité au moins $`1-\delta`$, simultanément pour les trois groupes des $h$ ancêtres,

```math
\left|
K_{v,r}
-
\mathbb E[K_{v,r}\mid\mathscr D,G_v]
\right|
\le
\sqrt{
\frac{m_{v,r}-\mathbf1_{\{G_v=r\}}}{2}
\log\frac{6h}{\delta}
}.
```

Ces intervalles se transportent exactement vers les quatre $`\Lambda_v^{ab}`$. Comme $`x\mapsto\phi_v(x)`$ est croissante sur $`(0,\infty)`$, l'arithmétique d'intervalles donne ensuite des bornes certifiées sur les quatre $`\Phi_u^{ab}`$, puis sur $`B_u,L_u`$ et $`\eta_u`$. Les configurations où une borne inférieure atteint zéro sont traitées avec le calcul en quatre états, sans division par un taux possiblement nul.

## 7. Désintégration exacte de la borne LCA

Définissons la fiabilité conditionnelle du nœud $u$ par

```math
\Gamma_u(\mathscr D)
=
\mathbb E[\eta_u\mid\mathscr D],
```

où l'espérance utilise le produit de lois binomiales groupées de la section précédente et le calcul exact des quatre poids ancestraux.

Dans le GSBM homogène,

```math
\boxed{
H_n^{\mathrm{LCA}}
=
\frac1n
+
\frac2{n^2}
\mathbb E_{\mathscr D}
\sum_{u\in\mathscr D}
|C_{u,1}||C_{u,2}|\,
\Gamma_u(\mathscr D).
}
```

La loi marginale de $`\mathscr D`$ est celle obtenue avec des temps d'arêtes indépendants de fonction de répartition $`q_p(t)`$. La borne se sépare donc exactement en :

1. une loi géométrique de percolation couplée dans le temps ;
2. un noyau de canal ancestral explicite conditionnellement à cette géométrie.

La version de bande est

```math
Q_n
\le
S_n(\beta)
+
\frac2{n^2}
\mathbb E_{\mathscr D}
\sum_{u:\,\beta<\beta_u\le1}
|C_{u,1}||C_{u,2}|\,
\Gamma_u(\mathscr D).
```

Cette identité est le bon point de départ pour chercher une borne de non-recouvrement strictement supérieure à $`0.794659\ldots`$.

## 8. Trois méthodes de calcul

### 8.1 Réalisation fixée — calcul exact

1. Construire le dendrogramme de partitions avec Kruskal, puis oublier les identités gagnantes.
2. Affecter chaque arête originale au LCA de ses extrémités.
3. Pour le nœud cible $u$, parcourir ses ancêtres et agréger les six statistiques $`(T_{v,r},\lambda_{v,r})`$.
4. Accumuler les quatre log-poids par `log-sum-exp`.
5. Ajouter le facteur local de $u$ et calculer $`m_u=\tanh(L_u/2)`$, puis $`\eta_u=m_u^2`$.

Pour une paire fixée, les buckets peuvent être remplis en parcourant les arêtes incidentes à $`C_1\cup C_2`$ et en les envoyant vers leur LCA. Une fois les agrégats disponibles, le message coûte $`O(h)`$.

### 8.2 Squelette fini — énumération exacte

Pour un cactus ou une bande courte, conditionner les tailles de groupes et les temps, énumérer les comptes binomiaux groupés indépendants, puis sommer $`\eta_u`$. Cette méthode est exponentielle dans le nombre brut de buckets, mais exacte et adaptée aux tests de petite taille.

Le script [computations/ancestral_lambda_chain.py](computations/ancestral_lambda_chain.py) implémente ce calcul sans dépendance externe. Les tests vérifient :

- les quatre taux après flip ;
- l'égalité entre la formule séparée et les sommes de parité directes ;
- la normalisation de la loi groupée ;
- sa marginale $`1+\mathrm{Bin}(m-1,s_p(\beta))`$ ;
- l'identité terminale $`\mathbb E\eta_u=1/m`$ à $`\beta=1`$ sans ancêtre.

### 8.3 Grandes chaînes — bornes certifiées

Pour de grands buckets :

1. remplacer les comptes par leurs intervalles de concentration simultanée ;
2. propager ces intervalles dans les quatre log-poids ;
3. tronquer les ancêtres dont le certificat sur $`|\Delta B_u|`$ est inférieur à la précision demandée ;
4. traiter exactement les petits buckets proches de zéro.

Cette stratégie donne une borne rigoureuse, contrairement à une substitution directe de $`\Lambda_v`$ par sa moyenne dans la fonction non linéaire $`\log\Lambda_v`$.

## 9. Stratégie vers des seuils de weak recovery

### Cactus de triangles — premier seuil exact accessible

Sur un cactus, le squelette ancestral possède une description récursive et les buckets ont une taille bornée. La loi complète des trois coefficients $`(h_1,h_2,J)`$ peut être propagée exactement. La contraction obtenue doit être comparée à la calibration non oracle

```math
\gamma_1=(2p-1)^2.
```

Une équation de reconstruction exacte ou un rayon spectral certifié est réaliste dans ce cadre.

### Bandes triangulaires — matrices de transfert certifiées

Pour une largeur fixée, l'état doit contenir :

- la connectivité de la frontière ;
- les temps de fusion actifs ;
- les quatre log-poids ancestraux, ou une discrétisation par intervalles de $`(h_1,h_2,J)`$.

Une matrice de transfert avec arithmétique d'intervalles peut produire des seuils numériques rigoureux. Il faut ensuite étudier leur monotonie en fonction de la largeur.

### Grille entière — deux verrous restants

La formule ancestrale ne suffit pas seule à identifier le seuil planaire exact. Il reste :

1. la loi proche-critique du squelette groupé
   $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)`$ ;
2. la cohérence signée après marginalisation du dendrogramme, c'est-à-dire le facteur $`\kappa_n`$ du critère de bande.

La bonne quantité à mesurer est désormais

```math
\frac2{n^2}
\mathbb E_{\mathscr D}
\sum_u
|C_{u,1}||C_{u,2}|\,
\Gamma_u(\mathscr D)\,
\delta_{\beta_u},
```

et non la seule distribution des $`\beta_u`$ ni la fiabilité locale obtenue en posant arbitrairement $`B_u=0`$.

Une seconde amélioration possible est un heat bath **collapsed** qui marginalise aussi les orientations le long du chemin $`u\leadsto\mathrm{racine}`$. Comme il conditionne sur moins d'information que le heat bath LCA simple, sa borne $L^2$ ne peut qu'être meilleure. Sur cactus et bandes, cette marginalisation se prête à une matrice de transfert ; sur la grille, elle demande un contrôle des cycles.

## 10. Statut exact des conclusions

- Décomposition des buckets en trois groupes et formule des quatre $`\Lambda_v^{ab}`$ : **établi, déterministe, volume fini**.
- Réduction à quatre log-poids et formule en $`(h_1,h_2,J)`$ : **établi lorsque les taux requis sont positifs ; calcul quatre états exact sinon**.
- Loi binomiale groupée conditionnelle au squelette non marqué dans le GSBM homogène : **établie annealed, volume fini**.
- Désintégration de $`H_n^{\mathrm{LCA}}`$ par $`\Gamma_u(\mathscr D)`$ : **établie conditionnellement à la finalisation formelle de la mesure jointe A1**.
- Concentration et certificats de troncature : **établis conditionnellement au squelette ; leur utilité asymptotique dépend de la géométrie des buckets**.
- Concentration de tous les LCA ponctuels dans une fenêtre proche-critique : **fausse pour la connectivité brute et non démontrée pour le poids informationnel**.
- Nouveau seuil strict sur la grille triangulaire : **à prouver**.
- Seuil exact $`p\simeq0.8358058`$ : **repère conjectural, non obtenu par les identités présentes**.

Le prochain calcul décisif est donc clair : évaluer $`\Gamma_u(\mathscr D)`$ exactement sur cactus, puis par intervalles sur des bandes de largeur croissante, avant toute extrapolation au plan.
