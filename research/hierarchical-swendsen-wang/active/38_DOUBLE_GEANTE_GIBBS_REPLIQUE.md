# Cible prioritaire : Gibbs exact répliqué sur la double géante

**Statut : identité finie exacte, réduction géométrique et désintégration
scalaire fermées ; cible positive exacte $`\mathcal D_L^\times`$ isolée ;
embedding récursif des ports et contraction TRI2 encore ouverts ; aucun
nouveau seuil triangulaire revendiqué.**

## 1. Décision

Pour chaque $`p>1/2`$, calculer le niveau de percolation formel
$`\beta_c(p)`$. Il appartient à l'horizon physique $`[0,1]`$ exactement
lorsque le modèle final atteint la criticité ; sinon on pose
$`\beta_\star=1`$ et l'on enregistre explicitement l'absence de traversée
critique. Dans chaque arbre final du dendrogramme :

1. couper l'arbre à $`\beta_c(p)`$ ;
2. utiliser les sous-arbres critiques comme blocs d'élimination ;
3. conserver exactement tous les facteurs, ports, buckets et rangs situés
   au-dessus de la coupe ;
4. tirer le Gibbs joint de l'arbre entier ;
5. utiliser des hasards indépendants entre arbres finaux, seulement lorsque
   la factorisation conditionnelle le permet.

Pour étudier la weak recovery, cette construction doit être faite dans deux
répliques postérieures indépendantes. Chaque réplique possède son propre
dendrogramme complet. Partager le même dendrogramme entre les deux répliques
révèle une information auxiliaire et ne calcule pas le carré de la
corrélation postérieure.

La cible n'est donc pas un canal hiérarchique approché :

```math
\boxed{
\text{deux Gibbs exacts sur deux forêts entières}
\ \longrightarrow\
\text{deux coupes à }\beta_c
\ \longrightarrow\
\text{double géante}
\ \longrightarrow\
\text{opérateur overlap}.
}
\qquad\text{(1.1)}
```

Cette formulation pousse le principe du chapitre 11 au-delà du recoloriage
Swendsen--Wang : l'[invariance postérieure](../../../ChapII.tex) reste le
garde-fou, mais la quantité propagée est le Gibbs joint de l'arbre entier.

### 1.1 Ce qui dépasse effectivement le chapitre 11

Le chapitre 11 fournit le resampling postérieur, la balance locale et une
borne fondée sur **un objet gelé à un pas**, dont les clusters sont
recoloriés indépendamment. Ces éléments ne constituent ici que le socle.
Le programme nouveau porte sur quatre objets absents de cette preuve :

1. le dendrogramme complet, avec tous ses rangs, au lieu d'une seule
   partition gelée ;
2. le Gibbs conditionnel exact de chaque arbre final, au lieu d'un
   recoloriage indépendant de ses clusters ;
3. deux augmentations postérieures indépendantes, nécessaires pour
   représenter exactement une quantité quadratique de weak recovery ;
4. l'intersection des deux composantes géantes et l'opérateur d'overlap
   défini par leurs deux corridors ancestraux.

Ainsi, répéter la preuve de percolation du chapitre 11 ne suffirait pas. Le
but est de montrer qu'après élimination exacte à $`\beta_c`$, les
interactions résiduelles du Gibbs entier se contractent dans le secteur à
deux répliques. C'est cette contraction, et non la seule absence de
percolation d'un recoloriage, qui pourrait améliorer la borne actuelle.

## 2. Une réplique : Gibbs exact d'un arbre entier

Pour une observation $O$, soit $`\mu_O`$ la postérieure. Notons
$`R_O(dD\mid\sigma)`$ la loi du dendrogramme complet conditionnellement à la
configuration. La mesure augmentée est

```math
\nu_O(d\sigma,dD)
=
\mu_O(d\sigma)R_O(dD\mid\sigma).
\qquad\text{(2.1)}
```

Après intégration du spin, posons

```math
\rho_O(dD)
=
\int\nu_O(d\sigma,dD),
\qquad
\pi_{O,D}(d\sigma)
=
\nu_O(d\sigma\mid D).
\qquad\text{(2.2)}
```

Sous l'a priori produit uniforme, la conditionnelle se factorise entre les
racines finales $`R\in\mathcal R(D)`$ :

```math
\pi_{O,D}
=
\bigotimes_{R\in\mathcal R(D)}
\pi_{O,D,R}.
\qquad\text{(2.3)}
```

Pour une racine $R$, sa Gibbs exacte a la forme

```math
\pi_{O,D,R}(\sigma_R)
=
\frac1{Z_{O,D,R}}
\prod_{\substack{
u\in D\\
u\subseteq R
}}
F_{u,p}^{D}(\sigma_R).
\qquad\text{(2.4)}
```

Le produit de (2.4) porte sur **tous** les nœuds de l'arbre. Il n'est jamais
tronqué à $`\beta_c`$.

## 3. La coupe critique comme séparateur exact

Pour un bloc $`A\in\Pi_{\beta_c}(D)`$, définissons ses ports supérieurs

```math
\partial_D^+A
=
\left\{
x\in A:
x\text{ appartient à un bucket }E_u
\text{ avec }\beta_u>\beta_c
\right\}.
\qquad\text{(3.1)}
```

Pour un état de ports $s_A$, le message interne exact est

```math
Z_{D,A}(s_A)
=
\sum_{\substack{
\sigma_A\\
\sigma_{\partial_D^+A}=s_A
}}
\prod_{\substack{
u\subseteq A\\
\beta_u\le\beta_c
}}
F_{u,p}^{D}(\sigma_A).
\qquad\text{(3.2)}
```

La Gibbs supérieure d'une racine devient alors

```math
\pi_{O,D,R}(s)
\propto
\left[
\prod_{\substack{
A\in\Pi_{\beta_c}(D)\\
A\subseteq R
}}
Z_{D,A}(s_A)
\right]
\left[
\prod_{\substack{
u\subseteq R\\
\beta_u>\beta_c
}}
F_{u,p}^{D}(s)
\right].
\qquad\text{(3.3)}
```

L'identité (3.3) est une élimination par blocs, pas une approximation.

- Les états de ports d'une même racine sont tirés conjointement.
- Les rangs postcritiques réels et les buckets multiports restent présents.
- Les intérieurs des blocs ne deviennent indépendants qu'après
  conditionnement par tous leurs ports.
- Une fois les ports tirés, les intérieurs peuvent être remplis avec des
  hasards indépendants.

C'est le sens exact de « Gibbs sur les sous-arbres coupés à
$`\beta_c`$ » compatible avec un Gibbs sur tout l'arbre.

## 4. Deux répliques : identité exacte de weak recovery

Conditionnellement à la même observation $O$, tirons

```math
(\sigma^{(1)},D^{(1)}),
(\sigma^{(2)},D^{(2)})
\overset{\mathrm{i.i.d.}}{\sim}
\nu_O.
\qquad\text{(4.1)}
```

Autrement dit, $`D^{(1)}`$ et $`D^{(2)}`$ sont indépendants sous
$`\rho_O`$, puis chaque spin est tiré selon son Gibbs entier
$`\pi_{O,D^{(r)}}`$. Les tirer conditionnellement à une vérité commune ne
suffirait pas : ils resteraient corrélés par cette vérité.

Pour $`f_{ij}(\sigma)=\sigma_i\sigma_j`$,

```math
\mu_O(f_{ij})^2
=
\mathbb E\left[
\pi_{O,D^{(1)}}(f_{ij})
\pi_{O,D^{(2)}}(f_{ij})
\mid O
\right].
\qquad\text{(4.2)}
```

En posant

```math
\tau_i
=
\sigma_i^{(1)}\sigma_i^{(2)},
\qquad
\tau_i\tau_j
=
f_{ij}(\sigma^{(1)})f_{ij}(\sigma^{(2)}),
\qquad\text{(4.3)}
```

la quantité pairwise exacte est

```math
Q_L(p)
=
\frac1{n_L^2}
\sum_{i,j}
\mathbb E[
\tau_i\tau_j
].
\qquad\text{(4.4)}
```

Le carré postérieur de (4.2) emploie deux hiérarchies indépendantes.
L'enveloppe à hiérarchie commune

```math
\mathbb E_D[
\pi_{O,D}(f_{ij})^2
\mid O
]
\qquad\text{(4.5)}
```

est plus grande par Jensen. Elle reste un diagnostic licite, mais ce n'est
pas la cible exacte.

## 5. Factorisation par racines et double géante

Conditionnellement à $`(O,D^{(1)},D^{(2)})`$, la parité moyenne de
$`\tau_i\tau_j`$ est nulle si $i$ et $j$ appartiennent à deux racines
différentes dans au moins une des répliques. Les contributions se regroupent
donc selon les intersections

```math
R_{1,a}\cap R_{2,b},
\qquad
R_{r,a}\in\mathcal R(D^{(r)}).
\qquad\text{(5.1)}
```

Sur l'événement où chaque forêt possède une unique racine géante, notée
$`R_{r,\star}`$, la seule intersection susceptible d'être macroscopique est

```math
G_{12}^\star
=
R_{1,\star}\cap R_{2,\star}.
\qquad\text{(5.2)}
```

La cible principale devient

```math
\mathcal E_{L,p}^{(2),\star}
=
\frac1{n_L^2}
\mathbb E
\left[
\mathbb E_{
\pi_{O,D^{(1)},R_{1,\star}}
\otimes
\pi_{O,D^{(2)},R_{2,\star}}
}
\left[
\left(
\sum_{i\in G_{12}^\star}\tau_i
\right)^2
\right]
\right].
\qquad\text{(5.3)}
```

Il reste à borner séparément les seconds moments des racines non géantes.
Cette étape se ferme géométriquement sans supposer les deux forêts
indépendantes sous la loi annealed. En effet,

```math
\sum_{\substack{
a,b\\
(a,b)\ne(\star,\star)
}}
|R_{1,a}\cap R_{2,b}|^2
\le
\sum_{a\ne\star}|R_{1,a}|^2
+
\sum_{b\ne\star}|R_{2,b}|^2.
\qquad\text{(5.4)}
```

Pour $a\ne\star$, on somme d'abord sur les intersections qui partitionnent
$R_{1,a}$ ; pour $a=\star,b\ne\star$, chaque intersection est contenue dans
$R_{2,b}$. L'inégalité est déterministe.

Chaque forêt a marginalement la loi de la percolation de paramètre
$`q_1(p)>q_c`$. L'unicité supercritique et la sous-extensivité de la plus
grande composante hors géante donnent donc

```math
\frac1{n_L^2}
\mathbb E
\left[
\sum_{a\ne\star}|R_{1,a}|^2
+
\sum_{b\ne\star}|R_{2,b}|^2
\right]
\longrightarrow0.
\qquad\text{(5.5)}
```

La corrélation entre les deux forêts via $O$ n'intervient pas dans cette
conclusion.

## 6. Raffinement commun des deux coupes critiques

Coupons séparément $`D^{(1)}`$ et $`D^{(2)}`$ à
$`\beta_c(p)`$. Les cellules naturelles de la double géante sont

```math
C_{A_1,A_2}
=
A_1\cap A_2\cap G_{12}^\star,
\qquad
A_r\in\Pi_{\beta_c}(D^{(r)}).
\qquad\text{(6.1)}
```

Leur contribution diagonale vérifie déterministement

```math
\sum_{A_1,A_2}
|C_{A_1,A_2}|^2
\le
\min
\left\{
\sum_{A_1}|A_1|^2,
\sum_{A_2}|A_2|^2
\right\}.
\qquad\text{(6.2)}
```

La sous-extensivité quadratique des blocs critiques élimine donc la
diagonale. Le seul terme nouveau est l'énergie entre cellules critiques
distinctes de $`G_{12}^\star`$.

Plus précisément, écrivons

```math
m_r^D(i,j)
=
\pi_{O,D^{(r)}}(\sigma_i\sigma_j)
\qquad(r=1,2)
\qquad\text{(6.3)}
```

et notons $`C(i)`$ la cellule de (6.1) contenant $i$. Définissons le reste
**signé**

```math
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
=
\frac1{n_L^2}
\mathbb E
\left[
\sum_{\substack{
i,j\in G_{12}^\star\\
C(i)\ne C(j)
}}
m_1^D(i,j)m_2^D(i,j)
\right].
\qquad\text{(6.4)}
```

Les inégalités (5.4) et (6.2), avec $`|m_r^D(i,j)|\le1`$, donnent

```math
\left|
Q_L(p)
-
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
\right|
\le
\mathcal G_L^{\mathrm{hors}}(p)
+
\mathcal G_L^{\mathrm{diag}}(p),
\qquad\text{(6.5)}
```

où

```math
\begin{aligned}
\mathcal G_L^{\mathrm{hors}}(p)
&=
\frac1{n_L^2}
\mathbb E
\left[
\sum_{a\ne\star}|R_{1,a}|^2
+
\sum_{b\ne\star}|R_{2,b}|^2
\right],
\\
\mathcal G_L^{\mathrm{diag}}(p)
&=
\frac1{n_L^2}
\mathbb E
\min\left\{
\sum_{A_1}|A_1|^2,
\sum_{A_2}|A_2|^2
\right\}.
\end{aligned}
\qquad\text{(6.6)}
```

Les deux termes géométriques tendent vers zéro : le premier est
supercritique hors géante, le second est critique. Ainsi,

```math
\boxed{
Q_L(p)\longrightarrow0
\quad\Longleftrightarrow\quad
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
\longrightarrow0.
}
\qquad\text{(6.7)}
```

L'équivalence concerne la limite nulle ; elle ne remplace pas le reste signé
par sa valeur absolue. Cette réduction ferme la porte géométrique TRI0. Pour
améliorer la borne actuelle jusqu'à $`p=0.81`$, le reste
$`\mathcal E_{\mathrm{off},L}^{(2),\star}`$ est l'objet intermédiaire exact.

La
[désintégration single-$D$](41_DESINTEGRATION_PALM_RESTE_SIGNE.md)
ferme maintenant la partie scalaire de TRI1. Pour un dendrogramme $D$,
notons $`g_D(i,j)`$ l'indicatrice que les deux endpoints appartiennent à sa
plus grande racine finale, $`s_D(i,j)`$ l'indicatrice qu'ils appartiennent
au même bloc critique et $`m_D(i,j)=\pi_{O,D}(\sigma_i\sigma_j)`$. Posons

```math
d_{O,ij}
=
\mathbb E_{D\mid O}
\left[
g_D(i,j)(1-s_D(i,j))m_D(i,j)
\right].
\qquad\text{(6.8)}
```

La cible positive exacte est

```math
\mathcal D_L^\times(p)
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
d_{O,ij}^2.
\qquad\text{(6.9)}
```

Si

```math
S_L^c(p)
=
\frac1{n_L^2}
\mathbb E_{O,D}
\sum_{\substack{
A\in\Pi_{\beta_c}(D)\\
A\subseteq R_\star(D)
}}
|A|^2,
\qquad\text{(6.9a)}
```

alors la différence de carrés exacte de la note 41 donne

```math
\left|
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
-
\mathcal D_L^\times(p)
\right|
\le
2\sqrt{S_L^c(p)}.
\qquad\text{(6.10)}
```

Sous les limites géométriques critiques et supercritiques déjà utilisées
dans (5.5)--(6.6),

```math
\boxed{
Q_L(p)\longrightarrow0
\quad\Longleftrightarrow\quad
\mathcal D_L^\times(p)\longrightarrow0.
}
\qquad\text{(6.11)}
```

Le raffinement commun de deux arbres n'est généralement pas un arbre. C'est
un DAG de cellules et de ports. Chercher à le remplacer par un arbre
indépendant ferait perdre des facteurs à $D$ fixé. La désintégration (6.8)
montre cependant qu'au niveau du scalaire final, les cancellations peuvent
être prises dans la moyenne signée single-$D$ avant le carré.

## 7. Palm cross--cross et opérateur restant

La note 41 construit la Palm non circulaire adaptée à
$`\mathcal D_L^\times`$. Pour $O,i,j$ fixés, posons

```math
\delta_{O,ij}
=
\mathbb E_{D\mid O}
\left[
g_D(i,j)(1-s_D(i,j))
\right],
\qquad
\alpha_L^\times
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
\delta_{O,ij}^2.
\qquad\text{(7.1)}
```

Si $`\alpha_L^\times=0`$, alors $`\mathcal D_L^\times=0`$. Sinon, on biaise
$`(O,i,j)`$ par $`\delta_{O,ij}^2`$, puis on tire les deux dendrogrammes
indépendamment conditionnellement à ces variables selon

```math
\rho_{O,ij}^{\times}(dD)
=
\frac{
g_D(i,j)(1-s_D(i,j))
}{
\delta_{O,ij}
}
\rho_O(dD).
\qquad\text{(7.2)}
```

Les triplets avec $`\delta_{O,ij}=0`$ ont une masse Palm nulle. Sous cette
Palm cross--cross, notons $`I,J`$ les endpoints aléatoires ainsi biaisés.
Alors

```math
\mathcal D_L^\times(p)
=
\alpha_L^\times
\mathbb E_{\mathbb P_L^\times}
\left[
m_{D^{(1)}}(I,J)m_{D^{(2)}}(I,J)
\right].
\qquad\text{(7.3)}
```

La loi de Palm scalaire est donc définie. TRI2 reste ouvert parce que les
deux corridors nus ne forment pas un état récursif suffisant. Il faut
conserver leurs rangs et buckets complets, tous leurs ports, les lois de bord
extérieures exactes et une statistique suffisante de l'observation commune.
Les espaces de ports variant avec $D$, les Jacobiennes sont en général des
opérateurs entre fibres, et non des scalaires.

Si un état $x$ représente plusieurs descendants $c$, leur sélection Palm
emploie les poids positifs de paires $`\omega_c(x)`$. L'opérateur signé
correctement normalisé a la forme

```math
(\mathcal L_p^{(2)}h)(x)
=
\mathbb E_{\mathbb P_L^\times}
\left[
\sum_c
\omega_c(x)
\left(
\mathbf J_c^{(1)}\otimes\mathbf J_c^{(2)}
\right)
h(X_c)
\ \middle|\
X=x
\right],
\qquad
\omega_c(x)\ge0,
\qquad
\sum_c\omega_c(x)=1.
\qquad\text{(7.4)}
```

Les produits signés de Jacobiennes ne sont ni des probabilités, ni remplacés
par leurs valeurs absolues. Après construction d'un état suffisant et d'une
limite homogène, la première porte falsifiable sera
$`\rho(\mathcal L_{0.81}^{(2)})<1`$. Avant cette construction, il faut
contrôler les produits inhomogènes ; un rayon spectral d'un opérateur
compressé prématurément ne suffit pas.

Avant cet investissement, l'enveloppe à un dendrogramme fournit un test
go/no-go plus fort mais plus tractable :

```math
\mathbb E
\left[
\frac1{n_L}
\lambda_{\max}
\left(
W_{R_L^\star}^{1/2}
M_{R_L^\star}^c
W_{R_L^\star}^{1/2}
\right)
\right]
\longrightarrow0.
\qquad\text{(7.5)}
```

Si (7.5) tient à $`p=0.81`$, la borne de Jensen de la note 36 suffit déjà.
S'il reste macroscopique comme sur le broadcast, la route à deux
dendrogrammes est nécessaire et devra exploiter les produits signés.

## 8. Calibration exacte sur le SBM

Dans le broadcast local du SBM,

```math
d
=
\frac{a+b}{2},
\qquad
\theta
=
\frac{a-b}{a+b}.
\qquad\text{(8.1)}
```

La dérivée d'un message Gibbs à travers une arête vaut $\theta$. Deux
hiérarchies indépendantes donnent donc

```math
\mathcal L_{\mathrm{SBM}}^{(2)}
=
d\theta^2
=
\frac{(a-b)^2}{2(a+b)}.
\qquad\text{(8.2)}
```

Le seuil linéaire est exactement Kesten--Stigum sur le broadcast. Cette
égalité résulte de la marginalisation exacte et vaut pour tout ordre de
coupe ; elle ne démontre ni un sweep dynamique, ni le transfert au SBM fini.
Dans ce dernier, la balance ou les non-arêtes constituent un port global qui
recouple les racines.

La coupe physique vérifie

```math
d\,p
\left(
1-e^{-u\beta_c}
\right)
=
1,
\qquad
p=\frac{1+\theta}{2},
\qquad
u=\log\frac p{1-p}.
\qquad\text{(8.3)}
```

Elle ne produit pas le carré. Elle organise l'élimination exacte qui conduit
à (8.2). Si la même coupe est révélée aux deux répliques, le transfert devient

```math
\eta_{\mathrm{partagée}}
=
\theta^2
+
\frac{
(1/d)(1-\theta)^2
}{
1-1/d
}
>
\theta^2.
\qquad\text{(8.4)}
```

Pour $`d=3`$ et $`\theta=1/2`$,

```math
d\theta^2
=
0.75,
\qquad
d\eta_{\mathrm{partagée}}
=
1.125.
\qquad\text{(8.5)}
```

Ce contre-test exact interdit de partager le dendrogramme dans l'analyse
quadratique.

La [note pilote SBM](37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) ferme le seuil de
reconstruction du broadcast par la densité d'évolution exacte et un
sandwich global. Elle ne prouve pas encore le transfert au graphe SBM fini.

## 9. Pourquoi cette stratégie est analytique

Il serait irréaliste de chercher la loi complète de l'arbre géant de
Kruskal, puis la loi complète de son Gibbs. La réduction ci-dessus demande
seulement :

- les seconds moments des racines hors géante ;
- la sous-extensivité quadratique des blocs à la coupe critique ;
- la loi locale de deux corridors vus depuis une paire de la double géante ;
- les messages exacts sur les ports rencontrés par ces corridors ;
- le rayon spectral et une enveloppe non linéaire pour la fonction de paire.

Cette liste reste difficile, mais chaque objet est falsifiable en volume
fini. Elle évite de postuler un mélange global ou une indépendance fausse des
sous-arbres.

## 10. Almost exact et exact recovery

Sur le GSBM triangulaire de degré six avec $`p<1`$ fixé, almost exact et
exact recovery sont impossibles. Même si tous les labels voisins d'un sommet
sont révélés, son test optimal ne voit que six canaux binaires
indépendants. Son erreur est

```math
\varepsilon_6(p)
=
\sum_{k=0}^{2}
\binom6k
p^k(1-p)^{6-k}
+
10p^3(1-p)^3
>
0.
\qquad\text{(10.1)}
```

Numériquement,

```math
\varepsilon_6(0.81)
=
0.0505275094\ldots,
\qquad
\varepsilon_6(0.835805792367)
=
0.0340799611\ldots.
\qquad\text{(10.2)}
```

En posant $`\delta=1-p`$, cette erreur s'écrit exactement

```math
\varepsilon_6(1-\delta)
=
10\delta^3
-15\delta^4
+6\delta^5
\sim
10\delta^3.
\qquad\text{(10.2a)}
```

L'oracle implique une fraction d'erreur espérée uniformément positive. En
sélectionnant une famille linéaire de sommets dont les étoiles d'arêtes sont
disjointes, leurs tests oracle sont indépendants ; la probabilité de les
classer tous correctement tend vers zéro. Pour poser une question
almost/exact triangulaire non triviale, il faut donc $`p_n\to1`$, un degré
divergent ou des observations répétées.

Plus précisément, dans le régime $`p_n\to1`$, cet oracle donne les conditions
nécessaires

```math
\text{almost exact}
\quad\Longrightarrow\quad
p_n\longrightarrow1,
\qquad
\text{exact}
\quad\Longrightarrow\quad
n\varepsilon_6(p_n)\longrightarrow0.
\qquad\text{(10.2b)}
```

La seconde implication suit d'un packing linéaire d'étoiles disjointes. Par
(10.2a), elle impose $`1-p_n=o(n^{-1/3})`$. Ce sont des obstructions oracle,
pas des preuves de suffisance. Elles sont reproduites par
[`triangular_recovery_regimes_diagnostic.py`](../computations/triangular_recovery_regimes_diagnostic.py).

Le SBM divergent reste un benchmark distinct :

| objectif | observable principale | fermeture attendue |
|---|---|---|
| weak recovery | carré du transfert signé $`\mathcal D_L^\times`$ | embedding des ports + opérateur (7.4) + enveloppe non linéaire |
| almost exact | affinité de Hellinger $`\mathbb E[e^{-L/2}]`$ | erreur locale tendant vers zéro |
| exact | $`n\,\mathbb P(LX\le0)`$ | erreur locale $`o(1/n)`$ |

Dans le SBM logarithmique
$`a_n=A\log n`$, $`b_n=B\log n`$, la frontière de premier ordre à retrouver
est

```math
\left(
\sqrt A-\sqrt B
\right)^2
=
2.
\qquad\text{(10.3)}
```

Ces deux extensions ne suivent pas de la contraction quadratique de (7.4).
Le lift Hellinger correct marginalise deux fonctions de partition
hiérarchiques séparément avant de prendre leur moyenne géométrique ; la
[note 37](37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) formalise cette distinction.

## 11. Portes go/no-go

| porte | énoncé à établir | décision si échec |
|---|---|---|
| R0 — Gibbs entier | vérifier (3.3) avec tous les ports et facteurs postcritiques | refuser toute réduction qui traite les blocs critiques comme indépendants |
| R1 — deux répliques | échantillonner $`D^{(1)},D^{(2)}`$ indépendamment conditionnellement à $O$ | ne pas employer l'enveloppe à dendrogramme commun comme cible |
| R2 — réduction géométrique | (5.4)–(6.7), établie sous les faits de percolation marginaux | éliminer les racines hors géante et la diagonale critique |
| R2b — enveloppe simple | à $`L=4,p=0.81`$, $`\lambda_{\max}/n=0.9507\ldots`$ ; limite ouverte | poursuivre les volumes et les cancellations à deux $D$ |
| R3s — désintégration scalaire | (6.8)–(7.3), établie dans la note 41 | travailler sur $`\mathcal D_L^\times`$, pas sur un carré single-$D$ quenched |
| R3o — état récursif | construire l'embedding commun des ports et messages extérieurs | ne pas définir $`\mathcal L_p^{(2)}`$ sur les seuls corridors |
| R4 — linéarisation | estimer puis certifier (7.4) après construction de l'état suffisant | abandonner cette famille si la contraction signée échoue avec une marge robuste |
| R5 — non-linéaire | dominer les grands messages avec le même biais de paire | ne revendiquer aucun seuil à partir du seul rayon spectral |
| R6 — fermeture | transformer la décroissance en $`Q_L(p)\to0`$ | annoncer une borne uniquement après cette étape |

## 12. Ordre de travail immédiat

1. comparer le [port global exact](39_PORT_GLOBAL_SBM_RECOVERY.md) du SBM
   fini au broadcast ;
2. prolonger en volume le
   [diagnostic single-D](../diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md),
   défavorable à $`L=4`$, avec audit de mélange si l'énumération exacte
   devient impossible ;
3. prolonger le diagnostic signé à deux dendrogrammes au-delà de $`L=4`$ :
   l'identité (6.5) y est auditée à l'erreur machine, mais le reste moyen est
   encore positif et macroscopique ;
4. estimer $`d_{O,ij}`$ en moyennant plusieurs dendrogrammes signés pour une
   même observation, puis estimer sans biais son carré
   $`\mathcal D_L^\times`$ et la masse critique $`S_L^c`$ ;
5. construire l'embedding commun des ports et conserver les messages
   extérieurs comme variables d'état ;
6. estimer les produits signés
   $`\mathbf J_c^{(1)}\otimes\mathbf J_c^{(2)}`$ sous la Palm cross--cross,
   avec les poids de paire positifs et sans valeur absolue prématurée ;
7. ne chercher une preuve multiscalaire que si le test spectral **répliqué**
   est favorable avec une marge robuste et une enveloppe non linéaire
   compatible.

## 13. Verdict

La stratégie « coupe critique puis Gibbs » est mathématiquement cohérente et
désormais falsifiable si le mot Gibbs désigne bien le Gibbs exact de chaque
arbre entier. Pour la weak recovery, sa forme correcte est répliquée :

```math
\boxed{
D^{(1)},D^{(2)}\text{ indépendants}
\quad+\quad
\text{deux Gibbs entiers}
\quad+\quad
\text{deux coupes à }\beta_c
\quad+\quad
\text{overlap sur la double géante}.
}
\qquad\text{(13.1)}
```

Le broadcast SBM valide exactement le bookkeeping local
$`d\theta^2`$, indépendamment du niveau de coupe ; le port fini est désormais
écrit exactement dans la note 39, mais sa comparaison au broadcast reste
ouverte. Sur la grille triangulaire, les termes hors double géante et
la diagonale critique sont éliminés par (6.5), puis la note 41 désintègre le
reste signé en la cible positive
$`\mathcal D_L^\times(p)`$, à une erreur au plus
$`2\sqrt{S_L^c(p)}`$. TRI1 est donc fermé au niveau scalaire et
géométrique. La difficulté restante, TRI2, est de construire puis contracter
le transfert signé moyenné en $D$ à observation fixée, dans un espace commun
qui conserve tous les ports, messages extérieurs et facteurs postcritiques.
