# Pont de volume fini : port global du SBM et trois régimes de recovery

**Statut : identités finies exactes pour le posterior et sa décomposition par
racines ; seuils classiques retrouvés comme benchmarks informationnels ;
aucune preuve de mélange hiérarchique ni nouvelle preuve arbre--graphe.**

Cette note complète le
[pilote broadcast](37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md). Elle répond à la
question suivante :

> après avoir coupé chaque arbre du dendrogramme et tiré son Gibbs, que faut-il
> encore ajouter pour représenter exactement le SBM fini ?

La réponse est un **port global scalaire**. Dans le modèle à labels i.i.d.,
ce port porte le potentiel des non-arêtes ; dans la bisection plantée, il
porte la contrainte de balance. Il interdit de tirer naïvement les arbres
finaux indépendamment, mais il admet une élimination exacte par convolution
des magnétisations de racines.

Le verdict scientifique est double.

1. Si l'on autorise le tirage exact de ce Gibbs global, la hiérarchie redonne
   automatiquement la postérieure et donc tous ses seuils statistiques. Ce
   constat est une identité d'augmentation de données, pas une explication
   nouvelle de ces seuils.
2. Pour obtenir une preuve hiérarchique, il reste à contrôler ce port, le
   nombre de sweeps et les queues du log-likelihood. La coupe critique fournit
   un ordre d'élimination exact ; elle ne sélectionne à elle seule ni
   Kesten--Stigum, ni l'almost exact, ni l'exact recovery.

## 1. Posterior fini exact

Considérons le SBM symétrique à deux classes avec labels
$`x_i\in\{-1,+1\}`$, $`0<b_n<a_n<n`$ et

```math
\mathbb P(A_{ij}=1\mid x_ix_j=+1)
=
\frac{a_n}{n},
\qquad
\mathbb P(A_{ij}=1\mid x_ix_j=-1)
=
\frac{b_n}{n},
\qquad
a_n>b_n>0.
\qquad\text{(1.1)}
```

Posons

```math
h_{1,n}
=
\frac12\log\frac{a_n}{b_n},
\qquad
h_{0,n}
=
\frac12
\log
\frac{1-a_n/n}{1-b_n/n}
<0,
\qquad
J_n=h_{1,n}-h_{0,n}>0.
\qquad\text{(1.2)}
```

Sous le prior i.i.d. uniforme, le posterior exact est

```math
\mu_A(x)
\propto
\exp\left[
h_{0,n}\sum_{i<j}x_ix_j
+
J_n\sum_{\{i,j\}\in E(A)}x_ix_j
\right].
\qquad\text{(1.3)}
```

Comme

```math
\sum_{i<j}x_ix_j
=
\frac12
\left[
\left(\sum_i x_i\right)^2-n
\right],
\qquad\text{(1.4)}
```

les non-arêtes se réduisent à un potentiel de la magnétisation totale
$`M(x)=\sum_i x_i`$ :

```math
\mu_A(x)
\propto
\exp\left[
\frac{h_{0,n}}2M(x)^2
\right]
\exp\left[
J_n\sum_{\{i,j\}\in E(A)}x_ix_j
\right].
\qquad\text{(1.5)}
```

Le facteur constant $`\exp(-nh_{0,n}/2)`$ a été omis.

Dans la bisection plantée exactement équilibrée, où $n$ est pair, le premier
facteur est constant, mais le prior impose

```math
\mathbf 1_{\{M(x)=0\}}.
\qquad\text{(1.6)}
```

Il n'existe donc pas de formulation usuelle qui donne simultanément un prior
produit, la suppression exacte des non-arêtes et l'indépendance des racines.

Le taux ferromagnétique exact du facteur porté par une arête présente est

```math
u_n^{\mathrm{fin}}
=
2J_n
=
\log
\frac{
a_n(1-b_n/n)
}{
b_n(1-a_n/n)
}.
\qquad\text{(1.7)}
```

Il diffère du taux du canal de broadcast
$`u_n^{\mathrm{br}}=\log(a_n/b_n)`$. Dans le régime
$`a_n,b_n=O(\log n)`$,

```math
u_n^{\mathrm{fin}}-u_n^{\mathrm{br}}
=
O\left(\frac{a_n-b_n}{n}\right).
\qquad\text{(1.8)}
```

Cette proximité locale ne supprime pas le port global de (1.5).

## 2. Compression exacte par magnétisation de racine

Augmentons chaque facteur d'arête présente par les liens et horloges du
dendrogramme, puis fixons un dendrogramme complet $D$. La représentation
Edwards--Sokal absorbe chaque facteur d'arête dans une racine finale de la
forêt. Notons $`\mathcal R(D)`$ l'ensemble de ces racines et
$`w_{D,R}(x_R)`$ leur poids interne exact, facteurs postcritiques compris.

Pour chaque racine, définissons le message de magnétisation

```math
W_{D,R}(m)
=
\sum_{\substack{
x_R\in\{-1,+1\}^{R}\\
\sum_{i\in R}x_i=m
}}
w_{D,R}(x_R).
\qquad\text{(2.1)}
```

Alors la fonction de partition conditionnelle du SBM i.i.d. vaut exactement

```math
Z_D^{\mathrm{iid}}
=
\sum_{(m_R)}
\left[
\prod_{R\in\mathcal R(D)}
W_{D,R}(m_R)
\right]
\exp\left[
\frac{h_{0,n}}2
\left(
\sum_R m_R
\right)^2
\right].
\qquad\text{(2.2)}
```

Pour la bisection plantée, elle vaut

```math
Z_D^{\mathrm{bal}}
=
\sum_{(m_R)}
\left[
\prod_{R\in\mathcal R(D)}
W_{D,R}(m_R)
\right]
\mathbf1_{\{\sum_Rm_R=0\}}.
\qquad\text{(2.3)}
```

Les équations (2.2)--(2.3) sont la factorisation correcte. Les racines ne sont
pas indépendantes, mais toute leur dépendance extérieure est transportée par
le seul entier $`\sum_Rm_R`$.

Dans le SBM ferromagnétique conditionné par le $D$ **complet**, chaque lien
ouvert impose l'égalité de ses deux spins. Chaque racine finale est donc
monochromatique et (2.1) se réduit à

```math
W_{D,R}(m)
=
c_{D,R}
\left(
\mathbf1_{\{m=|R|\}}
+
\mathbf1_{\{m=-|R|\}}
\right).
\qquad\text{(2.3a)}
```

Le port est alors un weighted subset-sum sur les tailles $`|R|`$, et non un
nouveau Gibbs interne difficile. La forme générale (2.1) reste utile pour
une hiérarchie partiellement marginalisée ou pour le GSBM signé, mais elle
n'est pas nécessaire au benchmark SBM full-$D$.

### 2.1 Élimination exacte du port

Numérotons les racines $`R_1,\ldots,R_r`$ et posons

```math
F_k(s)
=
\sum_{\substack{
m_1+\cdots+m_k=s
}}
\prod_{\ell=1}^k
W_{D,R_\ell}(m_\ell).
\qquad\text{(2.4)}
```

La récursion

```math
F_{k+1}(s)
=
\sum_m
F_k(s-m)W_{D,R_{k+1}}(m)
\qquad\text{(2.5)}
```

permet de contracter exactement le port. Le poids terminal est
$`\exp(h_{0,n}s^2/2)`$ dans le modèle i.i.d. et
$`\mathbf1_{\{s=0\}}`$ dans la bisection. Un passage arrière donne ensuite
un tirage exact des magnétisations de racines.

Dans le SBM full-$D$, les messages (2.3a) sont explicites et (2.5) est une
convolution unidimensionnelle, polynomiale en $n$ avec l'implémentation
naïve. Une difficulté exponentielle de calcul des messages internes peut
réapparaître dans le GSBM signé ou après marginalisation partielle de $D$ ;
elle ne doit pas être attribuée à ce benchmark ferromagnétique.

Le module
[`sbm_global_port_convolution.py`](../computations/sbm_global_port_convolution.py)
implémente exactement cette convolution conditionnellement aux tailles de
racines fournies :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_global_port_convolution.py \
  --root-sizes 3 2 1 --a 4 --b 1
```

Dans cet audit, les huit orientations donnent
$`N(0)=2`$, $`Z_D^{\mathrm{iid}}=2.851724288`$ après retrait des facteurs
internes communs, et une erreur nulle entre convolution et énumération
directe. Le script ne tire pas $D$ : il certifie l'étape algébrique (2.5),
pas la loi des tailles de racines ni un seuil de recovery.

### 2.2 Pourquoi le port n'est pas une petite correction à Kesten--Stigum

Sur le broadcast, la forêt auxiliaire finale a pour nombre moyen de
descendants ouverts

```math
d\theta.
\qquad\text{(2.6)}
```

À Kesten--Stigum, $`d\theta^2=1`$, donc

```math
d\theta
=
\frac1\theta
>1
\qquad
\text{pour }0<\theta<1.
\qquad\text{(2.7)}
```

La coupe critique est ainsi strictement incluse dans une forêt finale
supercritique. Des racines macroscopiques sont précisément le régime où une
contrainte sur $`\sum_Rm_R`$ peut produire une dépendance d'ordre un entre
orientations. On ne peut donc pas supprimer le port en invoquant seulement
la convergence locale vers un arbre.

### 2.3 Deux répliques : ne pas créer de port croisé

Pour la weak recovery, les deux postérieures sont indépendantes
conditionnellement à $A$. Chacune possède son propre port :

```math
\mu_A(x^{(1)})\mu_A(x^{(2)})
\propto
\exp\left[
\frac{h_{0,n}}2
\left(
M(x^{(1)})^2+M(x^{(2)})^2
\right)
\right]
\prod_{r=1}^2
\exp\left[
J_n\sum_{\{i,j\}\in E(A)}
x_i^{(r)}x_j^{(r)}
\right].
\qquad\text{(2.8)}
```

Il n'y a aucun terme $`M(x^{(1)})M(x^{(2)})`$. Dans la bisection, chaque
copie satisfait naturellement sa propre contrainte $`M(x^{(r)})=0`$. Ce qui
changerait l'expérience serait une convolution commune, un facteur croisé,
la contrainte $`M(x^{(1)})=M(x^{(2)})`$ ou une coupe partagée ; cette
dernière gonfle déjà le Jacobien d'une arête.

## 3. Coupe critique : deux taux à ne pas confondre

Sur le broadcast $`\mathrm{PGW}(d_n)`$, posons

```math
d_n=\frac{a_n+b_n}{2},
\qquad
p_n=\frac{a_n}{a_n+b_n},
\qquad
\theta_n=\frac{a_n-b_n}{a_n+b_n}.
\qquad\text{(3.1)}
```

La coupe géométrique du pilote vérifie

```math
d_np_n
\left(
1-e^{-u_n^{\mathrm{br}}\beta_{c,n}^{\mathrm{br}}}
\right)
=1.
\qquad\text{(3.2)}
```

Comme $`d_np_n=a_n/2`$,

```math
\beta_{c,n}^{\mathrm{br}}
=
\frac{
-\log(1-2/a_n)
}{
\log(a_n/b_n)
},
\qquad
a_n>2.
\qquad\text{(3.3)}
```

Cette coupe appartient à l'horizon $`[0,1]`$ exactement lorsque

```math
d_n\theta_n
=
\frac{a_n-b_n}{2}
\ge1.
\qquad\text{(3.4)}
```

Sur le SBM fini, remplacer simplement $`u_n^{\mathrm{br}}`$ par
$`u_n^{\mathrm{fin}}`$ ne rend pas (3.2) exacte : le port global corrèle les
satisfactions des arêtes. La quantité (3.3) est la coupe du modèle local de
broadcast et la limite asymptotique naturelle, pas un seuil de percolation
i.i.d. exact du posterior fini.

Lorsque $`a_n,b_n\to\infty`$ avec un rapport borné loin de zéro et un,
$`\beta_{c,n}^{\mathrm{br}}`$ est d'ordre $`1/a_n`$. Les blocs restent
critiques en nombre moyen de descendants, mais cette finesse géométrique ne
produit pas à elle seule la concentration du log-likelihood.

## 4. Weak recovery : ce qui est retrouvé et ce qui ne l'est pas

Le paramètre classique est

```math
\lambda_n
=
d_n\theta_n^2
=
\frac{(a_n-b_n)^2}{2(a_n+b_n)}.
\qquad\text{(4.1)}
```

Après marginalisation séparée de deux coupes, une arête du broadcast
transmet $`\theta_n^2`$ dans le secteur overlap. La ramification donne donc
$`\lambda_n`$. Pour des paramètres constants, la densité d'évolution sur
l'arbre retrouve exactement la frontière $`\lambda=1`$, conformément au
théorème de
[Mossel--Neeman--Sly](https://arxiv.org/abs/1311.4115).

Ce résultat a trois portées distinctes.

| niveau | conclusion correcte |
|---|---|
| identité statistique | le Gibbs augmenté, marginalisé exactement, est la postérieure |
| calibration locale | le Jacobien overlap du broadcast vaut $`d\theta^2`$ |
| théorème dynamique fini | encore ouvert pour un nombre explicite de sweeps et avec le port global |

Supposer disponible un heat bath global exact rendrait le premier niveau
tautologique : un tirage de la postérieure possède par définition ses seuils
optimaux. Pour expliquer ou prouver le seuil à l'aide de la dynamique, il faut
plutôt contrôler, pour une observable $`f_{ij}(x)=x_ix_j`$,

```math
\left\|
K_A^t f_{ij}
\right\|_{L^2(\mu_A)}^2
\qquad\text{(4.2)}
```

avec un $`t`$ et une complexité explicites, puis comparer cette quantité au
broadcast sans supprimer le port (2.2) ou (2.3).

Le premier lemme de transfert à viser est donc une comparaison locale
uniforme du type

```math
\frac1{n^2}
\sum_{i,j}
\mathbb E
\left[
\left\|
K_A^t f_{ij}
\right\|_{L^2(\mu_A)}^2
\right]
\le
\Phi_t(\lambda_n)+o(1),
\qquad\text{(4.3)}
```

où $`\Phi_t(\lambda)\to0`$ sous $`\lambda\le1`$. Écrire
$`\Phi_t`$ sans utiliser déjà le théorème de non-reconstruction est la porte
SBM-F.

La fermeture visée doit préciser l'ordre des limites :

```math
\lim_{t\to\infty}
\limsup_{n\to\infty}
\frac1{n^2}
\sum_{i,j}
\mathbb E
\left[
\left\|
K_A^t f_{ij}
\right\|_{L^2(\mu_A)}^2
\right]
=0.
\qquad\text{(4.4)}
```

## 5. Almost exact recovery : le bon exposant

Pour tester un label dont tous les autres labels pertinents sont révélés,
les deux hypothèses échangent des lois de Bernoulli de paramètres
$`a_n/n`$ et $`b_n/n`$. Si $`m_n`$ labels sont révélés dans chaque groupe,
leur affinité de Bhattacharyya exacte est

```math
H_{v,n}
=
\left[
\frac{\sqrt{a_nb_n}}n
+
\sqrt{
\left(1-\frac{a_n}{n}\right)
\left(1-\frac{b_n}{n}\right)
}
\right]^{2m_n}.
\qquad\text{(5.1)}
```

Pour $`2m_n=n+O(1)`$ et $`a_n,b_n=O(\log n)`$,

```math
-\log H_{v,n}
=
C_n+o(1),
\qquad
C_n
=
\frac12
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2.
\qquad\text{(5.2)}
```

L'erreur de Bayes oracle vérifie

```math
\frac{H_{v,n}^2}{4}
\le
P_{e,n}^{\mathrm{oracle}}
\le
\frac{H_{v,n}}2.
\qquad\text{(5.3)}
```

Ainsi, dans ce régime, l'erreur locale oracle tend vers zéro si et seulement
si $`C_n\to\infty`$. Ce critère est équivalent à
$`\lambda_n\to\infty`$, car

```math
C_n
\le
\lambda_n
\le
2C_n.
\qquad\text{(5.4)}
```

La hiérarchie retrouve donc le bon benchmark seulement si son élimination
séparée sous $`x_v=+1`$ et $`x_v=-1`$ reproduit (5.1). Une dynamique qui ne
transporte que le second moment d'overlap ne peut pas certifier cette queue.

Dans la bisection équilibrée, révéler tous les autres labels détermine
$`x_v`$ par le compte. Le benchmark oracle doit cacher une paire opposée ou
un sous-ensemble résiduel de taille divergente. Pour une paire opposée
cachée, l'affinité de l'expérience de permutation vaut, avec le coefficient
entre crochets de (5.1) noté $`\mathfrak b_n`$,

```math
H_{\mathrm{pair},n}
=
\mathfrak b_n^{\,2(n-2)},
\qquad\text{(5.5)}
```

et non l'affinité mono-sommet (5.1). Le doublement de l'exposant est compensé
par le nombre quadratique de paires candidates à l'échelle de l'exact
recovery ; la frontière finale reste la même, mais les deux expériences ne
doivent pas être identifiées.

Pour obtenir une achievability almost exact globale sans circularité, il faut
encore :

1. partir d'un initialiseur weak avec un avantage uniforme sur le hasard ;
2. amplifier cet avantage jusqu'à une fraction $`o(1)`$ d'erreurs, par
   bootstrap avec des données séparées ou un argument leave-one-out ;
3. borner uniformément la contamination du log-likelihood par les labels
   encore erronés.

Si l'on suppose dès le départ un initialiseur à $`o(1)`$ erreurs, la
hiérarchie n'analyse que l'étape de raffinement et ne prouve pas l'almost
exact recovery elle-même.

Aucune de ces trois étapes ne découle de l'identité de Gibbs. Le benchmark
est cohérent avec les résultats de degré divergent
d'[Abbe--Sandon](https://arxiv.org/abs/1506.03729).

## 6. Exact recovery : même géométrie, queue plus forte

Dans le régime

```math
a_n=A\log n,
\qquad
b_n=B\log n,
\qquad
A>B>0,
\qquad\text{(6.1)}
```

l'exposant local vaut

```math
C_n
=
\frac12
\left(
\sqrt A-\sqrt B
\right)^2
\log n.
\qquad\text{(6.2)}
```

Cette identité est exacte pour le proxy
$`C_n=\frac12(\sqrt{a_n}-\sqrt{b_n})^2`$. C'est l'exposant binomial fini
$`-\log H_{v,n}`$ qui vaut $`C_n+o(1)`$ sous les hypothèses de (5.2).

La frontière classique d'
[Abbe--Bandeira--Hall](https://arxiv.org/abs/1405.3267) est

```math
\boxed{
\left(
\sqrt A-\sqrt B
\right)^2
=2.
}
\qquad\text{(6.3)}
```

Au-dessus de la frontière stricte, une preuve hiérarchique devrait produire
une erreur locale $`o(1/n)`$ après initialisation et raffinement. En dessous,
l'affinité seule ne suffit pas à prouver l'impossibilité : il faut
l'asymptotique précise

```math
P_{e,n}^{\mathrm{oracle}}
=
n^{-C_n/\log n+o(1)}
\qquad\text{(6.4)}
```

et un second moment du nombre de sommets ambigus. Le cas d'égalité demande
des corrections de second ordre.

La coupe (3.3) est encore un ordre d'élimination exact, mais elle ne change
pas la frontière (6.3). Le lift correct n'est plus le produit de deux
corrélations : il marginalise séparément les fonctions de partition sous les
deux hypothèses, puis prend leur moyenne géométrique.

## 7. Trois fonctionnels et trois obligations de preuve

| objectif | fonctionnel local | obstruction restante |
|---|---|---|
| weak | second moment du message postérieur | dynamique finie et port global |
| almost exact | affinité $`\mathbb E[e^{-L/2}]`$ | concentration et contamination |
| exact | queue $`\mathbb P(L\le0)`$ à l'échelle $`1/n`$ | grande déviation précise et dépendances |

Le même dendrogramme peut servir d'ordre d'élimination aux trois lignes, mais
il ne permet pas de réutiliser la même inégalité probabiliste.

### 7.1 Ordre de travail falsifiable

1. Coupler l'implémentation vérifiée de (2.4)--(2.5) à des dendrogrammes
   SBM échantillonnés, en conservant le port dans chaque observable.
2. Comparer une coupe critique, une coupe légèrement sous-critique et une
   élimination sans coupe à coût interne fixé.
3. Mesurer $`\|K_A^tf\|_2^2`$ pour des sweeps explicites, pas seulement deux
   tirages déjà stationnaires.
4. Répéter l'élimination sous les deux hypothèses locales et contrôler
   directement l'affinité (5.1).
5. N'annoncer une preuve arbre--graphe qu'après une borne uniforme sur le
   port et les cycles.

### 7.2 Critères go/no-go

- Si le résultat dépend du partage d'une coupe ou d'un port entre deux
  répliques, il ne calcule pas le carré postérieur.
- Si $`\beta_c`$ disparaît après marginalisation et qu'aucune borne de coût
  ou de contraction ne l'utilise, la coupe est une organisation du calcul,
  pas le mécanisme du seuil.
- Si un diagnostic d'almost exact donne une erreur tendant vers zéro à degré
  borné, il contredit le test oracle.
- Si l'exact recovery est conclu par une simple union bound sur (5.3), la
  constante de seuil n'est pas certifiée.

## 8. Verdict

Le SBM classique confirme l'architecture à deux niveaux :

```math
\boxed{
\text{Gibbs exact dans chaque racine}
\quad+\quad
\text{convolution du port global}.
}
\qquad\text{(8.1)}
```

Cette architecture est exacte en volume fini et conserve la stratégie de
coupe critique. Elle explique aussi pourquoi « Gibbs indépendants pour chaque
arbre » doit être limité aux modèles à prior produit sans facteur global,
comme le GSBM triangulaire dans sa formulation actuelle, ou au broadcast
edge-only.

Elle ne fournit pas encore une nouvelle preuve des seuils. Le prochain gain
mathématique propre au SBM serait une comparaison quantitative du port
global avec la limite de broadcast pour (4.2), suivie d'un contrôle de
Hellinger pour (5.1). Pour le GSBM triangulaire, le port global n'existe pas ;
le verrou prioritaire reste le reste signé inter-cellules de la
[double géante](38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md).
