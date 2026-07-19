# Corrélation spin–spin au nœud de coalescence

> [!WARNING]
> **Oracle archivé comme stratégie globale.** Les calculs locaux restent
> consultables, mais le LCA seul ne capte pas la profondeur des deux bras.
> Voir le [programme actuel](../../CURRENT_STATUS.md).

Ce fichier développe la stratégie suivante : pour une paire $i,j$, aller directement au nœud
```math
u_{ij}=\mathrm{LCA}_D(i,j)
```
où les deux sommets deviennent membres du même cluster, puis étudier la parité du heat bath entre les deux fils de ce nœud.

L'intérêt est triple :

1. la parité du flip donne exactement la conservation ou l'inversion du produit $`\sigma_i\sigma_j`$ ;
2. le heat bath est une projection conditionnelle, ce qui fournit une borne rigoureuse et un estimateur de variance réduite de la corrélation longue portée ;
3. en rafraîchissant le dendrogramme et en répétant cette mise à jour LCA, on obtient une suite de bornes qui converge vers la corrélation exacte dès que le noyau pair-spécifique est ergodique.

Le résultat principal à un pas est
```math
Q_n
\le
H_n^{\mathrm{LCA}}
:=
\frac1{n^2}\,
\mathbb E\left[
n+2\sum_{u\in D}
|C_{u,1}|\,|C_{u,2}|\,\eta_u
\right],
\qquad
\eta_u=\tanh^2\left(\frac{L_u}{2}\right).
```
La percolation de Swendsen--Wang est exactement la version non pondérée de cette formule, obtenue en remplaçant chaque $`\eta_u`$ par $1$.

## 1. Quelle corrélation ?

Pour une observation fixée $O=(X,W)$, la corrélation de Gibbs est
```math
c_{ij}(O)
=
\left\langle\sigma_i\sigma_j\right\rangle_{\mu_O}.
```
Avec l'a priori i.i.d. uniforme, l'espérance annealed brute
```math
\mathbb E[\Sigma_i\Sigma_j]
```
est nulle pour $i\ne j$. La quantité de long-range order pertinente est donc $`c_{ij}(O)`$, ou sa version mise au carré. Pour deux répliques indépendantes,
```math
\mathbb E_O[c_{ij}(O)^2]
=
\mathbb E\left\langle
\sigma_i^{(1)}\sigma_j^{(1)}
\sigma_i^{(2)}\sigma_j^{(2)}
\right\rangle.
```
L'identité de Nishimori donne aussi
```math
\mathbb E_O[c_{ij}(O)^2]
=
\mathbb E\left[
\Sigma_i\Sigma_j\,c_{ij}(O)
\right].
```
La weak recovery binaire à probabilité positive est équivalente à la non-disparition de
```math
Q_n
=
\frac1{n^2}\sum_{i,j}
\mathbb E_O[c_{ij}(O)^2].
```
## 2. Le nœud $`u_{ij}`$

Supposons d'abord $i\ne j$ et que $i$ et $j$ appartiennent au même arbre de $D$. Leur plus petit ancêtre commun possède deux fils
```math
u_{ij}:C=C_1\mathbin{\dot\cup}C_2,
\qquad
i\in C_1,\quad j\in C_2.
```
Ce nœud est le seul nœud de la hiérarchie pour lequel le produit $`\sigma_i\sigma_j`$ est exactement l'orientation relative des deux fils, à des signes internes fixes près.

- À un ancêtre strict de $`u_{ij}`$, $i$ et $j$ appartiennent au même fils : un flip de ce fils les inverse ensemble et ne change pas leur produit.
- À un descendant de $`u_{ij}`$, un flip peut changer le produit, mais il modifie une orientation interne à $`C_1`$ ou $`C_2`$, et non directement l'orientation entre les deux clusters maximaux qui séparent $i$ et $j$.

Cela rend $`u_{ij}`$ canonique. Cela ne prouve pas encore qu'il minimise l'autocorrélation parmi tous les programmes possibles ; cette question est isolée plus bas.

Si $i$ et $j$ sont dans deux racines distinctes de $`\Pi_1`$, on remplace la mise à jour au LCA par deux recolorations globales, indépendantes et uniformes, de leurs racines. Sous a priori uniforme, la fiabilité relative associée vaut alors $0$.

Cette dernière convention utilise l'a priori binaire i.i.d. uniforme. Avec un potentiel général $`\mu_0`$, les orientations de deux racines peuvent rester corrélées ; il faut alors conserver le terme inter-racines produit par le heat bath exact de $`\mu_0`$.

## 3. Les quatre événements exacts

Soient $a,b\in\{0,1\}$ les indicatrices de flip de $`C_1,C_2`$, et
```math
Z_u=q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}.
```
Si la configuration de départ est la vérité $\Sigma$, les quatre événements sont :

| Événement | $i$ après mise à jour | $j$ après mise à jour | Probabilité |
|---|---|---|---:|
| $a=0,b=0$ | conforme à $`\Sigma_i`$ | conforme à $`\Sigma_j`$ | $`q_u^{00}/Z_u`$ |
| $a=0,b=1$ | conforme | inversé | $`q_u^{01}/Z_u`$ |
| $a=1,b=0$ | inversé | conforme | $`q_u^{10}/Z_u`$ |
| $a=1,b=1$ | inversé | inversé | $`q_u^{11}/Z_u`$ |

Comme $`i\in C_1`$ et $`j\in C_2`$,
```math
\sigma_i'\sigma_j'
=
\sigma_i\sigma_j(-1)^{a+b}.
```
Par conséquent,
```math
\mathbb P\left(
\sigma_i'\sigma_j'=\sigma_i\sigma_j
\mid O,\sigma,D
\right)
=
\frac{q_u^{00}+q_u^{11}}{Z_u},
```
et
```math
\mathbb P\left(
\sigma_i'\sigma_j'=-\sigma_i\sigma_j
\mid O,\sigma,D
\right)
=
\frac{q_u^{01}+q_u^{10}}{Z_u}.
```
Définissons la persistance signée
```math
m_u(\sigma,D)
:=
\frac{
q_u^{00}+q_u^{11}
-q_u^{01}-q_u^{10}
}{Z_u}.
```
Avec
```math
L_u
=
\log
\frac{q_u^{00}+q_u^{11}}
{q_u^{01}+q_u^{10}},
```
on a
```math
m_u=\tanh\frac{L_u}{2},
\qquad
\rho_u=|m_u|,
\qquad
\eta_u=\rho_u^2.
```
Point par point,
```math
\mathbb P(\text{produit conservé}\mid O,\sigma,D)
=
\frac{1+m_u}{2}.
```
Il faut conserver le signe de $`m_u`$ dans cette formule. Le carré $`\eta_u`$ est la quantité adaptée aux bornes de corrélation.

### Où intervient exactement l'arête de Kruskal ?

L'arête de Kruskal ne remplace jamais l'ensemble des liens entre les deux fils. Elle ne fournit que le temps de la fusion
```math
\beta_u
=
\min_{\substack{e\in E_u\\e\text{ satisfaite}}}\xi_e,
\qquad
E_u=\{e:\text{une extrémité dans }C_1,\text{ l'autre dans }C_2\}.
```
Tous les liens de $`E_u`$ interviennent par
```math
\Lambda_u=\sum_{e\in E_u}|W_e|\mathbf1_{\{e\text{ satisfaite}\}},
\qquad
T_u=\sum_{e\in E_u}|W_e|.
```

Posons $`F_v(x)=x e^{(1-\beta_v)x}`$, puis séparons le facteur du nœud $u$ des messages de l'a priori et des ancêtres :

```math
A_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succ u}F_v(\Lambda_v(\sigma^{ab})),
```

```math
M_u^+=A_u^{00}+A_u^{11},
\qquad
M_u^-=A_u^{10}+A_u^{01},
\qquad
B_u=\log\frac{M_u^+}{M_u^-}.
```
Comme la parité paire laisse $`\Lambda_u`$ inchangé tandis que la parité impaire le remplace par $`T_u-\Lambda_u`$, on obtient l'identité exacte
```math
\boxed{
L_u
=
B_u
+\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
}
```
Ainsi $`\beta_u=\xi_{e_u}`$, où $`e_u`$ est une arête gagnante latente dont l'identité est marginalisée dans le dendrogramme de partitions, intervient bien dans la règle de flip ; $`\Lambda_u`$ et $`T_u`$ utilisent en revanche **tous** les liens entre $`C_1`$ et $`C_2`$. Le terme $`B_u`$ n'est nul que lorsque l'a priori et les ancêtres se compensent dans le rapport des deux parités. Le calcul exact de tous les taux ancestraux est donné dans [08_ANCESTRAL_LAMBDA_CHAIN.md](../../foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md).

## 4. Le heat bath comme projection conditionnelle

Fixons $D$ et le nœud $`u=u_{ij}`$. Notons $`\mathcal G_u`$ l'information qui reste fixe lorsque l'on parcourt l'orbite
```math
\mathcal O_u(\sigma)
=
\{\sigma^{00},\sigma^{01},\sigma^{10},\sigma^{11}\}.
```
Concrètement, $`\mathcal G_u`$ contient :

- $D$ ;
- les spins hors de $`C_1\cup C_2`$ ;
- toutes les relations de spins à l'intérieur de $`C_1`$ et de $`C_2`$ ;
- mais pas les deux orientations globales de $`C_1,C_2`$.

Pour
```math
f_{ij}(\sigma)=\sigma_i\sigma_j,
```
le heat bath au nœud est l'opérateur de projection
```math
K_uf_{ij}
=
\mathbb E_{\nu_O}\left[
f_{ij}\mid\mathcal G_u
\right].
```
Posons
```math
r_{u,ij}
:=
\mathbb E_{\nu_O}\left[
f_{ij}\mid\mathcal G_u
\right].
```
Alors
```math
K_uf_{ij}=r_{u,ij},
\qquad
m_u=f_{ij}\,r_{u,ij},
\qquad
m_u^2=r_{u,ij}^2.
```
Le signe interne reliant $i$ à l'orientation de $`C_1`$, et $j$ à celle de $`C_2`$, peut dépendre de la paire. Son carré disparaît : $`\eta_u`$ est identique pour toutes les paires
```math
(i,j)\in C_1\times C_2.
```
Pour avoir un opérateur défini pour toute réalisation de $D$, notons $`\mathsf H_{ij}`$ la projection pair-spécifique suivante : heat bath au LCA lorsque $i,j$ sont dans le même arbre, recoloration indépendante de leurs deux racines sinon. Sa sigma-algèbre conservée est notée $`\mathcal G_{ij}`$. Posons
```math
g_{ij}:=\mathsf H_{ij}f_{ij},
\qquad
m_{ij}^{\mathrm{LCA}}:=f_{ij}g_{ij},
\qquad
\eta_{ij}^{\mathrm{LCA}}:=g_{ij}^2.
```
Sur l'événement où le LCA $`u_{ij}`$ existe,
```math
m_{ij}^{\mathrm{LCA}}=m_{u_{ij}},
\qquad
\eta_{ij}^{\mathrm{LCA}}=\eta_{u_{ij}}.
```
Sur l'événement de racines distinctes, ces deux quantités valent $0$ sous a priori uniforme.

Pour la diagonale, on adopte la convention triviale
```math
\mathsf H_{ii}=\mathrm{Id},
\qquad
f_{ii}=g_{ii}=m_{ii}^{\mathrm{LCA}}=\eta_{ii}^{\mathrm{LCA}}=1.
```
## 5. Identités finies au LCA

Toutes les espérances de cette section sont conditionnelles à $O$, puis prises sous la loi jointe stationnaire $`\nu_O(\sigma,D)`$.

### Théorème fini LCA — statut : établi conditionnellement à A1

L'algèbre suivante est complète dès que la loi jointe $`\nu_O`$ de [01_MATHEMATICAL_FRAMEWORK.md](../../foundations/01_MATHEMATICAL_FRAMEWORK.md) est formalisée avec toutes ses conventions de censure. La seule dépendance restante est donc le point A1 de la feuille de route, pas une hypothèse de factorisation des nœuds.

Pour toute paire $i,j$,
```math
c_{ij}(O)
=
\mathbb E_{\nu_O}\left[
f_{ij}m_{ij}^{\mathrm{LCA}}
\right]
=
\mathbb E_{\nu_O}[g_{ij}].
```
L'autocorrélation après un heat bath au LCA vaut
```math
A_{ij}^{\mathrm{LCA}}(O)
:=
\mathbb E_{\nu_O}\left[
f_{ij}(\sigma)f_{ij}(\sigma')
\right]
=
\mathbb E_{\nu_O}[g_{ij}^2].
```
Enfin,
```math
\boxed{
A_{ij}^{\mathrm{LCA}}(O)
=
\mathbb E_{\nu_O}[m_{ij}^{\mathrm{LCA}}]
=
\mathbb E_{\nu_O}[(m_{ij}^{\mathrm{LCA}})^2]
=
\mathbb E_{\nu_O}[\eta_{ij}^{\mathrm{LCA}}].
}
```
### Preuve

La projection pair-spécifique laisse $`\nu_O`$ invariante, donc
```math
\mathbb E[f_{ij}]
=
\mathbb E[\mathsf H_{ij}f_{ij}]
=
\mathbb E[g_{ij}]
=
\mathbb E[f_{ij}m_{ij}^{\mathrm{LCA}}].
```
Comme $`\mathsf H_{ij}`$ est une projection conditionnelle dans $`L^2(\nu_O)`$,
```math
\mathbb E[f_{ij}\mathsf H_{ij}f_{ij}]
=
\mathbb E[(\mathsf H_{ij}f_{ij})^2]
=
\mathbb E[g_{ij}^2].
```
Puis $`m_{ij}^{\mathrm{LCA}}=f_{ij}g_{ij}`$, avec $`f_{ij}^2=1`$. De plus,
```math
\mathbb E[m_{ij}^{\mathrm{LCA}}\mid\mathcal G_{ij}]
=
g_{ij}\,
\mathbb E[f_{ij}\mid\mathcal G_{ij}]
=
g_{ij}^2.
```
Cela donne toutes les identités.

## 6. Probabilité moyenne de rester conforme à la vérité

Il faut distinguer deux couplages, qui répondent à deux questions différentes.

### Couplage A — on démarre de la vérité

Par Nishimori, on peut identifier la réplique stationnaire initiale à $\Sigma$, puis tirer $D\mid\Sigma$ et effectuer le heat bath au LCA.

Conditionnellement à $O,\Sigma,D$, la probabilité de conserver le produit de vérité est
```math
\frac{1+m_{ij}^{\mathrm{LCA}}(\Sigma,D)}2.
```
Après moyenne stationnaire sur $\Sigma,D$,
```math
\boxed{
\mathbb P\left(
\Sigma_i\Sigma_j=\sigma_i'\sigma_j'
\mid O
\right)
=
\frac{1+A_{ij}^{\mathrm{LCA}}(O)}2.
}
```
Les événements « les deux spins sont corrects » et « les deux spins sont inversés » donnent le même **produit de la paire**. Ils ne sont pas équivalents pour le recouvrement global si $`C_1\cup C_2`$ n'est qu'un sous-cluster : une permutation globale des labels inverse toute la configuration, pas seulement ces deux fils.

Conditionnellement à l'information quotient $`\mathcal G_{ij}`$, l'ancienne et la nouvelle orientation relative sont deux tirages indépendants de la même loi binaire. On a donc plus précisément
```math
\boxed{
\mathbb P(
\text{la relation }i-j\text{ survit}
\mid O,\mathcal G_{ij}
)
=
\frac{1+\eta_{ij}^{\mathrm{LCA}}}{2}.
}
```
Point par point avant de moyenner l'orientation initiale, la probabilité est $`(1+m_{ij}^{\mathrm{LCA}})/2`$. Après cette moyenne stationnaire, elle devient $`(1+\eta_{ij}^{\mathrm{LCA}})/2`$. Sur l'événement où le LCA existe, ce sont respectivement $`(1+m_u)/2`$ et $`(1+\eta_u)/2`$. Il ne faut donc pas remplacer $`m_u`$ par $`\rho_u`$.

Le meilleur choix déterministe de la parité, c'est-à-dire le MAP local, réussit quant à lui avec probabilité $`(1+\rho_u)/2`$. Cette règle sert à l'inférence, mais ce n'est plus le heat bath qui conserve Gibbs.

### Couplage B — une réplique postérieure indépendante de la vérité

Si $`\sigma\sim\mu_O`$ est tirée indépendamment de $\Sigma$ conditionnellement à $O$, alors
```math
\Sigma\perp(\sigma,D,\sigma')\mid O.
```
Toute suite de flips invariants conserve $`\sigma'\sim\mu_O`$, et par conséquent
```math
\mathbb P(
\Sigma_i\Sigma_j=\sigma_i'\sigma_j'
\mid O
)
=
\frac{1+c_{ij}(O)^2}{2},
```
indépendamment du nœud choisi. Le LCA n'améliore donc pas magiquement une réplique déjà stationnaire. Son intérêt est analytique : il isole exactement la parité de pont et donne une projection calculable. Pour construire un estimateur, il faut utiliser le score signé ou le MAP, et non confondre l'autocorrélation du couplage avec une nouvelle information observée.

## 7. Borne rigoureuse sur la corrélation longue portée

Par Jensen,
```math
\boxed{
c_{ij}(O)^2
\le
A_{ij}^{\mathrm{LCA}}(O).
}
```
Sur une réalisation $D$ où les racines sont distinctes, la recoloration indépendante donne une contribution conditionnelle nulle à $`A_{ij}^{\mathrm{LCA}}`$.

Cette borne est exacte au sens où son écart est
```math
A_{ij}^{\mathrm{LCA}}(O)-c_{ij}(O)^2
=
\mathrm{Var}_{\nu_O}(g_{ij}).
```
La qualité du critère LCA se ramène donc à une question précise : le biais conditionnel $`g_{ij}`$ se concentre-t-il, pour les paires lointaines, autour de sa moyenne ?

## 8. Score LCA global : une percolation pondérée par l'information

Chaque paire non ordonnée $i\ne j$ dans un même arbre possède un unique LCA. Toutes les paires traversant les deux fils de $u$ ont la même valeur $`\eta_u`$. Par conséquent,
```math
\boxed{
Q_n
\le
H_n^{\mathrm{LCA}}
:=
\frac1{n^2}
\mathbb E\left[
n
+2\sum_{u\in D}
|C_{u,1}|\,|C_{u,2}|\,\eta_u
\right].
}
```
L'espérance porte sur $O,\sigma,D$. Les paires appartenant à des racines distinctes contribuent zéro.

Ainsi,
```math
H_n^{\mathrm{LCA}}\longrightarrow0
\quad\Longrightarrow\quad
Q_n\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
Le calcul de la somme ne demande que $O(|D|)$ termes, et non $O(n^2)$ paires.

Cette borne raffine exactement celle de Swendsen--Wang. Pour chaque arbre binaire,
```math
|C|^2
=
|C_1|^2+|C_2|^2+2|C_1||C_2|.
```
En télescopant sur la forêt,
```math
n+2\sum_{u\in D}|C_{u,1}||C_{u,2}|
=
\sum_{R\text{ racine de }D}|R|^2.
```
Comme $`0\le\eta_u\le1`$,
```math
\boxed{
H_n^{\mathrm{LCA}}
\le
\frac1{n^2}\,
\mathbb E\sum_{R\text{ racine}}|R|^2.
}
```
La borne de percolation ne retient que l'existence des fusions. La borne LCA affecte à chaque fusion le facteur de transmission $`\eta_u`$. Elle peut donc conclure même si des composantes sont grandes, lorsque les fusions responsables d'une fraction positive des paires sont peu fiables.

Le poids $`2|C_{u,1}||C_{u,2}|`$ est exactement le nombre de paires ordonnées dont le LCA est $u$. Une contribution macroscopique peut provenir de quelques séparations équilibrées, mais aussi de nombreux nœuds déséquilibrés dont les poids cumulés sont d'ordre $n^2$.

### Mesure des fusions informatives

En conservant le temps de coalescence, définissons
```math
\mathcal M_n(dt)
=
\frac2{n^2}\,
\mathbb E\left[
\sum_{u\in D}
|C_{u,1}||C_{u,2}|\eta_u\,
\delta_{\beta_u}(dt)
\right].
```
Alors
```math
H_n^{\mathrm{LCA}}
=
\frac1n+\int_0^1\mathcal M_n(dt).
```

La même mesure sans $`\eta_u`$ décrit seulement la géométrie single-linkage/percolation. Le facteur $`\eta_u`$ ajoute la quantité d'information transmise au moment précis où les deux clusters fusionnent.

Sous l'a priori binaire i.i.d. uniforme, pour une coupe déterministe $`\beta`$, notons $`\Pi_\beta`$ les composantes à ce niveau et

```math
S_n(\beta)
=
\frac1{n^2}\mathbb E\sum_{C\in\Pi_\beta}|C|^2.
```

La partie du score née avant la coupe est dominée par $`S_n(\beta)`$. On obtient donc

```math
\boxed{
Q_n
\le
S_n(\beta)
+
\mathcal M_n((\beta,1]).
}
```

Si $`S_n(\beta)\to0`$, toute weak recovery impose ainsi une masse informative strictement positive de fusions au-dessus de la coupe. Avec un potentiel $`\mu_0`$ général, il faut ajouter la contribution inter-racines. Le [critère de bande critique](../../diagnostics/07_CRITICAL_BAND_CRITERION.md) sépare ensuite la connexion quotient, la fiabilité locale et la cohérence signée nécessaires pour rendre cette intuition exacte.

### Version strictement longue portée

Pour une distance géométrique $R$, soit
```math
\mathcal P_R
=
\{\{i,j\}:d(i,j)\ge R\}.
```
Notons
```math
N_u(R)
=
\#\{
(i,j)\in C_{u,1}\times C_{u,2}:
d(i,j)\ge R
\}.
```
Alors
```math
\frac1{|\mathcal P_R|}
\sum_{\{i,j\}\in\mathcal P_R}
\mathbb E[c_{ij}(O)^2]
\le
\frac1{|\mathcal P_R|}
\mathbb E\left[
\sum_{u\in D}N_u(R)\eta_u
\right].
```
Cette formulation sépare la taille des clusters de la persistance informationnelle de leur fusion.

## 9. Répéter le LCA : une suite vers le critère exact

Le résultat à un pas se prolonge naturellement. Pour une paire fixée $i\ne j$, définissons le noyau marginal $`K_{ij}^{\mathrm{LCA}}`$ :

1. tirer un nouveau $`D\sim\nu_O(\cdot\mid\sigma)`$ ;
2. si $i,j$ ont un LCA, effectuer le heat bath exactement à ce nœud ;
3. sinon, sous a priori uniforme, recolorer indépendamment les deux racines ;
4. oublier $D$.

Sur l'espace augmenté, l'étape 2 est une projection conditionnelle. Si $`A:L^2(\mu_O)\to L^2(\nu_O)`$ désigne l'inclusion isométrique $`(Af)(\sigma,D)=f(\sigma)`$, alors, après marginalisation de $D$,
```math
K_{ij}^{\mathrm{LCA}}=A^*\mathsf H_{ij}A
```
est un noyau de Markov auto-adjoint, positif et contractant dans $`L^2(\mu_O)`$.

Pour $`f_{ij}=\sigma_i\sigma_j`$, posons
```math
A_{ij}^{(m)}(O)
=
\left\langle
f_{ij},
(K_{ij}^{\mathrm{LCA}})^m f_{ij}
\right\rangle_{\mu_O}.
```
Alors
```math
1=A_{ij}^{(0)}
\ge A_{ij}^{(1)}
\ge A_{ij}^{(2)}
\ge\cdots
\ge c_{ij}(O)^2,
```
et le premier terme non trivial est précisément
```math
A_{ij}^{(1)}(O)
=
\mathbb E_{\nu_O}\left[
\eta_{ij}^{\mathrm{LCA}}
\right].
```
Si ce noyau est ergodique pour l'observable $`f_{ij}`$,
```math
\boxed{
A_{ij}^{(m)}(O)\downarrow c_{ij}(O)^2.
}
```
Sur un graphe fini, avec tous les $`|W_e|<\infty`$ et un a priori uniforme de support plein, l'événement « aucune horloge n'ouvre avant $1$ » a une probabilité strictement positive. Les sommets sont alors des racines séparées et la parité de $i,j$ est rééchantillonnée. C'est une piste de contraction, mais pas encore une preuve d'ergodicité observable : il faut identifier l'espace fixe de $`K_{ij}^{\mathrm{LCA}}`$ et exclure les invariants extérieurs encore corrélés à $`f_{ij}`$. De plus, la constante obtenue peut décroître catastrophiquement avec $n$ et ne fournit aucun contrôle uniforme du mélange. Les interactions dures $`|W_e|=\infty`$ demandent un traitement séparé.

En effet, pour la mesure spectrale $`\mu_{ij}^{\mathrm{sp}}`$ de $`f_{ij}-c_{ij}`$, supportée dans $`[0,1]`$,
```math
A_{ij}^{(m)}(O)-c_{ij}(O)^2
=
\int_{[0,1]}\lambda^m\,d\mu_{ij}^{\mathrm{sp}}(\lambda).
```
Dans le couplage démarrant à la vérité,
```math
\mathbb P(
\sigma_i^{(m)}\sigma_j^{(m)}
=
\Sigma_i\Sigma_j
\mid O
)
=
\frac{1+A_{ij}^{(m)}(O)}2.
```
La limite est $`(1+c_{ij}(O)^2)/2`$, c'est-à-dire l'accord de deux répliques postérieures indépendantes. Cette suite répond exactement à la question « que devient la relation vraie après plusieurs flips LCA ? ».

En sommant sur les paires,
```math
Q_n
\le
H_n^{(m)}
:=
\frac1{n^2}\sum_{i,j}\mathbb E_O[A_{ij}^{(m)}(O)],
\qquad
H_n^{(m)}\downarrow Q_n
```
sous l'hypothèse d'ergodicité et les conditions permettant le passage à la limite. Seul $m=1$ se réduit immédiatement à la somme $O(|D|)$ sur les nœuds ; calculer efficacement les termes $m\ge2$ est un nouvel objectif.

Cette construction est pair-spécifique : le LCA choisi dépend de $(i,j)$. Elle fournit simultanément une famille de bornes scalaires, mais pas un unique parcours de flips réalisant tous les $`K_{ij}`$ à la fois.

## 10. Estimateur de corrélation par Rao–Blackwell

Pour un tirage stationnaire $(\sigma,D)$, définissons
```math
\widehat c_{ij}^{\mathrm{LCA}}
:=
\sigma_i\sigma_j\,m_{ij}^{\mathrm{LCA}},
```
Cette convention donne automatiquement la valeur $0$ pour deux racines distinctes et $1$ sur la diagonale.

Alors
```math
\mathbb E_{\nu_O}
\left[
\widehat c_{ij}^{\mathrm{LCA}}
\right]
=
c_{ij}(O).
```
Il s'agit du conditionnement de la variable brute $`\sigma_i\sigma_j`$ par l'orbite du heat bath. Sa variance vaut
```math
\mathrm{Var}
\left(
\widehat c_{ij}^{\mathrm{LCA}}
\right)
=
A_{ij}^{\mathrm{LCA}}(O)-c_{ij}(O)^2
\le
1-c_{ij}(O)^2.
```
Le LCA fournit donc simultanément :

- une borne supérieure sur $`c_{ij}^2`$ via $`\eta_u`$ ;
- un estimateur non biaisé de $`c_{ij}`$ via le score signé $`\sigma_i\sigma_jm_u`$.

Une valeur élevée de $`\eta_u`$ seule n'est pas une preuve de récupération : elle peut refléter une information conditionnelle qui change de signe. La suffisance doit utiliser le score signé.

## 11. Structure matricielle hiérarchique

Pour un tirage $(\sigma,D)$, la matrice
```math
\widehat C^{\mathrm{LCA}}_{ij}
=
\widehat c_{ij}^{\mathrm{LCA}}
```
est constante, à un facteur $`\sigma_i\sigma_j`$ près, sur chaque bloc $`C_{u,1}\times C_{u,2}`$.

Pour un vecteur $x$,
```math
(\widehat C^{\mathrm{LCA}}x)_i
=
x_i
+\sum_{u\in\mathrm{path}(i)}
m_u\sigma_i
\sum_{j\in\mathrm{sibling}_u(i)}
\sigma_jx_j.
```
Les sommes de sous-arbres se calculent de bas en haut, puis les contributions des frères se propagent de haut en bas. Une multiplication matrice–vecteur coûte ainsi $O(n+|D|)$, sans hypothèse d'équilibrage et sans matérialiser $n^2$ entrées.

La matrice obtenue pour un tirage unique est symétrique, mais n'est pas nécessairement positive semi-définie : chaque entrée utilise la projection propre à son LCA. Son espérance est bien la matrice de corrélation $`C_O`$, qui est positive semi-définie. Une mise en œuvre spectrale doit donc moyenner les opérateurs signés ou utiliser une procédure symétrisée, sans attribuer à chaque échantillon une positivité non démontrée.

Programme algorithmique :

1. échantillonner plusieurs couples $(\sigma,D)$ ;
2. calculer $`m_u`$ à tous les nœuds ;
3. moyenner implicitement les opérateurs $\widehat C^{\mathrm{LCA}}$ ;
4. extraire leur vecteur propre principal ;
5. comparer son signe à une réplique ou à la vérité dans les simulations.

## 12. Le LCA est-il optimal ?

Le LCA est **canonique**, mais son optimalité parmi tous les nœuds n'est pas automatique.

Pour tout heat bath admissible $K$, posons
```math
A_{ij}(K)
=
\left\langle f_{ij},Kf_{ij}\right\rangle_{\nu_O}.
```
Pour un heat bath qui est une projection conditionnelle,
```math
A_{ij}(K)
=
\left\|
\mathbb E[f_{ij}\mid\mathcal G_K]
\right\|_2^2.
```
Si deux sigma-algèbres sont emboîtées,
```math
\mathcal G_1\subseteq\mathcal G_2,
```
alors
```math
A_{ij}(K_1)\le A_{ij}(K_2).
```
Mais les sigma-algèbres d'un heat bath au LCA et d'un heat bath à un descendant ne sont généralement pas emboîtées. Il faut donc tester, et non supposer,
```math
K_{u_{ij}}
\in
\arg\min_{K\in\mathcal K_{ij}}
A_{ij}(K).
```
### Caractère canonique dans une classe restreinte

Le LCA est l'unique nœud qui :

- sépare $i$ et $j$ dans ses deux fils ;
- conserve toutes les relations internes à ces deux fils ;
- met à jour directement leur orientation relative.

Il est donc le candidat canonique dans cette classe restreinte de mouvements « orientation relative de deux clusters maximaux ». Parler d'optimalité demande encore de fixer un critère — autocorrélation, variance, coût ou précision MAP — puis de comparer les noyaux admissibles.

### Extension plus forte

Une mise à jour jointe de toutes les variables d'orientation sur le chemin entre $i$ et $j$, ou dans le sous-arbre enraciné en $`u_{ij}`$, conditionne sur moins d'information et peut fournir une borne plus petite. Elle est mathématiquement plus difficile, mais constitue l'étape suivante si le LCA seul reste trop peu précis.

Le [fichier 16](16_FLIP_PROBABILITIES_DESCENDANT_PATH.md) formalise cette
extension : il donne l'identité déterministe de parité sur les deux bras,
calcule l'oracle factorisé PATH-FAC et montre pourquoi la véritable
marginalisation descendante doit transporter un message indexé par la
frontière, plutôt qu'un seul scalaire par niveau.

Le [fichier 17](17_PATH_DECORRELATION_THRESHOLD.md) identifie ensuite le
critère de perte : l'atténuation factorisée doit diverger, ou, pour le
balayage joint, les normes $L^2$ des opérateurs tordus doivent avoir un
produit tendant vers zéro.

Plus généralement, si un « LCA collapsed » rééchantillonne la parité de $u$ en marginalisant aussi une partie des messages ancestraux, sa sigma-algèbre $`\mathcal G_u^{\mathrm{coarse}}`$ vérifie
```math
\mathcal G_u^{\mathrm{coarse}}\subseteq\mathcal G_u.
```
Jensen donne alors
```math
\left\|
\mathbb E[f_{ij}\mid\mathcal G_u^{\mathrm{coarse}}]
\right\|_2^2
\le
\left\|
\mathbb E[f_{ij}\mid\mathcal G_u]
\right\|_2^2.
```
Le programme d'optimalité devient donc : trouver la plus petite information hiérarchique que l'on peut marginaliser exactement à coût raisonnable. À l'extrême, un nouveau tirage postérieur indépendant donne déjà $`c_{ij}^2`$, mais cette « solution » est tautologique ; le LCA à $D$ fixé est l'extrémité locale calculable.

## 13. Formule locale dans le GSBM triangulaire

Dans le modèle homogène, tous les modules valent
```math
u_p=\log\frac p{1-p}.
```
Pour une coupe déterministe comportant $m$ liens, dont $k$ sont satisfaits, et sous l'hypothèse que l'a priori et les ancêtres sont neutres dans le rapport de parité,
```math
L_{m,k,\beta}^{\mathrm{loc}}
=
\log\frac{k}{m-k}
+(1-\beta)u_p(2k-m).
```
Ainsi,
```math
\mathbb P(\text{produit conservé}\mid m,k,\beta)
=
\frac1{1+e^{-L_{m,k,\beta}^{\mathrm{loc}}}},
```
et
```math
\eta_{m,k,\beta}^{\mathrm{loc}}
=
\tanh^2
\left(
\frac{L_{m,k,\beta}^{\mathrm{loc}}}{2}
\right).
```
Sous la loi annealed, pour une coupe déterministe indépendante des observations,
```math
k\sim\mathrm{Bin}(m,p).
```
Conditionnellement à $k\ge1$, la sous-densité d'une fusion à $`\beta\in[0,1]`$ est
```math
ku_pe^{-ku_p\beta}\,d\beta.
```
Une fiabilité locale moyenne, conditionnelle à une fusion avant $1$, est donc
```math
\bar\eta_m(p)
=
\frac{
\displaystyle
\sum_{k=1}^m
\binom mkp^k(1-p)^{m-k}
\int_0^1
\eta_{m,k,\beta}^{\mathrm{loc}}\,
ku_pe^{-ku_p\beta}\,d\beta
}{
\displaystyle
\sum_{k=1}^m
\binom mkp^k(1-p)^{m-k}
\left(1-e^{-ku_p}\right)
}.
```
Conventions :

- $k=0$ : aucune fusion issue de liens satisfaits ;
- $k=m$ : $L=+\infty$, donc $\eta=1$ ;
- $2k=m$ : le terme énergétique s'annule, et la fiabilité dépend seulement du préfacteur.

Cette formule est exacte pour une coupe déterministe sous les hypothèses locales. Pour un LCA sélectionné par Kruskal, la loi de $`k`$ conditionnellement au squelette et à $`(m,\beta)`$ est donnée ci-dessous ; le choix du nœud biaise la géométrie $`(E_u,m,\beta)`$, et les facteurs ancêtres doivent être réintroduits.

### Première correction exacte du biais au temps de fusion

Dans le couplage annealed homogène, notons
```math
q_p(t)=p(1-e^{-u_pt})
```
la probabilité qu'une arête soit ouverte à l'instant $t$. Conditionnellement au fait qu'un lien soit encore fermé juste avant $t$, la probabilité qu'il soit néanmoins satisfait, avec une horloge résiduelle supérieure à $t$, vaut
```math
s_p(t)
=
\frac{pe^{-u_pt}}{1-p+pe^{-u_pt}}
=
\mathrm{logistic}(u_p(1-t)).
```

La partie biaisée du vote est exactement la probabilité conditionnelle d'une arête de la bande :

```math
h_p(t)
=2s_p(t)-1
=
\tanh\left(\frac{u_p(1-t)}2\right).
```

Sur une coupe de taille $m$, l'existence d'une arête tardive est gouvernée par $`m h_p(t)`$, tandis que la stabilité statistique du vote l'est par $`m h_p(t)^2`$. Cette différence explique pourquoi une connexion de bande peut apparaître sans transmettre encore une information longue portée robuste.

Si une arête ouvrante fusionne à l'instant $t$ deux composantes séparées par une coupe fixée de $m$ liens, cette arête est satisfaite et les $m-1$ autres liens de la coupe sont fermés. Conditionnellement à la filtration juste avant $t$,

```math
\boxed{
k
\ \stackrel{d}{=}\
1+\mathrm{Bin}(m-1,s_p(t)).
}
```
Cette loi remplace $\mathrm{Bin}(m,p)$ au moment de l'ouverture. Elle est exacte dans ce cadre annealed conditionnellement au squelette non marqué, même lorsque la coupe $`E_u`$ est choisie par Kruskal. La géométrie aléatoire de la coupe reste biaisée et le message $`B_u`$ des ancêtres reste à intégrer ; leur désintégration exacte est donnée dans [08_ANCESTRAL_LAMBDA_CHAIN.md](../../foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md).

Elle donne le diagnostic suivant pour la contribution locale $`B_u=0`$ :

- si $m=1$, la fusion verrouille la parité et $`\eta_u=1`$ ;
- si $m=2$, le cas $k=1$ randomise parfaitement la parité, tandis que $k=2$ la verrouille ;
- pour $t<1$ fixé et $m\to\infty$, la majorité conditionnelle devient forte et $`\eta_u\to1`$ en probabilité ;
- à $`t=1`$, $`k=1+\mathrm{Bin}(m-1,1/2)`$ et, dans le modèle local sans
  message d'ancêtre, $`\eta_u=((2k-m)/m)^2`$ et
  $`\mathbb E\eta_u=1/m`$ exactement.

Cela suggère qu'un seul flip LCA peut rester proche de la borne percolative lorsque les coalescences longue portée surviennent tôt et à travers de grandes coupes. La suite itérée $`H_n^{(m)}`$, ou la propagation du message complet $`B_u`$, devient alors essentielle.

### Temps de coalescence de deux sommets

Pour deux sommets de la grille,
```math
\beta_{ij}
=
\inf\{t:i\leftrightarrow j\text{ dans }\Pi_t\}
=
\min_{\gamma:i\leadsto j}\max_{e\in\gamma}\xi_e.
```
Si $`\tau_{ij}(q)=\mathbb P_q(i\leftrightarrow j)`$ désigne la fonction de connexion de la percolation par arêtes,
```math
\mathbb P(\beta_{ij}\le t)=\tau_{ij}(q_p(t)).
```
Le problème triangulaire se décompose donc proprement en deux objets : la mesure géométrique des fusions, accessible par la percolation proche du seuil, et la loi conditionnelle de $`\eta_u`$, qui contient l'information de spin et le message des ancêtres.

## 14. Conditions nécessaire et suffisante : cible précise

La borne LCA donne immédiatement une condition nécessaire :
```math
H_n^{\mathrm{LCA}}\to0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
La réciproque est fausse sans hypothèse supplémentaire. Une cible possible est :

### Conjecture de netteté LCA

Dans une classe de modèles homogènes à préciser,
```math
\frac1{n^2}\sum_{i,j}
\mathbb E\left[
A_{ij}^{\mathrm{LCA}}(O)-c_{ij}(O)^2
\right]
\longrightarrow0.
```
Sous cette hypothèse,
```math
H_n^{\mathrm{LCA}}-Q_n\longrightarrow0,
```
et la non-disparition $`\liminf_n H_n^{\mathrm{LCA}}>0`$ devient équivalente à la weak recovery au sens fixé dans le dossier.

Cette conjecture demande que **le premier pas** soit déjà asymptotiquement net. La chaîne pair-spécifique fournit une voie plus robuste : pour chaque volume fini et chaque paire ergodique,
```math
H_n^{(m)}\downarrow Q_n
\qquad(m\to\infty).
```
Le prochain théorème nécessaire-et-suffisant à viser est donc une version uniforme permettant de choisir $`m=m_n`$ telle que
```math
H_n^{(m_n)}-Q_n\longrightarrow0.
```
Le problème se déplace alors vers un contrôle de mélange spectral du noyau LCA pair-spécifique. Cette convergence caractérise l'information ; elle ne fournit pas encore un algorithme polynomial unique, puisque les noyaux diffèrent selon la paire.

Une autre voie de suffisance consiste à montrer que la matrice signée
```math
\mathbb E[
\widehat C^{\mathrm{LCA}}\mid O
]
=
C_O
```
peut être estimée avec assez de précision pour que son vecteur propre principal fournisse un recouvrement non trivial.

## 15. Premiers lemmes à démontrer

La feuille de route détaillée est centralisée dans [05_PROOF_ROADMAP.md](../roadmaps/05_PROOF_ROADMAP.md). Pour la voie LCA, les quatre clôtures immédiates sont :

1. finaliser A1, puis rédiger le théorème pair-spécifique avec les cas LCA, racines distinctes et diagonale ;
2. formaliser la représentation spectrale du noyau itéré et séparer l'ergodicité à volume fini de son contrôle uniforme en $n$ ;
3. vérifier les identités sur une arête, un chemin, un triangle et un cactus ;
4. comparer quantitativement $`H_n^{\mathrm{LCA}}`$ à l'information-percolation avant toute annonce de nouveau seuil.
