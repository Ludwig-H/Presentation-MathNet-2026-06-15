# Deux dendrogrammes coupés à $\beta_c^{\mathrm{geom}}$

## 1. Limite locale : le canal de broadcast

Le voisinage d'un sommet uniforme du SBM clairsemé converge localement vers
un arbre de Galton--Watson $\mathrm{PGW}(d)$. Le label d'un enfant $t$,
conditionnellement au label de son parent $s$, suit

```math
P_\theta(t\mid s)
=
\frac{1+\theta st}{2},
\qquad s,t\in\{-1,+1\},
\qquad
\theta=\frac{a-b}{a+b}.
\qquad\text{(1.1)}
```

Écrivons

```math
p=\mathbb P(t=s)=\frac{1+\theta}{2},
\qquad
u=\log\frac p{1-p}=\log\frac ab.
\qquad\text{(1.2)}
```

Cette limite locale est le banc d'essai exact de la construction. Le
[SBM fini](04_DYNAMIQUE_HIERARCHIQUE.md#6-le-port-global-du-sbm-fini) garde
en plus un potentiel de magnétisation ou une contrainte de balance.

## 2. Une famille d'horloges, donc un vrai dendrogramme

Conditionnellement à une configuration $\sigma$, donnons à chaque arête
satisfaite une horloge

```math
\xi_e\sim\mathrm{Exp}(u),
\qquad
\xi_e=+\infty
\quad\text{si l'arête n'est pas satisfaite}.
\qquad\text{(2.1)}
```

Pour $0\le\beta\le1$, la partition $\Pi_\beta$ est formée des composantes
connexes des arêtes telles que $\xi_e\le\beta$. La même famille
$(\xi_e)_e$ est utilisée pour tous les niveaux : les partitions sont donc
emboîtées. Il ne faut pas redessiner une percolation indépendante pour
chaque $\beta$.

Le bit de coupe d'une arête est

```math
B_e(\beta)=\mathbf1_{\{\xi_e\le\beta\}}.
\qquad\text{(2.2)}
```

Sa probabilité annealed vaut

```math
q_\beta
=
\mathbb P(B_e(\beta)=1)
=
p(1-e^{-u\beta}),
\qquad
q_0=0,
\qquad
q_1=\theta.
\qquad\text{(2.3)}
```

La partition $\Pi_1$ est exactement la partition Swendsen--Wang du canal
ferromagnétique local. Si la forêt ne devient pas connexe, le dendrogramme
est simplement une forêt hiérarchique tronquée à $1$.

Deux auxiliaires seront utilisés et ne doivent pas être confondus :

- $B_\beta=(B_e(\beta))_e$ est la **projection à une coupe**. Lorsqu'on
  conditionne seulement par $B_\beta$, les horloges futures sont
  marginalisées et une arête fermée garde le canal résiduel calculé en
  section 5 ;
- $D$ est le **dendrogramme complet**, temps de fusion compris. À $D$ fixé,
  une coupe sert seulement à ordonner l'élimination de la conditionnelle
  $\nu_A(\cdot\mid D)$.

Les identités scalaires des sections 4–7 portent sur la projection
$B_\beta$, ou, ce qui revient au même, sur $D$ après marginalisation de
toute l'information autre que cette coupe. La dynamique full-$D$ est
décrite séparément dans la [note 04](04_DYNAMIQUE_HIERARCHIQUE.md).

## 3. La coupe géométrique $\beta_c^{\mathrm{geom}}$

La coupe géométrique est définie par

```math
q_{\beta_c^{\mathrm{geom}}}=\frac1d.
\qquad\text{(3.1)}
```

Lorsque $dp>1$,

```math
\boxed{
\beta_c^{\mathrm{geom}}
=
-\frac1u
\log\left(1-\frac1{dp}\right).
}
\qquad\text{(3.2)}
```

Elle appartient à l'horizon physique $[0,1]$ si et seulement si

```math
d\theta\ge1.
\qquad\text{(3.3)}
```

Au seuil de Kesten--Stigum, $d\theta^2=1$, donc
$d\theta=1/\theta>1$ dès que $0<\theta<1$. La coupe critique existe ainsi
bien avant la fin du dendrogramme dans le régime où la partition finale
Swendsen--Wang possède déjà une géante.

Cette séparation des deux seuils est fondamentale :

| phénomène | nombre de branchement |
|---|---:|
| blocs géométriques à la coupe | $q_{\beta_c^{\mathrm{geom}}}=1/d$ |
| géante Swendsen--Wang finale | $d\theta=1$ |
| reconstruction de l'information | $d\theta^2=1$ |

### Le temps informationnel $\beta_\chi$

Le secteur à deux répliques produit la rétention $\theta^2$. Pour la placer
sur la même échelle que la coupe géométrique, définissons

```math
q_{\beta_\chi}=\theta^2,
\qquad
\beta_\chi
=
-\frac1u
\log\left(
1-\frac{\theta^2}{p}
\right).
\qquad\text{(3.4)}
```

Comme $0<\theta^2\le q_1=\theta$, ce temps appartient toujours à
$[0,1]$. Dès que $\beta_c^{\mathrm{geom}}$ est défini, la croissance
stricte de $q_\beta$ donne

```math
\boxed{
\begin{aligned}
\beta_\chi<\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2<1,\\
\beta_\chi=\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2=1,\\
\beta_\chi>\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2>1.
\end{aligned}
}
\qquad\text{(3.5)}
```

L'égalité centrale est la représentation exacte du seuil demandée : le temps
de la percolation d'information rencontre la coupe géométrique précisément à
Kesten--Stigum. Il s'agit d'un couplage par égalité de lois marginales. Le
temps $\beta_\chi$ n'est ni une température, ni une coupe auxiliaire révélée
à l'estimateur, ni l'intersection littérale des deux forêts.

## 4. Couper sans perdre un seul facteur

Pour un bit $b\in\{0,1\}$, posons

```math
\Psi_{\theta,\beta}(s,t,b)
=
P_\theta(t\mid s)
\mathbb P(B_e(\beta)=b\mid s,t).
\qquad\text{(4.1)}
```

Par construction,

```math
\sum_{b=0}^1
\Psi_{\theta,\beta}(s,t,b)
=
P_\theta(t\mid s).
\qquad\text{(4.2)}
```

Sur tout arbre fini $T$, l'insertion des bits de coupe donne donc

```math
\prod_{e\in T}P_\theta(\sigma_{\mathrm{enf}(e)}
\mid\sigma_{\mathrm{par}(e)})
=
\sum_{B\in\{0,1\}^{E(T)}}
\prod_{e\in T}
\Psi_{\theta,\beta}(
\sigma_{\mathrm{par}(e)},
\sigma_{\mathrm{enf}(e)},
B_e).
\qquad\text{(4.3)}
```

On peut sommer d'abord les spins à l'intérieur des composantes de
$\Pi_{\beta_c^{\mathrm{geom}}}$, puis leurs états de bord, puis les bits de
coupe. C'est une
simple associativité de sommes finies : le Gibbs obtenu est encore le Gibbs
de l'arbre entier.

En particulier :

- les facteurs correspondant aux arêtes fermées à la coupe restent
  présents ;
- les composantes coupées ne sont pas déclarées indépendantes lorsqu'elles
  partagent des facteurs supérieurs ;
- la coupe est un ordre d'élimination, pas une approximation du modèle.

## 5. Canal résiduel de la coupe

Posons $q=q_\beta$. Conditionnellement à $B_e=1$, l'arête impose $s=t$ et
sa corrélation vaut $c_1=1$. Conditionnellement à $B_e=0$, sa corrélation
résiduelle vaut

```math
c_0
=
\frac{\theta-q}{1-q}.
\qquad\text{(5.1)}
```

Avec $\pi_1=q$ et $\pi_0=1-q$,

```math
\sum_{b=0}^1\pi_bc_b=\theta.
\qquad\text{(5.2)}
```

L'équation (5.2) dit que marginaliser la coupe restitue exactement le canal
initial.

À la coupe critique, $q=1/d$ et

```math
c_0^{\,c}
=
\frac{\theta-1/d}{1-1/d}.
\qquad\text{(5.3)}
```

## 6. La mauvaise construction : une coupe partagée

La construction partagée correspond à un autre couplage. On tire d'abord
$D$ sous sa marginale

```math
\rho_A(dD)
=
\int\mu_A(d\sigma)R_A(dD\mid\sigma),
\qquad\text{(6.1)}
```

puis $\sigma^{(1)},\sigma^{(2)}$ indépendamment sous
$\pi_{A,D}=\nu_A(\,\cdot\mid D)$. Chaque spin a encore pour marginale
$\mu_A$, mais les deux spins ne sont plus indépendants conditionnellement à
$A$ : ce couplage calcule
$\mathbb E_{D\mid A}[\pi_{A,D}(f)^2]$, et non
$\mu_A(f)^2$.

Sur une arête du broadcast, les deux facteurs voient alors le même bit
$B$. Le transfert d'overlap devient

```math
\eta_{\mathrm{partagée}}
=
\sum_b\pi_bc_b^2
=
\theta^2
+
\frac{q(1-\theta)^2}{1-q}
>
\theta^2
\qquad(q>0,\ \theta<1).
\qquad\text{(6.2)}
```

Il s'agit de la moyenne d'un carré. Le bit commun est une information
auxiliaire révélée aux deux copies ; l'expérience est plus informative que
le SBM original.

À $q=1/d$,

```math
\eta_{\mathrm{partagée}}
=
\theta^2
+
\frac{(1/d)(1-\theta)^2}{1-1/d}.
\qquad\text{(6.3)}
```

Par exemple, si $d=3$ et $\theta=1/2$,

```math
d\theta^2=0.75,
\qquad
d\eta_{\mathrm{partagée}}=1.125.
\qquad\text{(6.4)}
```

La coupe partagée déclare donc artificiellement le modèle supercritique
alors qu'il est sous Kesten--Stigum.

## 7. La bonne construction : deux dendrogrammes

Conditionnellement à la même observation $A$, tirons

```math
(\sigma^{(1)},D^{(1)}),
(\sigma^{(2)},D^{(2)})
\overset{\mathrm{i.i.d.}}{\sim}
\nu_A.
\qquad\text{(7.1)}
```

Les deux copies ont :

- le même niveau déterministe $\beta_c^{\mathrm{geom}}$ ;
- deux spins postérieurs indépendants ;
- deux familles d'horloges indépendantes conditionnellement à leurs spins ;
- deux partitions $\Pi_{\beta_c^{\mathrm{geom}}}^{(1)}$ et
  $\Pi_{\beta_c^{\mathrm{geom}}}^{(2)}$ qui peuvent être différentes.

Le transfert d'une arête est alors

```math
\begin{aligned}
\eta_{\mathrm{indépendante}}
&=
\sum_{b_1,b_2}
\pi_{b_1}\pi_{b_2}c_{b_1}c_{b_2}
\\
&=
\left(
\sum_b\pi_bc_b
\right)^2
=
\theta^2.
\end{aligned}
\qquad\text{(7.2)}
```

Il s'agit du produit de deux moyennes. Les identités (6.2) et (7.2) sont
des identités **annealed du canal local et de la projection de coupe**.
Elles ne sont pas des identités point par point pour le SBM fini
conditionné par son port global.

Une fois le produit $\theta^2$ obtenu, l'identité
$q_{\beta_\chi}=\theta^2$ le replace sur l'horloge géométrique. C'est à ce
moment seulement que (3.5) donne
$\beta_\chi=\beta_c^{\mathrm{geom}}$ au seuil. Construire à la place
l'intersection brute des arêtes ouvertes dans les deux forêts donnerait
$q_\beta^2$ et représenterait un autre objet.

Le bon ordre des opérations est donc

```math
\boxed{
\left(\sum_{D^{(1)}}W_1(D^{(1)})\right)
\left(\sum_{D^{(2)}}W_2(D^{(2)})\right),
}
\qquad\text{(7.3)}
```

et non

```math
\boxed{
\sum_D W_1(D)W_2(D).
}
\qquad\text{(7.4)}
```

La première expression calcule le carré postérieur ; la seconde calcule
une expérience augmentée à dendrogramme commun.

## 8. Raffinement commun sans hasard commun

Pour comparer les deux sorties, on peut former après coup les cellules

```math
C_{A_1,A_2}
=
A_1\cap A_2,
\qquad
A_r\in\Pi_{\beta_c^{\mathrm{geom}}}^{(r)}.
\qquad\text{(8.1)}
```

Ce raffinement commun est une opération géométrique déterministe sur deux
partitions déjà tirées. Il ne couple ni les horloges ni les spins. C'est la
bonne manière de dessiner une « double géante » ou des blocs communs sans
retomber dans le dendrogramme partagé.

## 9. La phrase à retenir

```math
\boxed{
\text{même observation et même niveau }\beta_c^{\mathrm{geom}},
\quad
\text{mais deux spins, deux horloges, deux dendrogrammes
et deux marginalisations.}
}
```

La [note suivante](03_PREUVE_DU_SEUIL_WEAK_RECOVERY.md) montre comment le
facteur $\theta^2$ devient exactement le seuil $d\theta^2=1$.
