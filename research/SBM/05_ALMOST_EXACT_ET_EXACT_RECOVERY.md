# Almost exact et exact recovery

## 1. Le changement de régime

Avec $a,b$ constants, une proportion positive de sommets ne reçoit pas
assez d'information locale. Almost exact et exact recovery sont donc
impossibles. Il faut laisser varier

```math
\mathbb P(A_{ij}=1\mid X_iX_j=+1)=\frac{a_n}{n},
\qquad
\mathbb P(A_{ij}=1\mid X_iX_j=-1)=\frac{b_n}{n}.
\qquad\text{(1.1)}
```

Les trois objectifs ne propagent pas le même fonctionnel :

| objectif | quantité naturelle | échelle critique |
|---|---|---|
| weak recovery | carré d'overlap | $\lambda_n=d_n\theta_n^2$ |
| almost exact | erreur locale ou affinité de Hellinger | divergence de Chernoff $\to\infty$ |
| exact | queue uniforme de l'erreur locale | erreur $o(1/n)$ |

La hiérarchie géométrique peut rester la même, mais son **lift
probabiliste** doit changer.

## 2. Expérience locale oracle

Pour classifier un sommet $v$, imaginons que les labels des autres sommets
soient connus sans révéler $X_v$ par une contrainte globale. Cette
expérience est littérale sous l'a priori i.i.d. ; pour une bisection
exactement équilibrée, il faut cacher au moins une paire opposée ou utiliser
un leave-block-out.

Posons

```math
p_n=\frac{a_n}{n},
\qquad
q_n=\frac{b_n}{n},
\qquad
\rho_n
=
\sqrt{p_nq_n}
+
\sqrt{(1-p_n)(1-q_n)}.
\qquad\text{(2.1)}
```

$\rho_n$ est le coefficient de Bhattacharyya d'une observation de
Bernoulli. Sous l'a priori i.i.d., l'expérience oracle observe exactement
les $n-1$ arêtes incidentes, indépendantes conditionnellement aux labels.
Son affinité locale vaut donc exactement

```math
H_{v,n}
=
\rho_n^{\,n-1}.
\qquad\text{(2.2)}
```

Définissons

```math
I_n=-2\log\rho_n.
\qquad\text{(2.3)}
```

Alors

```math
-\log H_{v,n}
=
\frac{(n-1)I_n}{2}.
\qquad\text{(2.4)}
```

Dans le régime clairsemé $a_n,b_n=o(n)$,

```math
\frac{nI_n}{2}
=
\left(1+o(1)\right)
\frac{
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
}{2}.
\qquad\text{(2.5)}
```

Si, plus précisément, $a_n,b_n=O(\log n)$ — ou si
$(a_n+b_n)^2/n=o(1)$ — le reste peut être pris additif :

```math
\frac{nI_n}{2}
=
\frac{
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
}{2}
+o(1).
\qquad\text{(2.6)}
```

Notons

```math
C_n
=
\frac{
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
}{2}.
\qquad\text{(2.7)}
```

Si $P_{e,n}^{\mathrm{oracle}}$ est l'erreur de Bayes locale, les relations
entre variation totale et affinité donnent

```math
\frac{H_{v,n}^2}{4}
\le
P_{e,n}^{\mathrm{oracle}}
\le
\frac{H_{v,n}}2.
\qquad\text{(2.8)}
```

## 3. Almost exact recovery

Le benchmark scalaire est

```math
\boxed{
nI_n\longrightarrow\infty
}
\qquad\text{ou, sous (2.5),}\qquad
\boxed{
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
\longrightarrow\infty.
}
\qquad\text{(3.1)}
```

Dans le SBM binaire symétrique, cette condition est équivalente à

```math
\lambda_n
=
\frac{(a_n-b_n)^2}{2(a_n+b_n)}
\longrightarrow\infty,
\qquad\text{(3.2)}
```

car

```math
\lambda_n
\le
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
\le
2\lambda_n.
\qquad\text{(3.3)}
```

Le double broadcast devient donc de plus en plus supercritique. Cette
observation explique qualitativement pourquoi l'overlap tend vers un, mais
elle ne donne pas à elle seule le taux optimal d'erreur.

Une preuve hiérarchique d'achievability doit encore :

1. produire une initialisation ayant $o(n)$ erreurs ;
2. séparer l'information utilisée pour l'initialisation de celle utilisée
   pour le raffinement ;
3. contrôler l'effet des labels encore faux sur les messages de blocs ;
4. montrer que l'erreur locale suit l'exposant (2.4).

## 4. Exact recovery

Dans le régime logarithmique

```math
a_n=A\log n,
\qquad
b_n=B\log n,
\qquad
A>B>0,
\qquad\text{(4.1)}
```

on a

```math
C_n
=
\frac{
\left(
\sqrt A-\sqrt B
\right)^2
}{2}
\log n
+o(\log n).
\qquad\text{(4.2)}
```

Le seuil de premier ordre est

```math
\boxed{
\left(
\sqrt A-\sqrt B
\right)^2=2.
}
\qquad\text{(4.3)}
```

Au-dessus, l'erreur locale optimale est assez petite pour contrôler
simultanément les $n$ sommets après une initialisation presque exacte. En
dessous, un nombre non négligeable de sommets reste localement ambigu.

Les inégalités strictes donnent la transition sans ambiguïté :

```math
\left(
\sqrt A-\sqrt B
\right)^2>2
\quad\Longrightarrow\quad
\text{exact recovery possible},
\qquad\text{(4.4)}
```

```math
\left(
\sqrt A-\sqrt B
\right)^2<2
\quad\Longrightarrow\quad
\text{exact recovery impossible}.
\qquad\text{(4.5)}
```

À l'égalité, les préfacteurs et les hypothèses de positivité de $A,B$
comptent ; ce dossier ne remplace pas cette analyse de fenêtre critique par
un argument de premier ordre.

Le carré global d'overlap ne suffit pas : $Q_n\to1$ tolère encore
$o(n)$ erreurs, et même un nombre borné d'erreurs détruit l'exact recovery.

## 5. Le lift Hellinger correct

Pour les deux hypothèses $X_v=+1$ et $X_v=-1$, soit
$W_\pm(o,D)$ le poids non normalisé d'une observation locale $o$ et d'un
dendrogramme auxiliaire. Les vraisemblances marginales sont

```math
L_\pm(o)
=
\int W_\pm(o,D_\pm)\,dD_\pm.
\qquad\text{(5.1)}
```

L'affinité pertinente est

```math
H_v
=
\frac1{\sqrt{Z_+Z_-}}
\sum_o
\sqrt{
\left[
\int W_+(o,D_+)\,dD_+
\right]
\left[
\int W_-(o,D_-)\,dD_-
\right]
}.
\qquad\text{(5.2)}
```

Chaque hiérarchie est donc marginalisée **avant** la racine carrée. La
quantité

```math
\frac1{\sqrt{Z_+Z_-}}
\int
\sum_o
\sqrt{W_+(o,D)W_-(o,D)}
\,dD
\qquad\text{(5.3)}
```

partage un dendrogramme entre les deux hypothèses. Elle décrit une
expérience augmentée plus informative et n'est pas l'affinité du SBM
observé.

Il faut aussi distinguer (5.2) des deux répliques ordinaires de weak
recovery : ici les deux côtés représentent deux **hypothèses locales**, pas
deux échantillons de la même postérieure.

## 6. Pourquoi la coupe critique tend naturellement vers zéro

Pour le canal de broadcast,

```math
d_n=\frac{a_n+b_n}{2},
\qquad
p_n^{=}
=
\frac{a_n}{a_n+b_n},
\qquad
u_n=\log\frac{a_n}{b_n}.
\qquad\text{(6.1)}
```

Comme $d_np_n^{=}=a_n/2$, la coupe géométrique vérifie

```math
\boxed{
\beta_{c,n}^{\mathrm{geom}}
=
-\frac{
\log(1-2/a_n)
}{
\log(a_n/b_n)
}.
}
\qquad\text{(6.2)}
```

Une solution finie de (6.2) existe seulement si $a_n>2$. Elle tombe dans
l'horizon du dendrogramme $[0,1]$ si et seulement si

```math
d_n\theta_n
=
\frac{a_n-b_n}{2}
\ge1,
\qquad\text{c'est-à-dire}\qquad
a_n-b_n\ge2.
\qquad\text{(6.3)}
```

La conclusion $\beta_{c,n}^{\mathrm{geom}}\to0$ ne se limite pas à un
rapport signal/bruit fixé. Dans le régime standard d'almost exact recovery,
la section 3 donne $\lambda_n=d_n\theta_n^2\to\infty$. Or
$u_n=2\,\mathrm{atanh}(\theta_n)\ge2\theta_n$ et
$p_n^{=}\ge1/2$, donc $p_n^{=}u_n\ge\theta_n$. Ainsi
$d_np_n^{=}u_n\ge d_n\theta_n\ge\lambda_n$, et
$d_np_n^{=}\to\infty$. Par conséquent,

```math
\beta_{c,n}^{\mathrm{geom}}
=
\frac{1+o(1)}{d_np_n^{=}u_n}
\le
\frac{1+o(1)}{\lambda_n}
\longrightarrow0.
\qquad\text{(6.4)}
```

Le cas proportionnel donne un taux particulièrement lisible. Supposons que

```math
a_n=A s_n,
\qquad
b_n=B s_n,
\qquad
A>B>0,
\qquad
s_n\longrightarrow\infty,
\qquad\text{(6.5)}
```

on obtient

```math
\beta_{c,n}^{\mathrm{geom}}
\sim
\frac{
2
}{
A\log(A/B)\,s_n
}
\longrightarrow0.
\qquad\text{(6.6)}
```

En particulier, pour $s_n=\log n$, la coupe pertinente devient
microscopique sans qu'il soit nécessaire de la fixer artificiellement à
zéro.

Cette observation confirme l'intuition :

```math
\text{signal proportionnel croissant}
\quad\Longrightarrow\quad
\text{blocs critiques formés très tôt}.
\qquad\text{(6.7)}
```

Elle ne dit pas que la preuve d'almost/exact recovery se réduit à la
géométrie de ces blocs : l'exposant local de Hellinger--Chernoff et ses
queues restent à contrôler.

## 7. Quel rôle pour $\beta=0$ et Glauber ?

À $\beta=0$, la projection de coupe $B_0$ est vide et déterministe. Après
avoir marginalisé — ou simplement ignoré — le reste du dendrogramme, la
conditionnelle des spins est donc la vraie postérieure $\mu_A$. Une mise à
jour séquentielle de chaque singleton selon sa conditionnelle sous $\mu_A$
est le bain thermique de Glauber.

Cette extrémité est utile comme **raffinement local** :

1. construire une initialisation par une méthode globale ou hiérarchique ;
2. descendre vers des blocs de plus en plus petits ;
3. terminer par des heat baths mono-site, leave-one-out ou leave-block-out ;
4. contrôler leur erreur par l'exposant de Chernoff (2.4).

Trois réserves sont nécessaires.

- Couper à zéro dans un $D$ complet déjà conditionné n'est pas Glauber :
  les temps et facteurs ancêtres restent révélés. Il faut revenir à la
  projection $B_0$ et marginaliser le reste de $D$.
- Même sous $B_0$, il faut choisir un programme mono-site séquentiel ;
  une recoloration parallèle n'est pas la chaîne de Glauber.
- Tempérer réellement la postérieure jusqu'à une température infinie
  supprimerait l'information et briserait Nishimori ; ce n'est pas ce que
  fait le temps $\beta$.
- Sous la contrainte de bisection exacte, un flip mono-site quitte l'espace
  d'états. Il faut utiliser des swaps, une dynamique de Kawasaki ou des
  mises à jour par paires équilibrées.

## 8. Programme de preuve proposé

Une extension hiérarchique crédible suit quatre étages.

### Étape A — initialisation

Obtenir un overlap positif au-dessus de Kesten--Stigum, puis un overlap
tendant vers un lorsque $\lambda_n\to\infty$.

### Étape B — blocs contrôlables

Pour éviter la susceptibilité divergente exactement à la coupe, choisir
$\beta_n^-$ par

```math
d_nq_{\beta_n^-}=1-\varepsilon_n,
\qquad
\varepsilon_n\downarrow0
\quad\text{assez lentement}.
\qquad\text{(8.1)}
```

Calculer les messages exacts de ces blocs sans supprimer les facteurs
supérieurs.

### Étape C — raffinement local

Employer des données séparées ou un leave-one-out pour que l'erreur de
l'initialisation ne biaise pas le test du sommet courant. Les heat baths aux
feuilles donnent l'interprétation dynamique de cette étape.

### Étape D — fermeture des queues

Pour almost exact, montrer $P_{e,n}\to0$. Pour exact, montrer

```math
\sup_v P_{e,n}(v)=o(1/n)
\qquad\text{(8.2)}
```

ou obtenir une borne globale équivalente. C'est ici que l'affinité
Hellinger--Chernoff, et non le seul carré d'overlap, devient indispensable.

Ce programme est compatible avec la dynamique hiérarchique ; sa fermeture
au seuil optimal reste un problème de recherche.
