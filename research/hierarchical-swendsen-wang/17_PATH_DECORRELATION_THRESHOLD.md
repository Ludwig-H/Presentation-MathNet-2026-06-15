# Seuil de décorrélation le long du chemin hiérarchique

> **Statut dans le programme.** Les résultats PATH-FAC de cette note sont un
> benchmark factorisé. Leur interprétation exacte est désormais le corridor
> produit mono-bit sous prior de parités indépendant ; le fichier 20 traite
> sa tensorisation abstraite avec prior corrélé. Elle ne décrit pas le
> corridor collapsed multiport réel.

Cette note cherche précisément quand, en fonction de $`p`$, deux sommets
lointains perdent leur corrélation sous la dynamique descendante. Elle part de
l'oracle PATH-FAC du fichier 16, puis isole un critère suffisant qui reste
valable pour le transfert joint exact.

Le verdict comporte quatre régimes distincts.

1. **Critère exact dans PATH-FAC.** La probabilité de conserver la relation
   de paire tend vers $`1/2`$ si et seulement si l'atténuation cumulée
   $`A_L=-\sum_w\log\Gamma_w`$ tend vers l'infini.
2. **Interfaces bornées.** Pour tout $`p<1`$ fixé, un nombre divergent de
   buckets non triviaux de taille bornée force $`A_L\to\infty`$. Le seuil
   PATH-FAC est alors le bord dégénéré $`p=1`$.
3. **Interfaces logarithmiques.** Si les buckets critiques ont une taille
   régulière $`m_L\sim\alpha\log H_L`$, un seuil non trivial apparaît. Il est
   donné par l'expression suivante.

```math
p_{\mathrm{path}}(\alpha)
=
\frac{
1+q_\triangle
+(1-q_\triangle)\sqrt{1-e^{-2/\alpha}}
}2.
```

En dessous, PATH-FAC se décorrèle ; au-dessus, il conserve la relation. Le cas
d'égalité exige une correction $`\log\log H_L`$.

4. **Interfaces hétérogènes.** Si toutes les tailles critiques s'échappent
   vers l'infini, la quantité géométrique exacte au niveau exponentiel est
   décrite par la quantité suivante.

```math
\Phi_L(I)
=
\sum_{w\in\mathcal P_L}
m_w^{-1/2}e^{-I m_w}.
```

Son abscisse de transition, lorsqu'elle existe, se transporte en un seuil
explicite en $`p`$ par $`I=I_c(p)`$. Pour le vrai profil descendant, il faut
remplacer $`I_c(p)`$ par $`I(t_w;p)`$ nœud par nœud.

Ces conclusions montrent immédiatement qu'il n'existe pas de seuil fonction
de $`p`$ seul avant d'avoir établi la loi géométrique des tailles de coupe le
long du chemin de Kruskal.

## 1. Atténuation exacte

Pour une paire $`i,j`$, notons

```math
\mathcal P_L(i,j)
```

l'ensemble des nœuds des deux bras qui séparent exactement un des deux
sommets, et posons

```math
\Gamma_w
:=
\Gamma_{m_w}(t_w;p).
```

### Lemme 1.1 — critère nécessaire et suffisant, statut : établi dans PATH-FAC

Définissons

```math
\boxed{
A_L(p)
:=
-\sum_{w\in\mathcal P_L(i,j)}\log\Gamma_w.
}
```

Alors

```math
\boxed{
P_L^{\mathrm{FAC}}(p)
=
\frac12\left(1+e^{-A_L(p)}\right).
}
```

Par conséquent,

```math
\boxed{
P_L^{\mathrm{FAC}}(p)\longrightarrow\frac12
\quad\Longleftrightarrow\quad
A_L(p)\longrightarrow+\infty.
}
```

Plus quantitativement, pour $`0<\varepsilon<1/2`$,

```math
P_L^{\mathrm{FAC}}(p)-\frac12\le\varepsilon
\quad\Longleftrightarrow\quad
A_L(p)\ge\log\frac1{2\varepsilon}.
```

#### Preuve

La proposition 8.2 du fichier 16 donne

```math
2P_L^{\mathrm{FAC}}-1
=
\prod_w\Gamma_w.
```

Prendre le logarithme du produit donne $`e^{-A_L}`$. Les deux équivalences
suivent par monotonie de l'exponentielle.

### Interprétation : longueur de corrélation

Pour un chemin régulier dont tous les canaux ont la même fiabilité
$`\Gamma`$, posons

```math
\xi(\Gamma)
:=
\frac1{-\log\Gamma}.
```

Alors, pour un chemin de longueur $`H`$,

```math
P_H^{\mathrm{FAC}}-\frac12
=
\frac12e^{-H/\xi(\Gamma)}.
```

La paire perd donc sa corrélation lorsque $`H\gg\xi`$ et la conserve lorsque
$`H\ll\xi`$.

## 2. Un nombre divergent de petites coupes suffit

Le corollaire 8.4 du fichier 16 demandait une densité positive de petites
coupes. Cette hypothèse peut être fortement affaiblie.

Pour $`M\ge2`$, définissons

```math
N_{L,M}
:=
\#\{w\in\mathcal P_L(i,j):2\le m_w\le M\}.
```

### Théorème 2.1 — perte pour interfaces bornées, statut : établi dans PATH-FAC

Fixons $`p\in(1/2,1)`$ et $`M<\infty`$. Il existe
$`\delta(p,M)>0`$ tel que, uniformément pour $`t\in[0,1]`$ et
$`2\le m\le M`$,

```math
\Gamma_m(t;p)\le1-\delta(p,M).
```

Ainsi,

```math
P_L^{\mathrm{FAC}}(p)-\frac12
\le
\frac12(1-\delta(p,M))^{N_{L,M}}.
```

En particulier,

```math
N_{L,M}\longrightarrow\infty
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow\frac12.
```

Il n'est pas nécessaire que $`N_{L,M}`$ soit proportionnel à la longueur du
chemin ; une croissance aussi lente que $`\log\log H_L`$ suffit.

#### Preuve

Pour $`m\ge2`$ et $`p<1`$, la loi

```math
K=1+\mathrm{Bin}(m-1,s_p(t))
```

attribue une masse strictement positive à au moins un compte dont le log-odds
est fini. On a donc $`\Gamma_m(t;p)<1`$. Cette fonction est continue en
$`t`$ sur $`[0,1]`$. Son maximum sur la réunion finie des compacts
$`\{m\}\times[0,1]`$, $`2\le m\le M`$, est strictement inférieur à $`1`$.
Chaque petite coupe contribue alors au moins
$`-\log(1-\delta)`$ à $`A_L`$.

### Corollaire 2.2 — condition nécessaire pour ne pas perdre, statut : établi dans PATH-FAC

À $`p<1`$ fixé, si PATH-FAC ne se décorrèle pas, alors, pour tout $`M`$
fixé,

```math
N_{L,M}=O(1)
```

le long d'une sous-suite où la corrélation reste minorée. Hormis les buckets
$`m=1`$, qui sont parfaits, les tailles de coupe doivent donc s'échapper vers
l'infini.

### Contre-audit géométrique

- Sur un arbre physique, chaque bucket a $`m=1`$ et
  $`\Gamma_1=1`$ : le chemin ne perd jamais sa relation.
- Si une grille produit seulement un nombre divergent de buckets $`m=2`$, la
  relation disparaît pour tout $`p<1`$ fixé.
- Entre ces extrêmes, le seuil dépend de la croissance des interfaces.

Ces trois géométries ont le même paramètre de bruit $`p`$ et des conclusions
différentes. Un « seuil en $`p`$ » sans hypothèse géométrique serait donc faux.

## 3. Taux de décorrélation dans un environnement de chemin

Soit $`H_L=|\mathcal P_L(i,j)|`$. Supposons qu'il existe un taux
$`\kappa(p)\in[0,+\infty]`$ tel que

```math
\frac{A_L(p)}{H_L}
=
-\frac1{H_L}\sum_{w\in\mathcal P_L}\log\Gamma_w
\longrightarrow
\kappa(p).
```

Alors

```math
\boxed{
P_L^{\mathrm{FAC}}(p)-\frac12
=
\frac12\exp[-H_L\kappa(p)+o(H_L)].
}
```

Dans un environnement stationnaire ergodique de loi $`\pi_p`$ sur
$`(m,t)`$, le théorème ergodique suggère, sous intégrabilité,

```math
\kappa(p)
=
-\int\log\Gamma_m(t;p)\,\pi_p(dm,dt).
```

Si $`\pi_p(m\ge2)>0`$ et les tailles restent finies sous $`\pi_p`$, alors
$`\kappa(p)>0`$ pour tout $`p<1`$ : la corrélation décroît exponentiellement
en $`H_L`$. La quantité réellement à estimer sur la grille triangulaire est
donc la loi $`\pi_p`$, qui dépend elle-même de la géométrie de percolation et
de la sélection de la paire critique.

## 4. Grandes coupes critiques régulières

Plaçons-nous désormais sur la grille triangulaire et posons

```math
q_\triangle=2\sin(\pi/18),
\qquad
p_{\mathrm{SW}}=\frac{1+q_\triangle}{2}.
```

Pour $`p>p_{\mathrm{SW}}`$, rappelons

```math
h_c(p)
=
\frac{2p-1-q_\triangle}{1-q_\triangle},
```

```math
I_c(p)
:=
-\frac12\log(1-h_c(p)^2).
```

Considérons l'oracle régulier formé de $`H_L`$ buckets indépendants, tous
critiques et tous de taille $`m_L\to\infty`$. Posons

```math
r_L:=m_L\bmod2.
```

Le fichier 15 donne

```math
1-\Gamma_{m_L}^c(p)
\sim
\frac{2C_{r_L}(p)}{\sqrt{m_L}}e^{-m_LI_c(p)}.
```

### Théorème 4.1 — fenêtre logarithmique aiguë, statut : établi dans PATH-FAC

Définissons la coordonnée

```math
\boxed{
z_L(p)
:=
m_LI_c(p)-\log H_L+\frac12\log m_L.
}
```

Alors, dans PATH-FAC :

- si $`z_L(p)\to-\infty`$, alors
  $`P_L^{\mathrm{FAC}}(p)\to1/2`$ ;
- si $`z_L(p)\to+\infty`$, alors
  $`P_L^{\mathrm{FAC}}(p)\to1`$ ;
- si $`z_L(p)\to z\in\mathbb R`$ et si $`r_L=r`$ est fixé, alors la limite
  suivante vaut :

```math
\boxed{
P_L^{\mathrm{FAC}}(p)
\longrightarrow
\frac12\left[
1+\exp\bigl(-2C_r(p)e^{-z}\bigr)
\right].
}
```

#### Preuve

Écrivons $`d_L=1-\Gamma_{m_L}^c(p)`$. Comme $`d_L\to0`$,

```math
-H_L\log\Gamma_{m_L}^c(p)
=
H_Ld_L(1+o(1)).
```

L'équivalent aigu donne

```math
H_Ld_L
\sim
2C_{r_L}(p)e^{-z_L(p)}.
```

Le lemme 1.1 termine les trois cas.

### Proposition 4.2 — spectre de tailles hétérogène, statut : établi dans PATH-FAC critique

Supposons que tous les buckets soient critiques, mais autorisons des tailles
$`m_{L,w}`$ différentes. Pour $`I\ge0`$, posons

```math
\Phi_L(I)
:=
\sum_{w\in\mathcal P_L}
m_{L,w}^{-1/2}e^{-I m_{L,w}}.
```

Si

```math
\min_{w\in\mathcal P_L}m_{L,w}\longrightarrow+\infty,
```

alors, pour tout $`p\in(p_{\mathrm{SW}},1)`$ fixé, il existe deux constantes
$`0<c_p<C_p<+\infty`$ telles que, pour $`L`$ assez grand,

```math
\boxed{
c_p\Phi_L(I_c(p))
\le
A_L(p)
\le
C_p\Phi_L(I_c(p)).
}
```

Par conséquent,

```math
\Phi_L(I_c(p))\longrightarrow+\infty
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow\frac12,
```

et

```math
\Phi_L(I_c(p))\longrightarrow0
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow1.
```

#### Preuve

Pour chacune des deux parités de $`m`$, l'équivalent du fichier 15 donne

```math
1-\Gamma_m^c(p)
\sim
2C_{m\bmod2}(p)m^{-1/2}e^{-mI_c(p)}.
```

Les deux constantes limites sont strictement positives et finies. Il existe
donc $`m_0(p)`$ et des constantes $`0<c'_p<C'_p<+\infty`$ qui encadrent le
rapport pour tout entier $`m\ge m_0(p)`$. Comme la taille minimale diverge,
cet encadrement vaut simultanément pour tous les nœuds du chemin. Enfin, si
$`d=1-\Gamma`$ et $`d\le1/2`$,

```math
d\le-\log(1-d)\le2d.
```

La somme de ces inégalités donne l'encadrement annoncé.

Cette proposition montre pourquoi une taille moyenne ne suffit pas : la
somme est exponentiellement dominée par les plus petites interfaces.

### Proposition 4.3 — spectre décoré par les niveaux descendants, statut : établi à $`p`$ fixé

Fixons $`p\in(p_{\mathrm{SW}},1)`$ et autorisons maintenant les niveaux
$`0\le t_{L,w}\le\beta_c(p)`$. Rappelons

```math
I(t;p)
=
\log\cosh\frac{u_p(1-t)}2.
```

Définissons

```math
\boxed{
\Phi_L^{\mathrm{desc}}(p)
:=
\sum_{w\in\mathcal P_L}
m_{L,w}^{-1/2}
e^{-m_{L,w}I(t_{L,w};p)}.
}
```

Si $`\min_w m_{L,w}\to+\infty`$, alors il existe
$`0<c_p<C_p<+\infty`$ tels que, pour $`L`$ assez grand,

```math
\boxed{
c_p\Phi_L^{\mathrm{desc}}(p)
\le A_L(p)\le
C_p\Phi_L^{\mathrm{desc}}(p).
}
```

La divergence de $`\Phi_L^{\mathrm{desc}}(p)`$ force donc la limite $`1/2`$
dans PATH-FAC, tandis que sa convergence vers zéro force la limite $`1`$.

#### Preuve de l'uniformité

Sur $`0\le t\le\beta_c(p)`$, le paramètre

```math
a_p(t)=u_p(1-t)
```

reste dans le compact $`[a_c(p),u_p]\subset(0,+\infty)`$. Dans la preuve par
point selle du fichier 16, les deux préfacteurs

```math
C_r(t;p)
=
\frac1{2s_p(t)\sqrt{2\pi}}
\sum_{j\in\mathbb Z}
\frac1{\cosh(a_p(t)(j+\varepsilon_r))}
```

sont donc continus, strictement positifs et uniformément bornés. Les termes
centraux de la somme harmonique ont une approximation de Stirling uniforme
sur ce compact. Ses queues sont dominées par une série géométrique uniforme,
car $`a_p(t)\ge a_c(p)>0`$. L'équivalent

```math
1-\Gamma_m(t;p)
\sim
2C_{m\bmod2}(t;p)m^{-1/2}e^{-mI(t;p)}
```

est ainsi uniforme en $`t\in[0,\beta_c(p)]`$. Les inégalités
$`d\le-\log(1-d)\le2d`$ terminent comme dans la proposition 4.2.

Cette uniformité est à $`p`$ fixé. Elle dégénère lorsque
$`p\downarrow p_{\mathrm{SW}}`$, car $`a_c(p)\downarrow0`$, et elle ne couvre
pas non plus la fenêtre simultanée $`p\uparrow1`$ du paragraphe 6.

### Lemme 4.4 — épaisseur de la bande descendante pertinente, statut : établi

Écrivons $`t=\theta\beta_c(p)`$ et fixons
$`p\in(p_{\mathrm{SW}},1)`$. Alors

```math
a_p(\theta\beta_c)
=
a_c(p)+b_c(p)(1-\theta),
```

où

```math
b_c(p)
:=
u_p\beta_c(p)
=
-\log\left(1-\frac{q_\triangle}{p}\right)>0.
```

Comme la dérivée de $`a\mapsto\log\cosh(a/2)`$ vaut
$`\tanh(a/2)/2`$,

```math
\boxed{
I(\theta\beta_c;p)
=
I_c(p)
+g_c(p)(1-\theta)
+O((1-\theta)^2),
}
```

avec

```math
g_c(p)
:=
\frac{h_c(p)b_c(p)}2>0.
```

Pour une coupe de taille $`m`$, le rapport entre son poids spectral au niveau
$`\theta\beta_c`$ et son poids critique vaut donc

```math
\exp\left[
-m g_c(p)(1-\theta)
+O(m(1-\theta)^2)
\right].
```

Ainsi, à $`p>p_{\mathrm{SW}}`$ fixé, seuls les niveaux vérifiant

```math
1-\theta=O(m^{-1}),
```

ou de façon équivalente $`\beta_c-t=O(m^{-1})`$, contribuent au même ordre
exponentiel qu'un bucket critique de même taille. Si
$`m(1-\theta)\to+\infty`$, leur contribution relative tend vers zéro.

Ce lemme rend précise l'expression « fusion très proche du seuil » : la
fenêtre pertinente dépend de la taille de l'interface, pas seulement de la
taille de la boîte.

### Corollaire 4.5 — seuil à hauteur relative régulière, statut : établi dans PATH-FAC

Supposons que les $`H_L`$ buckets aient
$`m_L\sim\alpha\log H_L`$ et le même niveau relatif
$`t_L=\theta\beta_c(p)`$, avec $`0\le\theta\le1`$. Posons

```math
I_\theta(p):=I(\theta\beta_c(p);p).
```

Les régimes stricts sont

```math
\alpha I_\theta(p)<1
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow\frac12,
```

```math
\alpha I_\theta(p)>1
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow1.
```

La fonction $`I_\theta(p)`$ est strictement croissante en $`p`$. Si

```math
\alpha I_\theta(p_{\mathrm{SW}})<1,
```

il existe donc un unique seuil intérieur $`p_{\mathrm{path}}(\alpha,\theta)`$.
En posant

```math
a_\alpha
:=
2\,\mathrm{arcosh}(e^{1/\alpha}),
```

il est l'unique solution dans $`(p_{\mathrm{SW}},1)`$ de

```math
\boxed{
(1-\theta)\log p
+\theta\log(p-q_\triangle)
-\log(1-p)
=
a_\alpha.
}
```

Si $`\alpha I_\theta(p_{\mathrm{SW}})\ge1`$, il n'existe pas de phase de
perte stricte dans l'intervalle accessible $`p>p_{\mathrm{SW}}`$ pour cet
oracle régulier. Enfin,

```math
p_{\mathrm{path}}(\alpha,0)
=
\frac{1+\sqrt{1-e^{-2/\alpha}}}{2},
```

lorsque cette valeur dépasse $`p_{\mathrm{SW}}`$, tandis que

```math
p_{\mathrm{path}}(\alpha,1)
=
p_{\mathrm{path}}(\alpha).
```

À $`\alpha`$ fixé, le seuil est croissant en $`\theta`$. Le modèle où tous
les descendants sont artificiellement placés à $`\beta_c`$ est donc
l'enveloppe la plus favorable à la **décorrélation**, et non une description
du vrai chemin. Le seul fait que le LCA soit critique ne justifie pas cette
substitution.

### Corollaire 4.6 — seuil spectral géométrique, statut : conditionnel à une abscisse aiguë

Définissons les deux abscisses étendues

```math
I_-
:=
\sup\{I\ge0:\Phi_L(I)\to+\infty\},
```

```math
I_+
:=
\inf\{I\ge0:\Phi_L(I)\to0\}.
```

La monotonie de $`\Phi_L`$ donne $`I_-\le I_+`$. Si elles coïncident en une
valeur $`I_*\in(0,+\infty)`$, le seuil de premier ordre est

```math
\boxed{
p_*
=
\frac{
1+q_\triangle
+(1-q_\triangle)\sqrt{1-e^{-2I_*}}
}2.
}
```

Plus précisément, $`p_{\mathrm{SW}}<p<p_*`$ entraîne la perte dans PATH-FAC,
et $`p>p_*`$ entraîne la conservation. Au point $`p=p_*`$, il faut connaître
la limite de $`\Phi_L(I_*)`$, les préfacteurs de parité et les corrections
sous-exponentielles.

Sans coïncidence de $`I_-`$ et $`I_+`$, la géométrie peut produire une bande
de crossover ou des sous-suites différentes ; annoncer un seuil unique en
$`p`$ serait alors injustifié.

### Corollaire 4.7 — seuil de premier ordre en $`p`$, statut : établi dans l'oracle régulier

Supposons

```math
\frac{m_L}{\log H_L}\longrightarrow\alpha\in(0,+\infty).
```

Pour les inégalités strictes :

```math
\alpha I_c(p)<1
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow\frac12,
```

```math
\alpha I_c(p)>1
\quad\Longrightarrow\quad
P_L^{\mathrm{FAC}}(p)\longrightarrow1.
```

L'équation $`\alpha I_c(p)=1`$ possède l'unique solution

```math
\boxed{
p_{\mathrm{path}}(\alpha)
=
\frac{
1+q_\triangle
+(1-q_\triangle)\sqrt{1-e^{-2/\alpha}}
}{2}.
}
```

Ainsi, sous cette hypothèse géométrique régulière,

```math
p_{\mathrm{SW}}<p<p_{\mathrm{path}}(\alpha)
\quad\Longrightarrow\quad
\text{perte de corrélation},
```

```math
p>p_{\mathrm{path}}(\alpha)
\quad\Longrightarrow\quad
\text{conservation dans PATH-FAC}.
```

Au point d'égalité, le seul rapport $`m_L/\log H_L`$ ne suffit pas : la
coordonnée $`z_L`$ montre que le terme $`\frac12\log m_L`$ et les corrections
d'ordre $`\log\log H_L`$ décident de la limite.

Ce corollaire est le cas particulier de la proposition 4.2 : si tous les
buckets ont la taille $`m_L\sim\alpha\log H_L`$, alors

```math
\Phi_L(I)
=
H_Lm_L^{-1/2}e^{-I m_L},
```

et son abscisse de transition vaut $`I_*=1/\alpha`$.

Lorsque $`\alpha\to\infty`$,

```math
p_{\mathrm{path}}(\alpha)
=
p_{\mathrm{SW}}
+\frac{1-q_\triangle}{\sqrt{2\alpha}}
+O(\alpha^{-3/2}).
```

## 5. Audit de la constante de Nishimori

La baseline information--percolation du dossier est

```math
p_{\mathrm{info}}
=
\frac{1+\sqrt{q_\triangle}}2
=
0.794659275831\ldots
```

Elle correspond, dans l'oracle régulier, au coefficient

```math
\boxed{
\alpha_{\mathrm{info}}
:=
\frac1{I_c(p_{\mathrm{info}})}
=
13.521628164595\ldots
}
```

Comme $`p_{\mathrm{path}}(\alpha)`$ est décroissante, PATH-FAC ne pourrait
améliorer cette baseline que si la géométrie donnait

```math
\alpha<\alpha_{\mathrm{info}}.
```

Cette comparaison est nécessaire avant toute revendication de meilleure
borne de weak recovery.

Soit

```math
p_{\mathrm N}^{(0)}
=
0.835805792367\ldots
```

la racine de l'équation triangulaire de Nishimori--Ohzeki étudiée dans le
fichier 13. La valeur de $`\alpha`$ qui force artificiellement
$`p_{\mathrm{path}}(\alpha)=p_{\mathrm N}^{(0)}`$ est

```math
\boxed{
\alpha_{\mathrm N}
:=
\frac1{I_c(p_{\mathrm N}^{(0)})}
=
7.053596192884\ldots
}
```

Ce calcul donne une cible géométrique très précise : la voie PATH-FAC
retrouverait le nombre de Nishimori au premier ordre si les interfaces
pertinentes vérifiaient

```math
m_L
\sim
7.053596192884\ldots\,\log H_L.
```

### Contre-audit

Cette coïncidence n'est pas une dérivation du seuil de Nishimori. La fonction
$`\alpha\mapsto p_{\mathrm{path}}(\alpha)`$ parcourt tout l'intervalle
$`(p_{\mathrm{SW}},1)`$ : chaque valeur de $`p`$ peut être reproduite en
choisissant $`\alpha=1/I_c(p)`$. Il faut donc démontrer le coefficient
géométrique $`\alpha_{\mathrm N}`$ à partir de la percolation critique, et non
l'ajuster après coup.

Pour obtenir une obstruction au moins jusqu'à $`p_{\mathrm N}^{(0)}`$, il
faudrait en particulier $`\alpha\le\alpha_{\mathrm N}`$, puis transporter la
décorrélation de PATH-FAC à la dynamique jointe. Au point d'égalité, la
fenêtre $`z_L`$ doit encore être contrôlée.

### Audit géométrique de la littérature planaire

Les résultats existants sur la percolation triangulaire proche-critique ne
ferment pas cette étape.

- Garban--Pete--Schramm construisent la mesure limite des points pivotaux,
  avec la normalisation critique gouvernée par l'exposant quatre bras. Cette
  mesure repère **où** une modification proche-critique change une connexion
  macroscopique ; elle ne donne pas la cardinalité $`m_w`$ de toute la coupe
  non marquée entre les deux fils d'un nœud de Kruskal.
- Leur construction de la limite d'échelle du MST encode la topologie et les
  niveaux proche-critiques des branchements macroscopiques. Elle ne fournit
  pas, dans les énoncés disponibles, la loi jointe
  $`(H_L,(m_w,t_w)_{w\in\mathcal P_L})`$ sous la loi de Palm imposant une paire
  lointaine dans la composante géante et un LCA critique.
- Les outlets de l'invasion percolation sont des arêtes record séparant des
  échelles successives. Leur nombre n'est ni $`H_L`$ ni $`N_{L,M}`$, et leur
  multiplicité n'est pas la taille $`m_w`$ d'un bucket. En particulier, un
  outlet unique correspond potentiellement à $`m_w=1`$, canal parfait qui ne
  contribue aucune atténuation.

Par conséquent, ni l'exposant pivotal $`3/4`$, ni une estimation du nombre
d'outlets, ni l'existence de la limite du MST ne justifient l'ansatz
$`m_L\sim\alpha\log H_L`$. Déduire $`\alpha=\alpha_{\mathrm N}`$ de l'un de
ces résultats sans théorème de transfert serait une erreur de catégorie.
Les références et leurs portées exactes sont récapitulées dans
[l'état de l'art](LITERATURE.md#bande-critique-pivots-et-sprinkling).

## 6. Buckets de taille fixe et fenêtre $`p\uparrow1`$

Le théorème 2.1 implique que, pour $`m\ge2`$ fixé et $`p<1`$ fixé, un chemin
régulier de longueur $`H_L\to\infty`$ se décorrèle toujours. Un crossover
subsiste si $`p=p_L\uparrow1`$.

Au niveau critique, posons

```math
\varepsilon_p
:=
1-s_c(p)
=
\frac{1-p}{1-q_\triangle}.
```

### Lemme 6.1 — déficit aigu à taille fixe, statut : établi

Pour $`m\ge2`$ fixé,

```math
\boxed{
1-\Gamma_m^c(p)
\sim
D_m\varepsilon_p^{d_m},
\qquad
p\uparrow1,
}
```

où

```math
d_m=\left\lceil\frac m2\right\rceil
```

et

```math
\boxed{
D_{2r}=\binom{2r-1}{r-1},
\qquad
D_{2r+1}=4\binom{2r}{r-1}
\quad(r\ge1).
}
```

#### Preuve

Avec $`P_+`$ et $`P_-`$ les deux lois de compte symétriques,

```math
1-\Gamma_m^c
=
2\sum_{k=0}^m
\frac{P_+(k)P_-(k)}{P_+(k)+P_-(k)}.
```

Pour $`s_c=1-\varepsilon`$ et $`1\le k\le m-1`$,

```math
P_+(k)
=
\binom{m-1}{k-1}
(1-\varepsilon)^{k-1}\varepsilon^{m-k},
```

```math
P_-(k)
=
\binom{m-1}{k}
(1-\varepsilon)^{m-k-1}\varepsilon^k.
```

Aux extrémités, $`P_+(0)=P_-(m)=0`$ et les termes harmoniques sont nuls.

Si $`m=2r`$, seul le compte central $`k=r`$ contribue au premier ordre ;
les deux masses y sont égales et donnent
$`D_{2r}=\binom{2r-1}{r-1}`$. Si $`m=2r+1`$, les deux comptes centraux
$`r,r+1`$ contribuent chacun, et le facteur extérieur $`2`$ donne
$`D_{2r+1}=4\binom{2r}{r-1}`$. Tous les autres comptes portent une puissance
strictement plus grande de $`\varepsilon`$.

### Corollaire 6.2 — fenêtre de décorrélation à taille fixe, statut : établi dans PATH-FAC

Pour $`H_L`$ canaux critiques identiques de taille $`m`$, posons

```math
Y_L
:=
H_LD_m
\left(\frac{1-p_L}{1-q_\triangle}\right)^{d_m}.
```

Alors :

- si $`Y_L\to+\infty`$, $`P_L^{\mathrm{FAC}}\to1/2`$ ;
- si $`Y_L\to0`$, $`P_L^{\mathrm{FAC}}\to1`$ ;
- si $`Y_L\to\lambda\in(0,+\infty)`$, alors
  $`P_L^{\mathrm{FAC}}\longrightarrow\frac12(1+e^{-\lambda})`$.

Le crossover est donc

```math
\boxed{
1-p_L
\asymp
(1-q_\triangle)H_L^{-1/d_m}.
}
```

La longueur de corrélation correspondante vérifie

```math
\boxed{
\xi_m^c(p)
\sim
\frac1{D_m}
\left(\frac{1-q_\triangle}{1-p}\right)^{d_m}.
}
```

### Niveaux descendants sous le LCA critique

Écrivons $`t=\theta\beta_c(p)`$, $`0\le\theta\le1`$. Comme

```math
e^{-u_p\beta_c(p)}=1-\frac{q_\triangle}{p},
```

on a exactement

```math
1-s_p(\theta\beta_c)
=
\frac{1-p}
{1-p+p(1-q_\triangle/p)^\theta},
```

et donc, lorsque $`p\uparrow1`$,

```math
1-s_p(\theta\beta_c)
\sim
\frac{1-p}{(1-q_\triangle)^\theta}.
```

Pour un bucket fixé de taille $`m`$ à hauteur relative $`\theta`$,

```math
1-\Gamma_m(\theta\beta_c;p)
\sim
D_m
\left[
\frac{1-p}{(1-q_\triangle)^\theta}
\right]^{d_m}.
```

Les niveaux descendants plus précoces, $`\theta<1`$, ont donc une longueur
de corrélation asymptotiquement plus grande que le nœud critique. Le bon
paramètre cumulé d'un chemin hétérogène est la somme de ces déficits, et non
le seul niveau maximal $`\beta_c`$.

## 7. Critère exact pour la dynamique jointe

Les seuils précédents sont exacts dans PATH-FAC. La proposition 10.2 du
fichier 16 donne cependant un moyen rigoureux d'obtenir la décorrélation sans
supposer l'indépendance.

Dans cette section, $`O,D`$ sont fixés ; les lois, opérateurs et normes sont
donc conditionnels à cet environnement.

Soit $`\lambda_r`$ la loi de l'état de frontière $`X_r`$ juste avant la
mise à jour $`r`$, et

```math
(\mathcal T_rf)(x)
=
\mathbb E[(-1)^{A_r\chi_r(i,j)}f(X_{r+1})\mid X_r=x,O,D]
```

l'opérateur tordu exact. Posons

```math
\kappa_r
:=
\|\mathcal T_r\|_{L^2(\lambda_{r+1})\to L^2(\lambda_r)}.
```

Par Jensen, $`0\le\kappa_r\le1`$.

### Proposition 7.1 — critère de contraction jointe, statut : établi

La corrélation exacte $`c_{ij}^{\mathrm{joint}}`$ après le balayage vérifie

```math
\boxed{
|c_{ij}^{\mathrm{joint}}|
\le
\prod_{r=1}^{H_L}\kappa_r.
}
```

Ainsi,

```math
\sum_{r=1}^{H_L}-\log\kappa_r\longrightarrow+\infty
\quad\Longrightarrow\quad
\mathbb P(\text{relation conservée})\longrightarrow\frac12.
```

#### Preuve

La récursion tordue donne

```math
c_{ij}^{\mathrm{joint}}
=
\lambda_1\mathcal T_1\cdots\mathcal T_{H_L}\mathbf1.
```

Cauchy--Schwarz, puis la sous-multiplicativité des normes, donnent

```math
|c_{ij}^{\mathrm{joint}}|
\le
\|\mathcal T_1\cdots\mathcal T_{H_L}\mathbf1\|_{L^2(\lambda_1)}
\le
\prod_r\kappa_r,
```

car $`\|\mathbf1\|_{L^2(\lambda_{H_L+1})}=1`$. L'inégalité
$`\kappa_r\le1`$ vient de

```math
|\mathbb E[Zf(X_{r+1})\mid X_r]|^2
\le
\mathbb E[f(X_{r+1})^2\mid X_r],
\qquad |Z|=1.
```

### Verrou restant

En général,

```math
\kappa_r\ne\Gamma_{m_r}(t_r;p).
```

La norme $`\kappa_r`$ voit les états rares de frontière et les dépendances
avec le futur, tandis que $`\Gamma_m`$ est une fiabilité locale moyenne. Elle
peut même valoir $`1`$ à cause d'un état déterministe alors que la fiabilité
moyenne est strictement inférieure à $`1`$. Le fichier 21 établit la
contraction de bloc sur cactus. Le résultat décisif restant sur les bandes
est une borne

```math
\|\mathcal T_{r:r+\ell}\|\le e^{-\mathfrak a_{r,\ell}(p)}
```

dont les atténuations $`\mathfrak a_{r,\ell}`$ s'additionnent le long du
chemin.

## 8. Audit et contre-audit

| Affirmation | Statut | Conclusion |
|---|---|---|
| $`P_L^{\mathrm{FAC}}\to1/2`$ ssi $`A_L\to\infty`$ | Établi dans PATH-FAC | critère exact |
| Une densité positive de petites coupes est nécessaire | Faux | un nombre divergent suffit |
| Pour des tailles bornées, il existe un seuil fixe $`p<1`$ | Faux dans PATH-FAC | toute valeur fixe $`p<1`$ se décorrèle |
| La taille moyenne des grandes coupes détermine le seuil | Faux en général | utiliser $`\Phi_L(I)`$, dominée par les petites tailles |
| $`m_L\sim\alpha\log H_L`$ produit un seuil non trivial | Établi dans l'oracle régulier | $`p_{\mathrm{path}}(\alpha)`$ explicite |
| Le LCA critique permet de placer tous les descendants à $`\beta_c`$ | Faux | seuls les niveaux dans $`\beta_c-t=O(1/m)`$ ont le même poids |
| L'égalité $`\alpha I_c(p)=1`$ décide seule du bord | Faux | correction $`\log\log H_L`$ et parité |
| La constante de Nishimori est retrouvée | Seulement si $`\alpha=7.053596\ldots`$ | exigence géométrique non démontrée |
| Les descendants précoces sont sans effet | Faux | ils changent la longueur de corrélation |
| Les $`\Gamma_w`$ se multiplient dans la dynamique complète | Non établi | remplacer par les opérateurs tordus |
| Une contraction tordue sommable implique la perte réelle | Établi | nouvelle cible rigoureuse |

## 9. Programme géométrique et dynamique

Pour transformer ces seuils en résultat sur le GSBM triangulaire, il faut
maintenant mesurer ou démontrer les objets suivants sous la loi de Palm d'une
paire critique lointaine.

1. La longueur $`H_L`$ des deux bras vers le LCA.
2. Les comptes $`N_{L,M}=\#\{w:2\le m_w\le M\}`$. Si l'un d'eux diverge pour
   un $`M`$ fixé, PATH-FAC est déjà décorrélé pour tout $`p<1`$ fixé.
3. Le profil empirique $`\{(m_w,t_w/\beta_c)\}_{w\in\mathcal P_L}`$.
4. Dans le régime de grandes interfaces, la fonction de partition
   $`\Phi_L(I)=\sum_w m_w^{-1/2}e^{-Im_w}`$, puis, dans le sous-cas régulier,
   le rapport $`m_w/\log H_L`$ et ses fluctuations d'ordre
   $`\log\log H_L`$.
5. Sur cactus et bandes, les normes de blocs tordus $`\kappa_r`$ et leur
   comparaison numérique avec $`\Gamma_{m_r}(t_r;p)`$.

Le premier test géométrique doit donc être formulé directement sur la
hiérarchie marquée. Pour chaque taille de boîte $`L`$, sous le conditionnement
de paire critique, enregistrer

```math
H_L,
\qquad
N_{L,M},
\qquad
\sum_{w\in\mathcal P_L}-\log\Gamma_{m_w}(t_w;p),
```

et non remplacer ces variables par un nombre de pivots ou d'outlets. Les
deux diagnostics décisifs sont :

```math
N_{L,M}\longrightarrow\infty
```

pour au moins un $`M`$ fixé, ce qui force la perte dans PATH-FAC pour tout
$`p<1`$ fixé, ou au contraire l'échappement de toutes les tailles et une loi
d'échelle de $`m_w`$ assez précise pour déterminer $`A_L(p)`$.

La réponse actuelle est donc précise : la perte de corrélation est gouvernée
par une **atténuation cumulée**. Un seuil explicite en $`p`$ apparaît une fois
fixée la croissance géométrique des interfaces ; sans cette donnée, parler
d'un seuil universel en $`p`$ confond le bruit du canal et la géométrie du
chemin.

## 10. Passage au critère global : correction de second moment

Pour une probabilité de relation d'une paire, le premier moment signé du
transfert est le bon objet. Pour interdire la weak recovery, il faut cependant
contrôler

```math
\mathbb E[H_S(I_n,J_n)^2],
```

et non seulement $`\mathbb E H_S(I_n,J_n)`$. Deux environnements portant des
corrélations $`+1`$ et $`-1`$ peuvent annuler le premier moment tout en
restant parfaitement informatifs. Le
[transfert répliqué sous Palm critique](18_CRITICAL_PALM_REPLICATED_TRANSFER.md)
établit la globalisation spectrale, définit le noyau répliqué partageant le
même environnement et documente l'ancienne domination HF-S2. Celle-ci est
fausse en général multiport ; le sweep complet doit être contrôlé directement
aux rangs réalisés.
