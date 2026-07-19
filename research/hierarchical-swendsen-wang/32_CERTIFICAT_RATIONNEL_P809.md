# Certificat rationnel less-noisy à $`p=0.809`$

> [!NOTE]
> Ce certificat reste un jalon exact et un contre-test lisible. La meilleure
> borne du dossier est désormais le [certificat
> P809439](34_CERTIFICAT_RATIONNEL_P809439.md), qui porte
> $`p_{\mathrm{WR}}\ge0.809439`$.

## 1. Énoncé exact

Fixons

```math
p_1=\frac{809}{1000},
\qquad
q_1=2p_1-1=\frac{309}{500},
```

et le canal d'effacement triangulaire

```math
(a_1,s_1,e_1)
=
\frac1{5000}(1660,559,1663).
\qquad\text{(1.1)}
```

Le même argument exact que pour le
[jalon A0](31_CERTIFICAT_RATIONNEL_A0.md) donne le résultat suivant.

**Théorème P809.** Pour tout a priori $`\mu`$ sur les quatre états relatifs
et toute fonction $`f`$,

```math
Q_{E_{a_1,s_1,e_1}}(\mu,f)
-Q_{Y_{q_1}}(\mu,f)
\ge
\frac1{50000}\,\mathrm{Var}_\mu(f).
\qquad\text{(1.2)}
```

Le certificat est exhaustif : il utilise des fractions exactes, quatre
suites de Sturm et une identité exacte de dominance diagonale. Il ne contient
ni maillage du simplexe ni calcul flottant.

Le lemme se relève au GSBM triangulaire complet.

**Théorème global P809.** Soit $`\mathbb T_L`$ le tore triangulaire de côté
$`L`$, muni de spins i.i.d. uniformes, et soit $`n_L=L^2`$. Pour tout
$`p\in[1/2,809/1000]`$ et tout estimateur, éventuellement randomisé,
$`\widehat X^{(L)}\in\{\pm1\}^{n_L}`$ construit depuis les observations,

```math
\lim_{L\to\infty}
\mathbb E\left[
\left(
\frac1{n_L}\sum_{i\in\mathbb T_L}X_i\widehat X_i^{(L)}
\right)^2
\right]
=0.
\qquad\text{(1.3)}
```

Il n'y a donc pas de weak recovery, au sens du
[fichier 03](03_HIERARCHICAL_WEAK_RECOVERY.md), sur cet intervalle. En
particulier, la borne d'impossibilité rigoureuse atteint
$`p=809/1000>0.8`$.

## 2. Marges de Chayes--Lei

Les contraintes du canal auxiliaire ont toutes une marge stricte :

```math
a_1+3s_1+e_1=1,
\qquad
e_1-a_1=\frac3{5000},
\qquad
1-(2a_1+3s_1)=\frac3{5000}.
\qquad\text{(2.1)}
```

L'inégalité d'association vaut

```math
a_1e_1-2s_1^2
=
\frac{1067809}{12500000}>0.
\qquad\text{(2.2)}
```

Enfin, la seconde condition auxiliaire de Chayes--Lei est vérifiée
exactement. Comme les deux membres sont positifs,

```math
a_1+e_1>\frac{2\sqrt2}{3+2\sqrt2}=6\sqrt2-8
```

est équivalente à l'inégalité rationnelle

```math
(a_1+e_1+8)^2-72
=
\frac{76882329}{25000000}>0.
\qquad\text{(2.3)}
```

## 3. Les quatre polynômes exacts

Le profil du canal physique s'écrit

```math
c_{q_1}(t)
=
\frac{381924t(536443-381924t)}{H_1(t)},
\qquad
H_1(t)=(36481+618000t)(654481-618000t).
\qquad\text{(3.1)}
```

On a $`H_1(t)>0`$ sur $`[0,1]`$. Après soustraction de la marge de
(1.2), posons

```math
\bar a_1=a_1-\frac1{50000}=\frac{16599}{50000},
\qquad
d_1(t)=\frac{\bar a_1-c_{q_1}(t)}{t}.
\qquad\text{(3.2)}
```

Dans le secteur non polarisé, le numérateur primitif de
$`\bar a_1+4s_1t-c_{q_1}(t)`$ est

```math
P_{\mathrm{np}}(t)
=396319738471239
-3370596266968040t
+9493561252800000t^2
-8539820640000000t^3.
\qquad\text{(3.3)}
```

Dans le secteur polarisé, les trois numérateurs primitifs nécessaires sont

```math
P_{\mathrm{tail}}(t)
=45640313188927373653
-297658582360966800000t
+485484077838348000000t^2,
\qquad\text{(3.4)}
```

```math
\begin{aligned}
P_{\mathrm{dec}}(t)
={}&3154192724499694324012093\\
&+100909346530592989224000000t\\
&-656024362891036720509600000t^2\\
&+994139601778209600000000000t^3\\
&-121418809934342400000000000t^4,
\end{aligned}
\qquad\text{(3.5)}
```

et

```math
P_{\mathrm{off}}(t)
=529787256879229
-2165830919071239t
+2723251793400000t^2
-953740612800000t^3.
\qquad\text{(3.6)}
```

Ici $`P_{\mathrm{tail}}`$ est obtenu après extraction du facteur positif
$`1-2t`$ dans $`d_1(t)-d_1(1/2)`$,
$`P_{\mathrm{dec}}`$ est le numérateur de $`-d_1'(t)`$, et
$`P_{\mathrm{off}}`$ celui de
$`d_1(t)+s_1/[t(1-t)]`$.

## 4. Comptage de Sturm

Les variations aux deux extrémités et les nombres de racines ouvertes sont :

| polynôme | intervalle | variations | racines | témoin positif |
|---|---:|---:|---:|---:|
| $`P_{\mathrm{np}}`$ | $`[0,1/2]`$ | $`2\to2`$ | $`0`$ | $`P_{\mathrm{np}}(1/4)=13583552529229`$ |
| $`P_{\mathrm{tail}}`$ | $`[0,1/2]`$ | $`1\to1`$ | $`0`$ | $`P_{\mathrm{tail}}(1/4)=1568422463582423653`$ |
| $`P_{\mathrm{dec}}`$ | $`[1/2,1]`$ | $`1\to1`$ | $`0`$ | $`P_{\mathrm{dec}}(3/4)=90807472665630430955362093`$ |
| $`P_{\mathrm{off}}`$ | $`[1/2,1]`$ | $`2\to2`$ | $`0`$ | $`P_{\mathrm{off}}(3/4)=139535521353199/4`$ |

Aucun polynôme ne s'annule à une extrémité. Les quatre signes sont donc
strictement positifs sur les intervalles indiqués.

## 5. Matrice polarisée et faces

Pour $`\mu_0>1/2`$, l'élimination de la coordonnée centrée $`g_0`$ donne la
même matrice $`3\times3`$ que dans le certificat A0. Ses termes hors
diagonale sont positifs par (3.6) : avec les notations de la
[section 5 du certificat A0](31_CERTIFICAT_RATIONNEL_A0.md#5-positivité-des-termes-hors-diagonale),
$`D_k\le\mu_0(1-\mu_0)`$. Pour
$`\{i,j,k\}=\{1,2,3\}`$, l'annulation exacte reste

```math
M_{ii}-M_{ij}-M_{ik}
=d_1(\mu_i)-d_1(\mu_0)\ge0.
\qquad\text{(5.1)}
```

Les signes de (3.4)--(3.5) prouvent la dernière inégalité. La matrice est
donc symétrique, à diagonale dominante et positive semi-définie. Les faces
du simplexe suivent par restriction au support, ou par continuité des formes
quadratiques des canaux finis. Sur la face binaire uniforme, la marge
renforcée est

```math
\bar a_1+2s_1-\frac{2q_1^2}{1+q_1^2}
=
\frac{49016699}{17274050000}>0.
\qquad\text{(5.2)}
```

Ceci démontre (1.2).

## 6. Relèvement du quotient aux trois spins

Les triangles montants forment une partition des arêtes de la grille. Pour
un tel triangle $`T=\{x_0,x_1,x_2\}`$, posons

```math
\rho_T(X)
=
(X_{x_0}X_{x_1},X_{x_0}X_{x_2})
\in\{\pm1\}^2.
\qquad\text{(6.1)}
```

Le facteur physique $`P_{T,q}`$ et le facteur auxiliaire $`Q_{T,E}`$ sont
les compositions respectives de $`\rho_T`$ avec $`Y_q`$ et
$`E_{a_1,s_1,e_1}`$. Le théorème P809 et le critère $`\chi^2`$ less-noisy de
[Makur--Polyanskiy](https://doi.org/10.1109/TIT.2018.2839743) donnent

```math
E_{a_1,s_1,e_1}\succeq_{\mathrm{ln}}Y_{q_1}.
\qquad\text{(6.2)}
```

Cette comparaison se relève sans hypothèse supplémentaire. En effet, pour
toute loi de $`(U,X_{x_0},X_{x_1},X_{x_2})`$, sa poussée en avant par
$`\rho_T`$ est une loi admissible de $`(U,\rho_T(X))`$ dans la définition de
l'ordre less-noisy. Ainsi

```math
Q_{T,E}\succeq_{\mathrm{ln}}P_{T,q_1}
```

comme canaux ayant les trois spins pour entrée.

Pour $`q\in[0,q_1]`$, écrivons chaque bruit du triangle physique sous la
forme $`Z_r^{(q)}=Z_r^{(q_1)}B_r`$, où les $`B_r`$ sont indépendants et
$`\mathbb E B_r=q/q_1`$. Le canal $`Y_q`$ est donc une dégradation explicite
de $`Y_{q_1}`$. Par transitivité,

```math
Q_{T,E}\succeq_{\mathrm{ln}}P_{T,q}
\qquad
(0\le q\le q_1).
\qquad\text{(6.3)}
```

Cela explique pourquoi le théorème global porte sur tout
$`p\in[1/2,p_1]`$, et pas seulement sur son extrémité supérieure.

## 7. Tensorisation facteur par facteur

Les observations de triangles montants distincts sont conditionnellement
indépendantes sachant tous les spins, même lorsque les triangles partagent un
sommet. Le théorème 3 de
[Polyanskiy--Wu](https://doi.org/10.4171/MSL/10) tensorise précisément l'ordre
less-noisy sur ce graphe biparti spins--facteurs. Pour tous sommets $`i,j`$,

```math
I(X_i;Y^{(p)}\mid X_j)
\le
I(X_i;E\mid X_j),
\qquad
\frac12\le p\le p_1.
\qquad\text{(7.1)}
```

Ici $`E`$ désigne la collection indépendante des états triangulaires
$`(a_1,s_1,s_1,s_1,e_1)`$. Les deux familles de sorties ne dépendent que de
relations de spins. La symétrie de flip global implique donc
$`I(X_i;Y^{(p)})=I(X_i;E)=0`$, tandis que les spins a priori sont
indépendants. Par conséquent, (7.1) équivaut à

```math
I(X_i;X_j,Y^{(p)})
\le
I(X_i;X_j,E).
\qquad\text{(7.2)}
```

Conditionnellement au sous-graphe révélé par $`E`$, les flips globaux de ses
composantes connexes restent indépendants et uniformes. La relation entre
$`X_i`$ et $`X_j`$ est donc connue exactement s'ils sont connectés, et reste
uniforme sinon. En logarithmes naturels,

```math
I(X_i;X_j,E)
=
(\log 2)\,\mathbb P_E(i\leftrightarrow j).
\qquad\text{(7.3)}
```

Les équations (7.2)--(7.3) sont la borne d'information-percolation adaptée
au canal triangulaire multi-état.

## 8. De l'information mutuelle à l'overlap

Notons $`Y=Y^{(p)}`$ et

```math
m_{ij}(Y)=\mathbb E[X_iX_j\mid Y].
\qquad\text{(8.1)}
```

La symétrie de flip global donne directement
$`\mathbb E[X_i\mid Y,X_j]=X_jm_{ij}(Y)`$. La minoration de Pinsker,
conditionnellement à $`X_j`$, donne donc

```math
\frac12\,\mathbb E[m_{ij}(Y)^2]
\le
I(X_i;Y\mid X_j)
=
I(X_i;X_j,Y).
\qquad\text{(8.2)}
```

En combinant (7.2), (7.3) et (8.2),

```math
\mathbb E[m_{ij}(Y)^2]
\le
2(\log 2)\,\mathbb P_E(i\leftrightarrow j).
\qquad\text{(8.3)}
```

Pour tout estimateur $`\widehat X`$, conditionner par $`Y`$, développer le
carré et utiliser Cauchy--Schwarz donne alors

```math
\mathbb E\left[
\left(\frac1{n_L}\sum_iX_i\widehat X_i\right)^2
\right]
\le
\frac1{n_L}
+\frac{\sqrt{2\log 2}}{n_L^2}
\sum_{i\ne j}
\sqrt{\mathbb P_E(i\leftrightarrow j)}.
\qquad\text{(8.4)}
```

La même preuve traite un estimateur randomisé en conditionnant aussi par son
hasard interne.

## 9. Sous-criticité et double limite sur le tore

Les calculs rationnels (2.1)--(2.3) vérifient toutes les hypothèses du
théorème de phase de
[Chayes--Lei](https://doi.org/10.1007/s10955-005-8078-7) : isotropie,
$`a_1e_1>2s_1^2`$, $`a_1+e_1>2\sqrt2/(3+2\sqrt2)`$ et $`a_1<e_1`$. Il existe
donc $`c,C>0`$ tels que, dans le modèle infini,

```math
\chi(R)
:=
\mathbb P_E(0\leftrightarrow\partial B_R)
\le
Ce^{-cR}.
\qquad\text{(9.1)}
```

Fixons $`R`$, puis prenons $`L`$ assez grand pour que la boule de rayon $`R`$
du tore soit isomorphe à celle de la grille infinie. Il y a au plus
$`C_0n_LR^2`$ paires ordonnées à distance au plus $`R`$. Pour toute autre
paire connectée, la composante du premier sommet atteint nécessairement
$`\partial B_R`$ ; cet événement ne dépend que des triangles rencontrés avant
la première sortie de la boule. Ainsi

```math
\frac1{n_L^2}
\sum_{i\ne j}\sqrt{\mathbb P_E(i\leftrightarrow j)}
\le
\frac{C_0R^2}{n_L}+\sqrt{\chi(R)}.
\qquad\text{(9.2)}
```

L'ordre des limites est essentiel. D'abord $`L\to\infty`$ à $`R`$ fixé,
puis $`R\to\infty`$. Les équations (9.1)--(9.2) montrent que le membre droit
de (8.4) tend vers zéro et prouvent (1.3).

Enfin, si une weak recovery existait avec des constantes
$`\varepsilon,\eta>0`$ au sens du fichier 03, elle imposerait

```math
\liminf_{L\to\infty}
\mathbb E\left[
\left(\frac1{n_L}\sum_iX_i\widehat X_i^{(L)}\right)^2
\right]
\ge
\varepsilon^2\eta>0,
```

en contradiction avec (1.3). Cela achève la preuve du théorème global P809.

## 10. Position dans le programme

La borne $`p=0.809`$ améliore strictement le jalon A0 à $`0.805`$ et se situe
à moins de $`10^{-3}`$ du candidat tangent
$`0.809909289251\ldots`$ du [fichier 11](11_TRIANGLE_BLOCK_SDPI.md). Elle est
entièrement indépendante de la dynamique hiérarchique : c'est désormais le
benchmark rigoureux que cette dynamique devrait dépasser ou expliquer.

## 11. Vérification autonome

Depuis la racine du dépôt :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/rational_a0_less_noisy_certificate.py \
  --candidate p809
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_rational_a0_less_noisy_certificate.py' -v
```

La sortie du certificateur contient

```text
candidate: p809
status: CERTIFIED_PSD
scope: exhaustive
unresolved_regions: 0
variance_gap: 1/50000
binary_margin: 49016699/17274050000
```
