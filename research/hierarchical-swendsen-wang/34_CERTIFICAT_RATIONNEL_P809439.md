# Certificat rationnel et borne globale à $`p=0.809439`$

## 1. Résultat

Fixons

```math
p_2=\frac{809439}{1000000}=0.809439,
\qquad
q_2=2p_2-1=\frac{309439}{500000},
```

et

```math
(a_2,s_2,e_2)
=
\frac1{500000000}(166642280,55571811,166642287),
\qquad
\delta_2=\frac1{50000000}.
\qquad\text{(1.1)}
```

**Théorème local P809439.** Pour tout a priori $`\mu`$ sur les quatre états
relatifs et toute fonction $`f`$,

```math
Q_{E_{a_2,s_2,e_2}}(\mu,f)
-Q_{Y_{q_2}}(\mu,f)
\ge
\delta_2\,\mathrm{Var}_\mu(f).
\qquad\text{(1.2)}
```

**Théorème global P809439.** Sur le tore triangulaire $`\mathbb T_L`$ à
$`n_L=L^2`$ spins i.i.d. uniformes, pour tout
$`p\in[1/2,809439/1000000]`$ et tout estimateur, éventuellement randomisé,
$`\widehat X^{(L)}`$ construit depuis les observations,

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

Il n'y a donc pas de weak recovery au sens du
[fichier 03](03_HIERARCHICAL_WEAK_RECOVERY.md) sur cet intervalle. La borne
rigoureuse devient

```math
p_{\mathrm{WR}}\ge\frac{809439}{1000000}=0.809439.
\qquad\text{(1.4)}
```

## 2. Hypothèses de Chayes--Lei avec marges

La normalisation et la sous-criticité stricte valent

```math
a_2+3s_2+e_2=1,
\qquad
e_2-a_2
=
1-(2a_2+3s_2)
=
\frac7{500000000}>0.
\qquad\text{(2.1)}
```

L'association a la marge exacte

```math
a_2e_2-2s_2^2
=
\frac{10796599147227459}{125000000000000000}>0.
\qquad\text{(2.2)}
```

La condition de densité se vérifie sans approximation de $`\sqrt2`$ :

```math
(a_2+e_2+8)^2-72
=
\frac{777355138600377489}{250000000000000000}>0.
\qquad\text{(2.3)}
```

Comme les membres pertinents sont positifs, (2.3) équivaut à

```math
a_2+e_2
>
\frac{2\sqrt2}{3+2\sqrt2}.
```

Les hypothèses d'isotropie, d'association et de sous-criticité du théorème de
[Chayes--Lei](https://doi.org/10.1007/s10955-005-8078-7) sont donc toutes
strictes.

## 3. Réduction aux quatre polynômes

Le profil exact du canal physique est

```math
c_{q_2}(t)
=
\frac{383009978884t(537257484163-383009978884t)}{H_2(t)},
```

```math
H_2(t)
=
(36313494721+618878000000t)(655191494721-618878000000t).
\qquad\text{(3.1)}
```

Les deux facteurs de $`H_2`$ sont strictement positifs sur $`[0,1]`$. Posons

```math
\bar a_2=a_2-\delta_2=\frac{16664227}{50000000},
\qquad
d_2(t)=\frac{\bar a_2-c_{q_2}(t)}t.
\qquad\text{(3.2)}
```

Le secteur non polarisé est contrôlé par

```math
P_{\mathrm{np}}(t)
=1982400847413554713024317119535
-16886556652293525468619903939898t
+47330451139550895979364000000000t^2
-42569116315311277848000000000000t^3.
\qquad\text{(3.3)}
```

Après extraction du facteur $`1-2t`$ dans
$`d_2(t)-d_2(1/2)`$, le polynôme primitif est

```math
P_{\mathrm{tail}}(t)
=137084007706052200483756974963331002694947
-894938462691318390938821881638944400000000t
+1460628246646072375716946304372000000000000t^2.
\qquad\text{(3.4)}
```

Le numérateur primitif de $`-d_2'(t)`$ est

```math
\begin{aligned}
P_{\mathrm{dec}}(t)
={}&9433172315345871934969498145675028982913299390874787\\
&+303711722682995718660328951892648021559576000000000000t\\
&-1974331516564512256114987698060786064097177024800000000t^2\\
&+2992214636015683935852477362179332800000000000000000000t^3\\
&-364727750098333998369989835623555200000000000000000000t^4.
\end{aligned}
\qquad\text{(3.5)}
```

Enfin, le numérateur primitif de
$`d_2(t)+s_2/[t(1-t)]`$ est

```math
P_{\mathrm{off}}(t)
=5286982498277589323500182269121
-21742080055560441029094634239070t
+27299948009212567865774000000000t^2
-9522669648479236262728000000000t^3.
\qquad\text{(3.6)}
```

## 4. Deux certificats de signe indépendants

Les suites de Sturm exactes donnent :

| polynôme | intervalle | variations | racines | témoin positif |
|---|---:|---:|---:|---:|
| $`P_{\mathrm{np}}`$ | $`[0,1/2]`$ | $`2\to2`$ | $`0`$ | $`P_{\mathrm{np}}(1/4)=107544876270731256409182269121/2`$ |
| $`P_{\mathrm{tail}}`$ | $`[0,1/2]`$ | $`1\to1`$ | $`0`$ | $`P_{\mathrm{tail}}(1/4)=4638657448602126231360648576844902694947`$ |
| $`P_{\mathrm{dec}}`$ | $`[1/2,1]`$ | $`1\to1`$ | $`0`$ | $`P_{\mathrm{dec}}(3/4)=273593896149620185631544922647609396597933222940874787`$ |
| $`P_{\mathrm{off}}`$ | $`[1/2,1]`$ | $`2\to2`$ | $`0`$ | $`P_{\mathrm{off}}(3/4)=638533907674300355677413179637/2`$ |

Aucun polynôme ne s'annule aux extrémités. Les quatre signes sont donc
strictement positifs.

Un second audit exact, indépendant du comptage de Sturm, convertit chaque
polynôme dans la base de Bernstein et subdivise dyadiquement son intervalle
jusqu'à obtenir uniquement des coefficients strictement positifs. Il ferme
avec les nombres de boîtes et profondeurs maximales suivants :

| polynôme | boîtes | profondeur maximale |
|---|---:|---:|
| $`P_{\mathrm{np}}`$ | $`3`$ | $`2`$ |
| $`P_{\mathrm{tail}}`$ | $`12`$ | $`11`$ |
| $`P_{\mathrm{dec}}`$ | $`1`$ | $`0`$ |
| $`P_{\mathrm{off}}`$ | $`11`$ | $`10`$ |

Les tests reconstruisent ces couvertures avec des fractions exactes.

## 5. Dominance diagonale et bords

Dans le secteur $`\max_x\mu_x\le1/2`$, le signe de (3.3) et le lemme des
trois projections du [fichier 11](11_TRIANGLE_BLOCK_SDPI.md) prouvent (1.2).

Dans le secteur polarisé, ordonnons les masses de sorte que $`\mu_0>1/2`$.
Les signes de (3.4)--(3.5) impliquent

```math
d_2(\mu_i)\ge d_2(1/2)\ge d_2(\mu_0)
\qquad
(i=1,2,3).
\qquad\text{(5.1)}
```

Pour $`\{i,j,k\}=\{1,2,3\}`$, posons

```math
D_k=(\mu_0+\mu_k)(1-\mu_0-\mu_k).
```

Comme $`D_k\le\mu_0(1-\mu_0)`$, le signe de (3.6) rend tous les termes hors
diagonale strictement positifs. L'identité exacte est alors

```math
M_{ii}-M_{ij}-M_{ik}
=d_2(\mu_i)-d_2(\mu_0)\ge0.
\qquad\text{(5.2)}
```

La matrice symétrique est à diagonale dominante et positive semi-définie.
Les faces suivent par restriction au support ou par continuité des formes de
canaux finis. La face binaire uniforme conserve la marge

```math
\bar a_2+2s_2-\frac{2q_2^2}{1+q_2^2}
=
\frac{73167609074569033}{43219061840125000000}>0.
\qquad\text{(5.3)}
```

Ceci démontre le théorème local P809439.

## 6. Passage au théorème global

Pour un triangle montant $`T=\{x_0,x_1,x_2\}`$, définissons le quotient

```math
\rho_T(X)
=
(X_{x_0}X_{x_1},X_{x_0}X_{x_2}).
\qquad\text{(6.1)}
```

Le critère $`\chi^2`$ de
[Makur--Polyanskiy](https://doi.org/10.1109/TIT.2018.2839743) transforme
(1.2) en

```math
E_{a_2,s_2,e_2}\succeq_{\mathrm{ln}}Y_{q_2}.
\qquad\text{(6.2)}
```

Toute loi de $`(U,X_{x_0},X_{x_1},X_{x_2})`$ se pousse en avant par
$`\rho_T`$ ; (6.2) se relève donc aux canaux ayant les trois spins pour
entrée. Pour $`0\le q\le q_2`$, le canal $`Y_q`$ est une dégradation de
$`Y_{q_2}`$ par trois BSC indépendants de moyenne $`q/q_2`$.

Les triangles montants partitionnent les arêtes et leurs sorties sont
conditionnellement indépendantes sachant les spins. Le théorème 3 de
[Polyanskiy--Wu](https://doi.org/10.4171/MSL/10) tensorise alors la
comparaison. Pour tous $`i,j`$ et tout $`p\in[1/2,p_2]`$,

```math
I(X_i;X_j,Y^{(p)})
\le
I(X_i;X_j,E)
=
(\log 2)\,\mathbb P_E(i\leftrightarrow j).
\qquad\text{(6.3)}
```

La dernière égalité vient des flips indépendants et uniformes des
composantes du sous-graphe exactement révélé par $`E`$. Avec
$`m_{ij}(Y)=\mathbb E[X_iX_j\mid Y]`$, la symétrie globale et Pinsker donnent

```math
\frac12\,\mathbb E[m_{ij}(Y)^2]
\le
I(X_i;Y\mid X_j)
=
I(X_i;X_j,Y).
\qquad\text{(6.4)}
```

Par conséquent,

```math
\mathbb E\left[
\left(\frac1{n_L}\sum_iX_i\widehat X_i\right)^2
\right]
\le
\frac1{n_L}
+\frac{\sqrt{2\log 2}}{n_L^2}
\sum_{i\ne j}\sqrt{\mathbb P_E(i\leftrightarrow j)}.
\qquad\text{(6.5)}
```

Chayes--Lei et (2.1)--(2.3) donnent, dans le modèle infini,
$`\chi(R):=\mathbb P_E(0\leftrightarrow\partial B_R)\le Ce^{-cR}`$. Pour
$`L` assez grand à $`R`$ fixé, la géométrie locale du tore est celle de la
grille et

```math
\frac1{n_L^2}
\sum_{i\ne j}\sqrt{\mathbb P_E(i\leftrightarrow j)}
\le
\frac{C_0R^2}{n_L}+\sqrt{\chi(R)}.
\qquad\text{(6.6)}
```

Prendre d'abord $`L\to\infty`$, puis $`R\to\infty`$, prouve (1.3). Si une
weak recovery existait avec avantage $`\varepsilon`$ et probabilité
$`\eta`$, le membre gauche de (1.3) aurait une limite inférieure au moins
$`\varepsilon^2\eta`$. Le théorème global est donc démontré.

## 7. Plafond exact du schéma de dominance diagonale

La proximité du plafond n'est pas évaluée par un scan. Fixons

```math
p^+=\frac{809439019}{1000000000}=0.809439019,
\qquad
q^+=\frac{309439019}{500000000},
```

et supposons que la même preuve par séparation de $`d`$, positivité hors
diagonale et sous-criticité stricte puisse fermer avec $`\delta\ge0`$. Posons
$`A=a-\delta`$. Au point rationnel $`t_0=6127079/20000000`$, la condition
$`d(t_0)\ge d(1/2)`$ impose

```math
A\ge
A_0
:=
\frac{c_{q^+}(t_0)-2t_0c_{q^+}(1/2)}{1-2t_0}.
\qquad\text{(7.1)}
```

Au point $`t_1=56557217/100000000`$, combiner
$`d(t_1)+s/[t_1(1-t_1)]>0`$ avec $`2a+3s<1`$ impose

```math
A<
A_{\max}
:=
\frac{1-3(1-t_1)c_{q^+}(t_1)}{3t_1-1}.
\qquad\text{(7.2)}
```

Or le calcul rationnel exact donne

```math
A_0-A_{\max}
=
\frac{
14065222857405417517792295205369626373869821016600100484402322184760321915204457808629397159792953
}{
2601657833694544556860086082446992935025579384744178481259486514689111374965381788339862985625000000000000
}>0.
\qquad\text{(7.3)}
```

Cette incompatibilité persiste pour tout $`p\ge p^+`$. En effet, avec
$`z=q^2`$, le numérateur primitif de $`A_0-A_{\max}`$ est

```math
\begin{aligned}
P(z)={}&-1562500000000000000000000000000000000\\
&-3723156998744256875000000000000000000z\\
&+4474518809248178047275159638775000000z^2\\
&+25105318581993286191938002251843285299z^3\\
&+34353196848684139192702888251025000000z^4\\
&+20387125116385832068338520637956714701z^5\\
&+4504071041702544145298475000000000000z^6.
\end{aligned}
\qquad\text{(7.4)}
```

Sa suite de Sturm a une variation aux deux extrémités de
$`[(q^+)^2,1]`$, donc aucune racine, et $`P((q^+)^2)>0`$. Tous les
dénominateurs éliminés sont positifs. La positivité stricte en $`p^+`$ donne
aussi, par continuité, un voisinage gauche sans certificat. Si
$`p_{\mathrm{DD}}^{\sup}`$ désigne le plafond de cette preuve suffisante,

```math
0.809439
\le
p_{\mathrm{DD}}^{\sup}
<
0.809439019.
\qquad\text{(7.5)}
```

Ce no-go concerne uniquement la séparation univariée suivie de dominance
diagonale ; il ne réfute ni l'ordre less-noisy lui-même ni une autre
certification de la matrice sur le simplexe.

## 8. Vérification autonome

Depuis la racine du dépôt :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/rational_a0_less_noisy_certificate.py \
  --candidate p809439
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_rational_a0_less_noisy_certificate.py' -v
```

La sortie contient

```text
candidate: p809439
status: CERTIFIED_PSD
scope: exhaustive
unresolved_regions: 0
variance_gap: 1/50000000
binary_margin: 73167609074569033/43219061840125000000
order_slack: 7/500000000
self_dual_slack: 7/500000000
```
