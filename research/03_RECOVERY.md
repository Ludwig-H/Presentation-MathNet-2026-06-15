# La coupe critique : ce qu'elle apporte, ce qu'il reste à prouver

## Quelle criticité ?

Sur la grille triangulaire, considérons le cas binaire symétrique du GSBM :
$`Y_e=\Sigma_i\Sigma_j Z_e`$, avec bruits indépendants et
$`\mathbb P(Z_e=1)=p\in(1/2,1)`$.
Alors $`W_e=u_pY_e`$, où $`u_p=\log(p/(1-p))`$.
Depuis la vérité — ou une réplique postérieure après moyenne sur les
observations — les liens gelés sont indépendants, de probabilité

```math
q_p(t)=p(1-e^{-u_pt}),\qquad
t_c(p)=\frac{-\log(1-q_c/p)}{u_p},\qquad
q_c=2\sin(\pi/18).
```

La valeur $`q_c`$ est le seuil de percolation par arêtes du réseau triangulaire
([Grimmett, §3](https://www.statslab.cam.ac.uk/~grg/papers/UScornell5.pdf)).
La coupe critique appartient à $`[0,1]`$ si et seulement si
$`p\ge(1+q_c)/2\simeq0{,}673648`$. **À observation fixée, les liens ne sont
pas une percolation indépendante de paramètre $`q_p(t)`$.**

Cette coupe donne une organisation géométrique. Elle n'est ni une température
critique de la postérieure, ni un seuil de récupération démontré.

## La quantité à faire tendre vers zéro

Pour $`n`$ sommets, posons

```math
m_{ij}(O)=\mathbb E[\Sigma_i\Sigma_j\mid O],\qquad
Q_n=\frac1{n^2}\sum_{i,j}\mathbb E[m_{ij}(O)^2].
```

Pour tout estimateur $`\tau(O)`$ à valeurs dans $`\{-1,1\}^n`$,
deux applications de Cauchy–Schwarz donnent

```math
\mathbb E\left[\left(\frac1n\sum_i\Sigma_i\tau_i\right)^2\right]
\le\frac1{n^2}\sum_{i,j}\mathbb E|m_{ij}(O)|\le\sqrt{Q_n},
\qquad
\mathbb E[\mathrm{ov}_n]\le\frac12+\frac12Q_n^{1/4}.
```

Ainsi $`Q_n\to0`$ interdit la weak recovery. Une configuration obtenue après
un pas invariant depuis la vérité a la bonne marginale postérieure,
mais n'est pas une réplique indépendante de cette vérité.

## Ce que la coupe permet de séparer exactement

Supposons $`p\ge(1+q_c)/2`$, pour que la coupe soit dans l'horizon.
Tirer $`A=A_{t_c}`$ selon sa loi jointe avec les spins. Notons
$`m_{ij}(O,A)=\mathbb E[\Sigma_i\Sigma_j\mid O,A]`$ et

```math
d_{ij}(O)=\mathbb E_{A\mid O}
\left[\mathbf1_{\{C_A(i)\ne C_A(j)\}}m_{ij}(O,A)\right],
\qquad
S_n=\frac1{n^2}\mathbb E\sum_{C\in\Pi_{t_c}}|C|^2.
```

On dispose de la borne finie

```math
\left|Q_n-\frac1{n^2}\sum_{i,j}\mathbb E[d_{ij}(O)^2]\right|\le2S_n.
```

**Preuve.** Écrire $`m_{ij}(O)=b_{ij}(O)+d_{ij}(O)`$, où $`b`$ est la
contribution des paires dans un même amas. Alors
$`|b_{ij}|\le\mathbb P(C_A(i)=C_A(j)\mid O)`$ et
$`|m_{ij}^2-d_{ij}^2|\le2|b_{ij}|`$. Sommer donne la borne.

Sur des tores triangulaires croissants à criticité, $`S_n\to0`$ : la fraction
du plus grand amas tend vers zéro, et la somme des carrés est majorée par
cette fraction. Cela découle de l'absence d'amas infini critique, en
contrôlant d'abord la probabilité d'atteindre un rayon fixé.

**Le problème restant est donc précis : contrôler les corrélations entre
amas distincts, sous leur Gibbs résiduel signé.** La moyenne en $`A`$ doit
précéder le carré. Inverser cet ordre perd les compensations de signes.

## La prochaine étape utile

Écrire ces corrélations sur les triangles et leurs bords, conserver les
interactions extérieures, puis chercher une contraction uniforme quand le
volume grandit. Tester plusieurs coupes : rien ne démontre encore que
$`t_c`$ donne la meilleure borne.

La référence élémentaire à retrouver est l'information-percolation :
$`\eta=(2p-1)^2<q_c`$ interdit la recovery, soit
$`p<(1+\sqrt{q_c})/2\simeq0{,}794659`$
([Abbe–Boix, théorème 3.6](https://arxiv.org/abs/1806.03227)).
Le certificat historique à 0,809439 fixe un objectif plus exigeant.
**Aucune amélioration de ce seuil n'est démontrée ici.**

L'[audit de la coupe dans la géante](04_GEANTE_CRITIQUE.md) précise comment
conserver la géante finale sans fausser cette conditionnelle, et démontre
l'échec de la majoration arête par arête après gel critique.
