# Couper la géante au niveau critique : audit

**Verdict : la réduction à la géante est exacte. La coupe seule ne donne
pas de meilleur seuil.** Deux raccourcis sont réfutés ci-dessous ; une cible
précise reste ouverte pour dépasser la borne d'impossibilité à 0,809439.

## 1. Pourquoi cette coupe est naturelle

On garde le modèle de la [note précédente](03_RECOVERY.md), sur des tores
triangulaires de $`n=L^2`$ sommets, avec $`2p-1>q_c`$. Le niveau $`c=t_c(p)`$ est celui du
gel critique, pas une nouvelle température de la postérieure.
À $`p=0{,}81`$ :

| Quantité | Valeur |
|---|---:|
| densité de liens à la coupe critique | 0,347296 |
| densité finale Swendsen–Wang | 0,62 |
| niveau critique $`c`$ | 0,386168 |
| poids résiduel par arête $`h=(1-c)u_p`$ | 0,890063 |

La géante est celle de la **coupe finale à un**. À la coupe critique, la
fraction du plus grand amas tend vers zéro.

Plus précisément, $`c`$ est la plus grande coupe **fixe** pour laquelle
$`S_n(t)=n^{-2}\mathbb E\sum_{C\in\Pi_t}|C|^2\to0`$.
Au-dessus, une géante apparaît et ce terme ne disparaît plus.
La coupe critique regroupe donc autant que possible sans laisser un amas
macroscopique dans cette borne géométrique. Cela ne prouve pas son optimalité
pour l'information ; des coupes légèrement sous-critiques peuvent être plus
faciles à contrôler. Une suite déterministe $`c_L\to t_c`$ conserve aussi
$`S_n(c_L)\to0`$, même en approchant par au-dessus : la fenêtre proche de
la criticité reste donc à explorer. Les lois géométriques sont ici
moyennées sur les observations. Les résultats de percolation utilisés sont ceux de
[Grimmett–Manolescu](https://arxiv.org/abs/1105.5535) et, pour l'unicité de
la géante sur les graphes finis, d'[Easo–Hutchcroft](https://arxiv.org/abs/2112.12778).

## 2. Garder la géante change la loi

Notons $`A=A_c`$ les arêtes précoces, $`\mathcal R=\Pi_1`$ la partition
finale et $`G`$ sa plus grande composante, avec départage fixé par les sommets.
On conserve $`H=(A,\mathcal R)`$, **sans conserver les arêtes finales**.
Les garder imposerait toutes les parités à l'intérieur de $`G`$ : recouper
ne libérerait aucun spin relatif.

La conditionnelle exacte, avec $`s_e(\sigma)`$ l'indicatrice d'arête satisfaite,
est

```math
\pi_O(\sigma\mid H)\propto
\mathbf1_{\{A\text{ satisfait}\}}
\sum_{\substack{B\subseteq E\setminus A\\\Pi(A\cup B)=\mathcal R}}
\prod_{e\in B}\left[(e^h-1)s_e(\sigma)\right].
```

C'est la somme sur toutes les façons de connecter les amas précoces pour
obtenir exactement les racines finales. Elle remplace le simple Gibbs
résiduel lorsque la géante finale est connue.

Ce conditionnement sert ici à la preuve. Une MCMC qui conserve aussi
$`\mathcal R`$ doit employer cette loi ; ses mises à jour à $`c=0`$ ne sont
plus celles de Glauber décrites dans la [note 1](01_HIERARCHIE.md).

**Preuve.** Multiplier les probabilités des trois états d'une arête par
son poids postérieur donne respectivement $`(e^{u_p}-e^h)s_e`$ pour une
arête précoce, $`(e^h-1)s_e`$ pour une arête tardive, et $`1`$ pour une
arête toujours fermée. Fixer $`A`$, puis sommer les arêtes tardives $`B`$,
donne la formule. La somme se factorise entre les racines de $`\mathcal R`$.

## 3. La réduction exacte au cœur de la géante

Posons $`m^H_{ij}=\mathbb E[\Sigma_i\Sigma_j\mid O,H]`$ et

```math
d_{ij}(O)=\mathbb E_{H\mid O}
\left[\mathbf1_{\{i,j\in G,\ C_A(i)\ne C_A(j)\}}m^H_{ij}\right],
\qquad
F_n=\frac1{n^2}\sum_{i,j}\mathbb E[d_{ij}(O)^2].
```

Pour le $`Q_n`$ de la note précédente, on a exactement

```math
|Q_n-F_n|\le2(S_n^c+S_n^{\mathrm{hors}}),\qquad
S_n^c=\frac1{n^2}\mathbb E\sum_{C\in\Pi_c}|C|^2,\qquad
S_n^{\mathrm{hors}}=\frac1{n^2}\mathbb E\sum_{R\in\mathcal R,\ R\ne G}|R|^2.
```

**Preuve.** Les orientations globales des racines sont indépendantes et
uniformes, donc $`m^H_{ij}=0`$ entre racines distinctes. Les paires restantes
non comptées dans $`d`$ sont dans un même amas précoce ou une petite racine
finale. Leur masse est majorée par les deux sommes ci-dessus ; appliquer
$`|x^2-y^2|\le2|x-y|`$ pour $`|x|,|y|\le1`$ conclut.

À $`c=t_c`$, le premier terme est contrôlé par la fraction du plus grand
amas critique ; le second par celle de la deuxième racine finale.
Ils tendent donc tous deux vers zéro dans le régime annoncé.
Ainsi **prouver $`F_n(0{,}81)\to0`$ suffirait à améliorer le seuil**.
La moyenne sur $`H`$ doit précéder le carré. Deux copies de $`H`$ utilisées
pour ce carré sont indépendantes conditionnellement à $`O`$.

## 4. Pourquoi la percolation arête par arête échoue

Révéler $`A_t`$ donne une parité exacte sur les arêtes ouvertes. Sur une
arête fermée, la fiabilité devient
$`r_t=(2p-1-q)/(1-q)`$, où $`q=q_p(t)`$.
La comparaison d'[information-percolation, théorème 3.6](https://arxiv.org/abs/1806.03227)
appliquée après ce gel donne

```math
\eta_{\mathrm{gel}}=q+(1-q)r_t^2
=(2p-1)^2+\frac{4q(1-p)^2}{1-q}.
```

**C'est une majoration moins bonne qu'avant le gel.** À la coupe critique,
$`\eta_{\mathrm{gel}}>q_c`$ dès que $`c<1`$.
À $`p=0{,}81`$, elle vaut $`0{,}461234`$, contre $`0{,}3844`$ avant gel.
La seule contraction des amas suivie de cette comparaison ne peut donc pas
améliorer le seuil. Ce constat ne réfute pas une analyse conservant la
moyenne signée qui définit $`F_n`$.

## 5. Où un gain peut encore se trouver

Pour une racine composée de deux amas précoces, fixons leurs configurations
internes. Si $`k_+`$ et $`k_-`$ comptent les arêtes entre amas favorisant
chacune des deux orientations relatives, leurs poids exacts sont

```math
\mathbb P(z_1z_2=\pm1\mid H,O)\propto e^{h k_\pm}-1.
```

Le terme $`-1`$ impose une connexion finale. Deux liens opposés de même
poids donnent une corrélation nulle ; des liens tous de même signe donnent
une corrélation de module un. Ces deux situations existent déjà dans des
triangles, selon leur frustration. Il n'y a donc **aucune contraction
stricte valable pour toutes les fusions**.

La piste utile consiste à contrôler statistiquement ces compensations,
avec les connexions extérieures conservées, puis à établir $`F_n\to0`$.
Le calcul d'un triangle isolé ne suffit pas. **Aucune nouvelle borne de
weak recovery n'est démontrée par cet audit.**

Les identités finies et les contre-exemples sont vérifiés par
[`check_giant_cut.py`](check_giant_cut.py), par énumération exhaustive avec
poids flottants ; aucune extrapolation numérique en volume n'est utilisée.
