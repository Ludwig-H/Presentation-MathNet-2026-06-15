# Projections de heat bath et comparaison favorable à $`p=4/5`$

Cette note poursuit le programme du transfert répliqué avec un objectif
précis : séparer ce qui est désormais démontré du travail encore nécessaire
pour obtenir une impossibilité de weak recovery à

```math
p_0=\frac45=0.8.
```

Le bilan est le suivant.

1. Si $`i,j`$ sont dans deux racines distinctes à la coupe $`1`$, tout sweep
   complet top-down ou bottom-up donne exactement $`H_S(i,j)=0`$.
2. À dendrogramme fixé, un sweep est un produit de projections orthogonales.
   Son second moment répliqué est la norme $`L^2`$ de ce produit, ou encore
   l'autocorrélation d'un sweep palindromique.
3. À coupe et message ancestral fixés, le nœud critique est plus persistant
   que le nœud tardif sauf sur un événement d'**anti-alignement** explicite.
4. La monotonie point par point sans cette exception est fausse, y compris à
   $`p=0.8`$ pour un squelette admissible.
5. Sous la loi complète du compte et à taille fixée, le bucket critique
   domine pourtant tout bucket tardif au sens de Blackwell.
6. Sur les petits tores, la moyenne du second moment est également plus
   grande pour la fenêtre critique que pour les fusions tardives. Cela
   soutient HF-S2, sans la prouver.
7. Le verrou restant pour $`p=0.8`$ se décompose en deux lemmes quantitatifs :
   domination critique/postcritique et accumulation de blocs critiques
   contractants.

## 1. Un sweep est un produit de projections

Fixons $`O,D`$ et écrivons

```math
\pi_D(d\sigma)=\nu_O(d\sigma\mid D).
```

Pour une mise à jour $`u`$, soit $`\mathcal G_u`$ la tribu engendrée par les
variables non rééchantillonnées. Son opérateur de heat bath est

```math
(P_ug)(\sigma)
=
\mathbb E_{\pi_D}[g(\widetilde\sigma)\mid\mathcal G_u](\sigma).
```

Dans $`L^2(\pi_D)`$,

```math
P_u^2=P_u,
\qquad
P_u^*=P_u,
\qquad
\|P_ug\|_2\le\|g\|_2.
```

Ainsi chaque heat bath est une projection orthogonale. Pour un programme

```math
S=(u_1,\ldots,u_M),
```

l'opérateur de Markov du sweep vaut, avec la convention d'action sur les
fonctions,

```math
K_{S,D}=P_{u_1}\cdots P_{u_M}.
```

Cette représentation est standard pour le Gibbs sampler. Ici elle est
appliquée aux blocs d'orientations définis par le dendrogramme, et non aux
seuls spins individuels.

### Théorème 1.1 — identité $`L^2`$ du transfert répliqué, statut : établi

Pour

```math
f_{ij}(\sigma)=\sigma_i\sigma_j,
```

on a, conditionnellement à $`O,D`$,

```math
\boxed{
\mathbb E_{\sigma\sim\pi_D}
[H_S(i,j)^2]
=
\|K_{S,D}f_{ij}\|_{L^2(\pi_D)}^2
=
\langle
f_{ij},K_{S,D}^*K_{S,D}f_{ij}
\rangle_{\pi_D}.
}
\tag{1.1}
```

#### Preuve

À état initial fixé,

```math
H_S(i,j)
=
f_{ij}(\sigma)(K_{S,D}f_{ij})(\sigma).
```

Comme $`f_{ij}^2=1`$, son carré est
$`(K_{S,D}f_{ij})^2`$. L'intégration sous $`\pi_D`$ donne la première
égalité ; l'identité d'adjoint donne la seconde.

Le produit

```math
K_{S,D}^*K_{S,D}
=
P_{u_M}\cdots P_{u_1}P_{u_1}\cdots P_{u_M}
```

est le sweep palindromique. Il est positif même si le sweep systématique
$`K_{S,D}`$ n'est pas auto-adjoint. Le transfert tensoriel à deux répliques
du fichier 18 et cette formule calculent donc le même second moment.

### Corollaire 1.2 — dissipation pythagoricienne, statut : établi

Posons $`g_0=f_{ij}`$ et appliquons les projections dans l'ordre dans lequel
elles agissent sur la fonction :

```math
g_r=P_{u_{M-r+1}}g_{r-1},
\qquad 1\le r\le M.
```

Alors

```math
\boxed{
1-\|K_{S,D}f_{ij}\|_2^2
=
\sum_{r=1}^M
\|(I-P_{u_{M-r+1}})g_{r-1}\|_2^2.
}
\tag{1.2}
```

Chaque terme est une variance conditionnelle. Pour obtenir $`p=0.8`$, il
faut donc montrer que le chemin de la paire dissipe asymptotiquement toute
la norme, même après la modification des fonctions $`g_r`$ par les mises à
jour précédentes.

## 2. Racines distinctes : annulation exacte

Ajoutons à la fin d'un programme quelconque la recoloration globale,
indépendante et uniforme, de chaque racine finale de $`D`$. Cette opération
est un heat bath valide sous l'a priori uniforme.

### Lemme 2.1 — effacement des racines, statut : établi

Pour le programme ainsi complété,

```math
\boxed{
\beta_{ij}>1
\quad\Longrightarrow\quad
H_S(i,j)=0.
}
\tag{2.1}
```

Si $`i,j`$ sont dans une même racine, ces recolorations ne changent pas leur
relation et ne modifient donc pas $`H_S(i,j)`$.

#### Preuve

Lorsque les racines sont distinctes, le bit global de la racine contenant
$`i`$ est équitable et indépendant de celui de la racine contenant $`j`$.
Il retourne exactement un des deux sommets avec probabilité $`1/2`$ ; la
parité finale a donc une espérance nulle. Dans une même racine, le flip
global retourne simultanément les deux sommets et laisse $`f_{ij}`$ invariant.

Ce lemme formalise exactement le cas évoqué dans la question : une paire qui
n'est même pas connectée à $`\beta=1`$ est éliminée sans approximation. Il
autorise à concentrer tout le travail sur les paires du même arbre.

### Corollaire 2.2 — sweep complet dans les deux ordres, statut : établi

Le même résultat vaut sans ajouter une étape au programme : tout sweep
top-down ou bottom-up qui met à jour chaque racine finale contient déjà sa
recoloration globale équitable. En particulier,

```math
\boxed{
\beta_{ij}>1
\quad\Longrightarrow\quad
H_{S_{\rm TD}}(i,j)=H_{S_{\rm BU}}(i,j)=0.
}
\tag{2.2}
```

En effet, sous l'a priori uniforme, la densité conditionnelle se factorise
sur les racines finales. Dans une racine $`R`$, le heat bath de son dernier
nœud contient la projection

```math
\mathcal R_R g(\sigma)
=\frac12\{g(\sigma)+g(\sigma^R)\}.
```

Tous les heat baths descendants sont équivariants par le flip global
$`\sigma\mapsto\sigma^R`$ ; ils préservent donc les secteurs pair et impair
de cette involution. Si $`i\in R`$ et $`j\notin R`$, alors $`f_{ij}`$ est
impair dans $`R`$ et $`\mathcal R_R f_{ij}=0`$. En bottom-up, la projection
de racine tue directement $`f_{ij}`$ à la fin de son action sur les
fonctions ; en top-down, les opérateurs descendants préservent d'abord son
caractère impair, puis la même projection le tue. Cette preuve ne suppose
pas l'indépendance des décisions internes aux deux arbres.

## 3. Les quatre classes de paires

Fixons une fenêtre gauche en coordonnée de percolation,

```math
[q_\triangle-\delta_q,q_\triangle].
```

Pour une paire lointaine, les quatre classes sont :

```math
\begin{array}{c|c}
\text{classe}&q_p(\beta_{ij})\cr
\hline
\text{précoce}&<q_\triangle-\delta_q\cr
\text{critique}&\in[q_\triangle-\delta_q,q_\triangle]\cr
\text{tardive}&\in(q_\triangle,2p-1]\cr
\text{racines distinctes}&>2p-1.
\end{array}
\tag{3.1}
```

La première classe a une masse qui disparaît lorsque la distance diverge et
que la fenêtre est ensuite resserrée, par décroissance sous-critique et RSW.
La dernière a un transfert nul par le lemme 2.1. Le problème global est donc
réduit à

```math
\boxed{
\text{comparer la classe tardive à la classe critique,
puis contracter la classe critique.}
}
\tag{3.2}
```

Définissons

```math
A_{L,S}^{\mathrm c}(p)
=
\mathbb E[H_S(I_L,J_L)^2\mid\text{classe critique}],
```

et $`A_{L,S}^{\mathrm{late}}(p)`$ de même. La forme minimale de HF-S2 est

```math
A_{L,S}^{\mathrm{late}}(p)
\le
A_{L,S}^{\mathrm c}(p)+\varepsilon_{L,S}^{\mathrm{HF}}(p),
\qquad
\varepsilon_{L,S}^{\mathrm{HF}}(p)\longrightarrow0.
\tag{3.3}
```

## 4. Comparaison exacte d'un nœud critique et tardif

Fixons un bucket de taille $`m`$ contenant $`K`$ arêtes satisfaites. Posons

```math
d=2K-m
```

et, lorsque $`1\le K\le m-1`$,

```math
\ell_t
=
\log\frac K{m-K}
+(1-t)u_p d.
\tag{4.1}
```

Soit $`B`$ le message extérieur complet. La persistance locale est

```math
\eta_t(B,K,m)
=
\tanh^2\left(\frac{B+\ell_t}{2}\right).
\tag{4.2}
```

Pour $`d\ne0`$, écrivons

```math
c_t=\mathrm{sgn}(d)\ell_t>0,
\qquad
b=\mathrm{sgn}(d)B.
```

Si $`\beta_c\le t\le1`$, alors $`c_{\beta_c}\ge c_t`$.

### Lemme 4.1 — critère exact d'anti-alignement, statut : établi

Pour $`1\le K\le m-1`$ et $`d\ne0`$,

```math
\boxed{
\eta_t(B,K,m)>\eta_{\beta_c}(B,K,m)
\quad\Longleftrightarrow\quad
b<-\frac{c_{\beta_c}+c_t}{2}.
}
\tag{4.3}
```

En dehors de cet événement, la fusion critique est bien la plus persistante.
En particulier, $`dB\ge0`$ suffit.

De plus, avec $`c_{\mathrm{rel}}=2/(3\sqrt3)`$,

```math
\boxed{
[\eta_t-\eta_{\beta_c}]_+
\le
c_{\mathrm{rel}}u_p|d|(t-\beta_c)
\mathbf1_{\mathcal A_t},
}
\tag{4.4}
```

où $`\mathcal A_t`$ est l'événement de droite dans (4.3).

#### Preuve

La fonction $`x\mapsto\tanh^2(x/2)`$ est croissante en $`|x|`$. Or

```math
|b+c_t|>|b+c_{\beta_c}|
```

équivaut, après différence des carrés, à
$`b<-(c_{\beta_c}+c_t)/2`$. La borne (4.4) vient de la constante de
Lipschitz exacte et de

```math
c_{\beta_c}-c_t=u_p|d|(t-\beta_c).
```

Les cas $`K=0,m`$ sont déterministes, et $`d=0`$ ne dépend pas de $`t`$.

### Conséquence moyenne

Sous un couplage conservant $`(m,K,B)`$ entre les expériences critique et
tardive,

```math
\mathbb E\eta_t
\le
\mathbb E\eta_{\beta_c}
+c_{\mathrm{rel}}u_p
\mathbb E\left[
|2K-m|(t-\beta_c)\mathbf1_{\mathcal A_t}
\right].
\tag{4.5}
```

La comparaison favorable locale est donc ramenée à une queue d'anti-alignement.
Pour la grille, il faut encore payer l'erreur de couplage des géométries
$`(m,K,B)`$. Pour le sweep complet, $`B`$ devient un état-frontière dynamique
et (4.5) doit être appliquée par blocs, pas multipliée comme si les nœuds
étaient indépendants.

### Théorème 4.2 — ordre de Blackwell des buckets, statut : établi

Le contre-exemple pointwise ne survit pas lorsque le compte $`K`$ est moyenné
avec sa loi bayésienne correcte. Posons $`s=s_p(t)`$ et introduisons l'état
latent $`X\in\{+,-\}`$ de la parité vraie. Pour un bucket non marqué de taille
$`m`$, l'expérience binaire observée est

```math
\begin{aligned}
P_s^+(K=k)
&=\binom{m-1}{k-1}s^{k-1}(1-s)^{m-k},\\
P_s^-(K=k)
&=\binom{m-1}{k}s^{m-1-k}(1-s)^k.
\end{aligned}
\tag{4.6}
```

Ici $`0\le k\le m`$ et un coefficient binomial hors support vaut zéro.

Autrement dit,

```math
K\mid X=+\ \overset{\mathrm d}=\ 1+\mathrm{Bin}(m-1,s),
\qquad
K\mid X=-\ \overset{\mathrm d}=\ \mathrm{Bin}(m-1,1-s).
\tag{4.7}
```

Si $`1/2\le s_2\le s_1\le1`$, alors l'expérience
$`\mathcal E_{m,s_1}`$ domine $`\mathcal E_{m,s_2}`$ au sens de Blackwell.
Il existe donc un noyau stochastique $`G_{s_1\to s_2}`$, indépendant de
$`X`$, tel que

```math
P_{s_2}^x=P_{s_1}^xG_{s_1\to s_2},
\qquad x\in\{+,-\}.
\tag{4.8}
```

#### Preuve

Le rapport de vraisemblance vaut, pour $`1\le k\le m-1`$,

```math
\frac{P_s^+(k)}{P_s^-(k)}
=
\frac{k}{m-k}
\left(\frac{s}{1-s}\right)^{2k-m};
\tag{4.9}
```

il est strictement croissant en $`k`$ lorsque $`s\ge1/2`$. Les tests de
Neyman--Pearson sont donc les seuils $`K\ge r`$. À seuil entier fixé,

```math
\mathrm{TPR}_s(r)
=\mathbb P\{1+\mathrm{Bin}(m-1,s)\ge r\}
```

croît avec $`s`$, tandis que

```math
\mathrm{FPR}_s(r)
=\mathbb P\{\mathrm{Bin}(m-1,1-s)\ge r\}
```

décroît. Chaque sommet de la courbe ROC à $`s_1`$ est donc au nord-ouest du
sommet correspondant à $`s_2`$. Après randomisation aux seuils, la concavité
des courbes ROC donne leur domination en tout point. Pour deux états, le
théorème de comparaison de Blackwell identifie cette domination à la
dégradation (4.8).

Puisque $`s_p(t)`$ décroît avec $`t`$, on obtient pour tout
$`\beta_c\le t\le1`$ :

```math
\boxed{
\mathcal E_{m,s_p(\beta_c)}
\succeq_{\rm Blackwell}
\mathcal E_{m,s_p(t)}.
}
\tag{4.10}
```

Cette conclusion reste vraie en présence d'une information latérale
$`B`$ arbitraire dès que $`B\perp K\mid X`$. En particulier, pour toute loi
de $`X\mid B`$,

```math
\mathbb E\!\left[
\mathbb E[X\mid B,K_{\beta_c}]^2
\right]
\ge
\mathbb E\!\left[
\mathbb E[X\mid B,K_t]^2
\right].
\tag{4.11}
```

Conditionnellement à un squelette non marqué fixé, les buckets disjoints ont
des marques indépendantes. Ainsi (4.11) enlève complètement la queue
d'anti-alignement pour un update top-down dont l'état latéral n'a utilisé que
les buckets ancestraux. Elle ne compare pas encore deux géométries de
Kruskal différentes et ne suffit pas à composer le sweep : après l'update,
les descendants utilisent le bucket courant dans leurs propres messages
ancestraux.

La restriction à taille fixée est essentielle : le contre-lemme 4.2 du
fichier 20 certifie qu'à $`p=t=4/5`$ un bucket critique $`m=4`$ et un bucket
tardif $`m=2`$ sont incomparables au sens de Blackwell. « Avancer le temps »
ne permet donc pas de remplacer librement la géométrie de l'interface.

Le script `favorable_time_comparison.py` contre-audite (4.10) par le critère
équivalent de convexité des croyances postérieures, sans solveur externe.

## 5. Contre-exemple point par point à $`p=0.8`$

La monotonie sans événement exceptionnel est fausse dans le modèle exact.
Prenons

```math
p=0.8,
\qquad
u_p=\log4,
\qquad
m=3,
\qquad
K=2.
```

Ajoutons un ancêtre au niveau $`0.81`$ dont les trois groupes ont les tailles
pondérées

```math
u_p(1,6,6)
```

et les comptes satisfaits

```math
u_p(1,0,6).
```

Ces comptes sont admissibles et donnent

```math
B=-1.509523361867\ldots
```

Le niveau critique vaut

```math
\beta_c(0.8)=0.410716539196\ldots
```

et les messages locaux sont

```math
\ell_{\beta_c}=1.510067519374\ldots,
\qquad
\ell_{0.8}=0.970406052784\ldots
```

Par conséquent,

```math
\eta_{\beta_c}=7.4026844\,10^{-8},
\qquad
\eta_{0.8}=0.069281670285\ldots
\tag{5.1}
```

La fusion tardive est ici plus persistante parce que le renforcement local
critique annule presque exactement un ancêtre opposé. Le script
[`favorable_time_comparison.py`](computations/favorable_time_comparison.py)
reproduit ces valeurs. Ce contre-exemple ne réfute pas HF-S2 **en moyenne** ;
il réfute seulement toute preuve pathwise ignorant la loi du message
ancestral. Il ne contredit pas non plus le théorème 4.2 : celui-ci moyenne
$`K`$ sous les deux lois de (4.6), alors que (5.1) fige le même compte
anti-aligné dans les deux expériences.

## 6. Premier bloc quantitatif à $`p=0.8`$

Un bucket critique de taille deux a

```math
K=1+\mathrm{Bernoulli}(s_c),
\qquad
s_c(0.8)=0.693582222752\ldots
```

Si $`K=2`$, la parité est forcée. Si $`K=1`$, le message local est nul et la
persistance vaut $`\tanh^2(B/2)`$. Sur l'événement de screening
$`|B|\le b`$, le coefficient moyen exact est donc

```math
\boxed{
\kappa_2(b)
=
s_c+(1-s_c)\tanh^2(b/2).
}
\tag{6.1}
```

À $`p=0.8`$ :

| $`b`$ | $`\kappa_2(b)`$ |
|---:|---:|
| $`0`$ | $`0.693582222752`$ |
| $`0.5`$ | $`0.711962739449`$ |
| $`1`$ | $`0.759018433743`$ |
| $`1.5`$ | $`0.817195502689`$ |
| $`2`$ | $`0.871312395742`$ |

### Corollaire 6.1 — dégradation explicite pour $`m=2`$

L'expérience de compte à deux arêtes est un canal d'effacement : sous
$`X=+`$, les sorties $`K=2`$ et $`K=1`$ ont les masses $`s`$ et $`1-s`$ ;
sous $`X=-`$, les sorties $`K=0`$ et $`K=1`$ ont les mêmes masses. Pour
$`s_t\le s_c`$, le noyau de Blackwell de (4.8) est explicite :

```math
0\longmapsto
\begin{cases}
0,&s_t/s_c,\\
1,&1-s_t/s_c,
\end{cases}
\qquad
2\longmapsto
\begin{cases}
2,&s_t/s_c,\\
1,&1-s_t/s_c,
\end{cases}
\qquad
1\longmapsto1.
\tag{6.2}
```

À $`p=0.8`$, entre $`\beta_c`$ et $`t=1`$, la probabilité d'effacement
supplémentaire vaut exactement numériquement

```math
1-\frac{1/2}{s_c}
=0.279104937240\ldots.
\tag{6.3}
```

Ainsi la partie « comparaison critique/tardive » d'un bloc $`m=2`$ est
entièrement résolue. Ce noyau n'efface pas le verrou de composition : un
descendant voit aussi ce compte comme message ancestral.

Cela donne un objectif vérifiable. Si l'on peut extraire sous Palm critique
un nombre $`N_L^{(b)}\to\infty`$ de blocs de taille deux :

1. dont le message-frontière est borné par $`b`$ dans les deux répliques ;
2. dont le noyau tardif est une dégradation BSC du noyau critique ;
3. dont les erreurs de découplage ont une somme $`o(1)`$ ;

alors

```math
\mathbb E[H_S(I_L,J_L)^2\mid\text{critique}]
\le
\mathbb E[\kappa_2(b)^{N_L^{(b)}}]+o(1)
\longrightarrow0.
\tag{6.4}
```

L'indépendance des marques de buckets disjoints sachant le squelette aide
pour (6.1), mais elle ne donne pas seule les points 1--3 : le sweep modifie
l'état-frontière et les ancêtres sont partagés.

## 7. Diagnostic fini de HF-S2 à $`p=0.8`$

Le script
[`pair_favorability_diagnostic.py`](computations/pair_favorability_diagnostic.py)
classe exactement toutes les paires lointaines, choisit une paire uniforme
dans chaque classe non vide, puis pondère l'environnement par le nombre exact
de paires de cette classe. Il estime le carré conditionnel sans biais de
plug-in. La fenêtre utilisée est
$`q_p(\beta_{ij})\in[q_\triangle-0.05,q_\triangle]`$ et la distance est au
moins $`L/4`$.

Pour un sweep :

| $`L`$ | répétitions | critique BU | tardif BU | $`\Delta`$ BU | critique TD | tardif TD | $`\Delta`$ TD |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 200 | $`0.807\pm0.028`$ | $`0.543\pm0.036`$ | $`0.264\pm0.044`$ | $`0.793\pm0.032`$ | $`0.559\pm0.036`$ | $`0.234\pm0.046`$ |
| 6 | 120 | $`0.799\pm0.032`$ | $`0.606\pm0.034`$ | $`0.192\pm0.043`$ | $`0.812\pm0.031`$ | $`0.602\pm0.036`$ | $`0.210\pm0.043`$ |
| 8 | 60 | $`0.766\pm0.048`$ | $`0.688\pm0.054`$ | $`0.077\pm0.071`$ | $`0.766\pm0.050`$ | $`0.703\pm0.050`$ | $`0.063\pm0.068`$ |

Les trois lignes se reproduisent depuis la racine du dépôt avec :

```bash
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 4 --repetitions 200 --sweeps 200 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 6 --repetitions 120 --sweeps 160 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 8 --repetitions 60 --sweeps 120 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
```

Les erreurs marginales et celles de
$`\Delta=A^{\mathrm c}-A^{\mathrm{late}}`$ sont des jackknifes
delete-one-environnement du ratio pondéré ; le contraste conserve
l'appariement des environnements. Les six contrastes ponctuels sont positifs.
À $`L=8`$, ils ne sont toutefois qu'à environ une erreur standard de zéro :
ce volume ne fournit donc pas une évidence séparée forte. La classe
« racines distinctes » donne un second moment compatible avec zéro ; sa
valeur théorique est exactement zéro par le corollaire 2.2.

Ce diagnostic soutient le sens de HF-S2, mais il révèle aussi le second
verrou : le second moment critique ne décroît pas encore avec $`L`$ sur ces
tailles. La comparaison favorable paraît donc plus accessible que la
contraction critique elle-même.

## 8. Arbre de preuve pour atteindre $`p=0.8`$

Le résultat recherché découlerait des quatre certificats suivants.

### C0 — racines, statut : démontré

```math
\beta_{ij}>1\Longrightarrow H_S(i,j)=0.
```

### C1 — géométrie précoce, statut : démontré à fenêtre fixe

La masse des paires lointaines fusionnant sous
$`q_\triangle-\delta_q`$ disparaît. Il faut ensuite prendre
$`\delta_q\downarrow0`$ dans l'ordre annoncé.

### C2 — domination HF-S2, statut : ouvert mais réduit

Le sweep top-down est la route séquentielle la plus propre. Pour la preuve
d'impossibilité, le [corridor collapsed](20_COLLAPSED_CORRIDOR_BLACKWELL.md)
est désormais prioritaire : il tensorise exactement le théorème 4.2 à
squelette fixé et marginalise le feedback des descendants. Il reste à
construire un couplage critique/tardif des corridors et à montrer que

```math
\varepsilon_{L}^{\mathrm{geom}}
+\varepsilon_L^{\mathrm{frontière}}.
\tag{8.1}
```

tend vers zéro. Le premier terme compare les géométries et les tailles de
groupes ; le second compresse et transporte l'état-frontière du transfert
répliqué. Pour le sweep top-down non collapsed, il faut ajouter une erreur de
dépendance dynamique ; à défaut, (4.4) fournit une route plus faible par les
queues d'anti-alignement.

Sur une chaîne de cactus triangulaires, le fichier 21 ferme cette comparaison
sans hypothèse HF-S2 : les coefficients connecté et pivotal décroissent
explicitement avec le rang de fusion. Ce résultat ne couple pas encore les
corridors de la grille.

### C3 — contraction critique, établie sur cactus, ouverte sur la grille

Extraire des blocs Palm critiques dont le transfert signé répliqué possède
un coefficient strictement inférieur à un. Les buckets $`m=2`$ screenés
donnent le premier candidat avec les constantes (6.1). Le certificat final
doit porter sur un bloc collapsed complet, pas sur une fiabilité marginale.
Le fichier 20 donne la formule produit exacte lorsque les parités du corridor
sont factorisées et isole l'état de bord comme unique obstruction à cette
factorisation. Le fichier 21 résout cette obstruction sur le cactus : sous la
densité Palm fixant le LCA au rang $`q`$,

```math
A_h^{\rm LCA}(p,q)
=
\kappa_{\rm flux}(p,q)\kappa_{\rm conn}(p,q)^{h-1}
\longrightarrow0.
```

Le fichier 22 montre que le LCA seul donnerait seulement
$`\kappa_{\rm flux}(p,q)`$, indépendamment de la distance. La contraction
asymptotique exige donc les deux bras descendants ; elle ne provient pas du
seul nœud critique.

Sous C0--C3,

```math
\mathbb E[H_S(I_L,J_L)^2]\longrightarrow0,
```

puis le théorème spectral du fichier 18 interdit la weak recovery à
$`p=0.8`$ et, par dégradation BSC, pour tout $`p\le0.8`$.

## 9. Calcul certifié recommandé

Le prochain objet fini n'est pas une nouvelle simulation de grande taille.
Le cactus est maintenant calculé exactement dans le fichier 21. Il faut
construire, sur une bande triangulaire de largeur deux, le noyau de bloc

```math
Q_B(x,\epsilon,dy),
```

où $`x,y`$ contiennent la partition de frontière et
$`\epsilon\in\{-1,+1\}`$ la parité transmise. Pour chaque bloc critique et
tardif, deux calculs indépendants sont requis.

1. Énumération directe de toutes les marques et de toutes les sorties des
   heat baths.
2. Projection collapsed du corridor, recalculée par matrice de transfert et
   arithmétique d'intervalles à $`p=4/5`$ ; son rayon spectral doit être
   certifié strictement inférieur à un.

La domination peut alors être cherchée comme une dégradation du canal de
parité : conditionnellement à la frontière, la sortie tardive doit s'obtenir
en passant la sortie critique dans un BSC supplémentaire de fiabilité dans
$`[0,1]`$. Cette propriété est un programme linéaire fini. Si elle échoue,
le certificat doit être agrandi à deux ou trois niveaux ; elle ne doit pas
être remplacée par une comparaison des seules moyennes.

Les inégalités de censoring pour systèmes attractifs suggèrent cette forme
de comparaison, mais elles ne s'appliquent pas automatiquement ici : les
signes frustrés et le contre-exemple (5.1) détruisent la monotonie globale.

## 10. Audit et contre-audit

| affirmation | statut | contre-audit |
|---|---|---|
| Deux racines finales donnent $`1/2,1/2`$ | Établi exactement pour les deux ordres | le sweep doit contenir le heat bath de chaque racine |
| Le transfert répliqué exige toujours une matrice tensorielle | Faux | après moyenne stationnaire en $`\sigma`$, il vaut $`\|Kf\|_2^2`$ |
| Le niveau critique est meilleur que tout niveau tardif à coupe fixée | Vrai sans message opposé | faux point par point avec ancêtres |
| Le bucket critique domine le tardif sous sa loi complète | Établi à taille fixée par Blackwell | exige $`B\perp K\mid X`$ ; le sweep réutilise ensuite $`K`$ |
| Le niveau critique compense tout changement de taille | Faux | contre-certificat cross-size exact dans le fichier 20 |
| Blackwell se compose sur le corridor collapsed | Établi dans le fichier 20 | même squelette ; la géométrie Palm reste à coupler |
| Le LCA critique est le cas favorable sur le cactus | Établi dans le fichier 21 | repose sur les articulations |
| La conformité Nishimori du cactus tend vers $`1/2`$ | Établi dans le fichier 21 | ne prouve pas la grille |
| Le LCA seul suffit à exploiter la distance | Faux | coefficient $`\kappa_{\rm flux}`$ constant, fichier 22 |
| Le corridor complet contracte plus que le LCA seul | Établi pour le collapsed | comparaison top-down en un passage non automatique |
| Toute violation vient d'un anti-alignement quantifiable | Établi pour un nœud fixé | le sweep demande l'état-frontière dynamique |
| Les petits tores prouvent HF-S2 | Faux | ils donnent six tests compatibles avec son sens |
| HF-S2 suffit à $`p=0.8`$ | Faux | il faut aussi faire tendre le second moment critique vers zéro |
| Un bucket $`m=2`$ screené contracte strictement | Établi localement | sa multiplication le long du vrai sweep reste à certifier |
| Ajouter des updates améliore toujours la contraction | Faux hors systèmes monotones | comparer chaque programme par (1.1) |

Le progrès principal est donc une réduction plus fine : le mot
« favorable » n'est plus une hypothèse monolithique. Il se sépare en une
annulation exacte des racines, une domination de Blackwell locale, un
couplage de géométrie et de dépendance dynamique, puis un certificat de
blocs. La queue d'anti-alignement reste la route de secours si l'indépendance
conditionnelle nécessaire à Blackwell est perdue. C'est à ce niveau que
$`p=4/5`$ peut être attaqué sans perdre la dynamique hiérarchique.
