# Transfert répliqué sous Palm critique et obstruction spectrale

Cette note pousse la voie du parcours complet de la hiérarchie jusqu'à
l'objet qui peut réellement interdire la weak recovery. Le point principal
est une correction de niveau logique : la probabilité **annealed** qu'une
paire conserve sa relation ne suffit pas. Il faut contrôler le second moment
conditionnel du transfert de paire, puis le convertir en une borne sur la
plus grande valeur propre de la matrice de persistance du sweep.

Le programme obtenu est précis.

1. Le conditionnement de deux sommets distants connectés à la coupe critique
   localise leur LCA vers $`\beta_c`$.
2. Un sweep top-down ou bottom-up donne exactement un produit d'opérateurs
   tordus sur les deux bras vers le LCA.
3. Pour la weak recovery, le bon objet est le **transfert tordu répliqué** :
   les deux copies ont le même environnement $`(O,D)`$, mais des aléas de
   heat bath indépendants.
4. Une contraction de ce transfert sous la loi de Palm critique, complétée
   par une domination favorable des paires postcritiques, entraîne
   $`\mathbb E[\lambda_{\max}(H_S)/n]\to0`$.

La première cible numérique recommandée est

```math
\boxed{p_0=\frac45=0.8.}
```

Elle bat strictement la borne rigoureuse actuelle
$`p_{\mathrm{info}}=0.794659275831\ldots`$, tout en restant assez éloignée du
point de Nishimori--Ohzeki conjecturé
$`p_{\mathrm N}^{(0)}=0.835805792367\ldots`$ pour permettre un certificat
de blocs avec marge.

## 1. Conditionnement critique et localisation du LCA

Sur la grille triangulaire, posons

```math
q_p(t)=p(1-e^{-u_pt}),
\qquad
u_p=\log\frac p{1-p},
```

et

```math
q_\triangle=2\sin(\pi/18),
\qquad
\beta_c(p)=q_p^{-1}(q_\triangle).
```

Pour deux sommets $`i,j`$ à distance $`R`$, écrivons

```math
\mathcal F_R
=
\{i\leftrightarrow j\text{ dans }\Pi_{\beta_c}\}.
```

Leur niveau de coalescence est

```math
\beta_{ij}
=
\inf\{t:i\leftrightarrow j\text{ dans }\Pi_t\}.
```

### Lemme 1.1 — localisation critique, statut : établi à partir des bornes standards de percolation

Pour tout $`\varepsilon>0`$ fixé,

```math
\boxed{
\mathbb P\left(
\beta_{ij}\le\beta_c-\varepsilon
\middle|
\mathcal F_R
\right)
=
\frac{
\tau_{ij}(q_p(\beta_c-\varepsilon))
}{
\tau_{ij}(q_\triangle)
}
\longrightarrow0.
}
```

#### Preuve

L'identité est la définition de $`\beta_{ij}`$ et le couplage monotone en
coordonnée de percolation. Le paramètre du numérateur est strictement
sous-critique ; sa fonction de connexion décroît donc exponentiellement en
$`R`$. Au point critique, RSW sur des anneaux dyadiques donne une minoration
polynomiale de $`\tau_{ij}(q_\triangle)`$. Le quotient tend vers zéro.

Ce lemme justifie rigoureusement l'expression « LCA au seuil » sous le
conditionnement de paire. Il ne dit pas que tous les descendants ou tous les
ancêtres du LCA ont un niveau critique.

Sur un tore de côté $`L`$, avec $`d(i,j)\ge\rho L`$, la version uniforme
requiert les estimations RSW et sous-critiques adaptées aux conditions de
bord. C'est un passage technique séparé ; la simulation du paragraphe 8
emploie directement la loi finie correcte.

## 2. Corrélation exacte d'un parcours complet

Fixons $`O,D,\sigma`$ et un programme $`S`$ qui visite chaque feuille et
chaque nœud interne une fois, soit bottom-up, soit top-down. Posons

```math
\zeta_x=\sigma_x\sigma_x'
```

et

```math
\boxed{
H_S(i,j)
=
\mathbb E_S[\zeta_i\zeta_j\mid O,D,\sigma].
}
```

Si $`A_r`$ est la décision de flip à la mise à jour $`r`$ et si
$`\chi_r(i,j)`$ indique que la décision retourne exactement un des deux
sommets, alors

```math
H_S(i,j)
=
\mathbb E_S\left[
(-1)^{\sum_r A_r\chi_r(i,j)}
\middle|O,D,\sigma
\right].
```

En particulier,

```math
\mathbb P_S(\text{relation conservée}\mid O,D,\sigma)
=
\frac{1+H_S(i,j)}2.
```

Les nœuds qui apparaissent explicitement dans la parité sont ceux des deux
bras de $`i,j`$ vers leur LCA. Les branches latérales et les ancêtres stricts
restent toutefois dans la loi des décisions $`A_r`$ par les quatre taux
$`\Lambda_v^{ab}`$.

Avec un état-frontière suffisant $`X_r`$ et le noyau exact
$`Q_r(x,a,dx')`$, définissons

```math
(\mathcal T_r f)(x)
=
\sum_{a=0}^1
(-1)^{a\chi_r(i,j)}
\int f(x')Q_r(x,a,dx').
```

On a alors, pour toute hiérarchie finie,

```math
\boxed{
H_S(i,j)
=
\lambda\mathcal T_1\cdots\mathcal T_{H_R}\mathbf1.
}
```

L'ordre top-down ou bottom-up change les noyaux et leur produit, mais pas
l'identité de parité.

## 3. Pourquoi le premier moment annealed ne suffit pas

La probabilité demandée pour une paire aléatoire sous le conditionnement
critique est

```math
\frac12\left[
1+
\mathbb E(H_S(i,j)\mid\mathcal F_R)
\right].
```

Son retour vers $`1/2`$ ne suffit pas à interdire la weak recovery.

### Contre-exemple minimal

Supposons que l'environnement prenne deux valeurs équiprobables et que

```math
H_S(i,j)=+1
\quad\text{dans la première,}
\qquad
H_S(i,j)=-1
\quad\text{dans la seconde.}
```

Alors

```math
\mathbb E H_S(i,j)=0,
```

mais

```math
\mathbb E H_S(i,j)^2=1.
```

L'annulation du premier moment vient seulement d'une compensation entre
environnements parfaitement informatifs. Le théorème d'obstruction du
fichier 03 exige de faire disparaître le second moment.

## 4. Globalisation spectrale exacte

Soient $`I_n,J_n`$ indépendants et uniformes dans $`V_n`$, indépendamment du
hasard du sweep. Rappelons que, conditionnellement à $`O,D,\sigma`$,
$`H_S`$ est la matrice de Gram

```math
H_S=\mathbb E_S[\zeta\zeta^{\mathsf T}\mid O,D,\sigma].
```

Elle est donc positive semi-définie.

### Théorème 4.1 — critère pairwise $`L^2`$, statut : établi

Pour tout sweep admissible,

```math
\boxed{
\mathbb E\left[\frac{\lambda_{\max}(H_S)}n\right]
\le
\sqrt{
\mathbb E\left[H_S(I_n,J_n)^2\right]
}.
}
```

Ainsi,

```math
\boxed{
\mathbb E[H_S(I_n,J_n)^2]\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
}
```

#### Preuve

La norme de Frobenius donne point par point

```math
\frac{\lambda_{\max}(H_S)}n
\le
\frac{\|H_S\|_{\mathrm F}}n
=
\left(
\frac1{n^2}\sum_{i,j}H_S(i,j)^2
\right)^{1/2}.
```

Après espérance, Jensen donne

```math
\mathbb E\left[\frac{\lambda_{\max}(H_S)}n\right]
\le
\left(
\frac1{n^2}\sum_{i,j}\mathbb E[H_S(i,j)^2]
\right)^{1/2}.
```

Le terme entre parenthèses est exactement
$`\mathbb E[H_S(I_n,J_n)^2]`$. Le théorème d'obstruction du fichier 03
termine.

### Corollaire 4.2 — convergence en probabilité d'une paire uniforme

Comme $`|H_S(i,j)|\le1`$, il suffit de montrer

```math
H_S(I_n,J_n)\longrightarrow0
```

en probabilité sous la loi jointe de l'environnement et de la paire uniforme.
La bornitude donne alors la convergence $`L^2`$.

## 5. Le transfert tordu répliqué

Pour un environnement fixé,

```math
H_S(i,j)^2
=
(\lambda\otimes\lambda)
(\mathcal T_1\otimes\mathcal T_1)
\cdots
(\mathcal T_{H_R}\otimes\mathcal T_{H_R})
(\mathbf1\otimes\mathbf1).
```

Cette identité a une interprétation probabiliste simple : on exécute deux
sweeps indépendants sachant le **même** $`(O,D,\sigma)`$, et l'on tord par le
produit de leurs deux parités. L'environnement doit être partagé entre les
copies ; le remplacer par deux hiérarchies indépendantes calculerait
$`(\mathbb E H_S)^2`$ au lieu de $`\mathbb E H_S^2`$.

### 5.1 Route A — contraction uniforme de blocs

Pour un bloc $`B`$ de mises à jour, écrivons

```math
Q_B(x,\epsilon,dy)
```

pour la loi jointe de l'état final et de la parité
$`\epsilon\in\{-1,+1\}`$. Le noyau signé vaut

```math
T_B(x,dy)
=
Q_B(x,+1,dy)-Q_B(x,-1,dy).
```

Son coefficient de contraction $`L^\infty`$ exact est

```math
\boxed{
\kappa_B^{(\infty)}
=
\sup_x
\int|T_B|(x,dy)
\le1.
}
```

Si le chemin peut être partitionné en blocs $`B_1,\ldots,B_K`$, alors

```math
|H_S(i,j)|
\le
\prod_{k=1}^K\kappa_{B_k}^{(\infty)}.
```

Une proposition suffisante, directement certifiable sur un cactus ou une
bande, est donc : il existe $`\delta>0`$ tel que le nombre de blocs vérifiant
$`\kappa_B^{(\infty)}\le1-\delta`$ diverge en probabilité sous la loi de Palm
critique, tandis que les blocs exceptionnels sont négligeables. Cette route
contrôle automatiquement les deux moments.

Le contre-audit est important : un seul état-frontière déterministe peut
forcer $`\kappa_B^{(\infty)}=1`$. Il faut alors agrandir le bloc ou employer
une norme pondérée ; moyenner les fiabilités locales ne résout pas ce défaut.

### 5.2 Route B — rayon spectral annealed répliqué

La route réaliste sur la grille entière consiste à inclure dans l'état le
squelette Palm local et les deux états-frontières des sweeps. Si cet
environnement augmenté possède une limite stationnaire Markovienne, notons
$`\mathscr U_{p,S}`$ son opérateur répliqué à un bloc.

La quantité décisive est

```math
\boxed{
\rho_{2,S}(p)
=
r(\mathscr U_{p,S}),
}
```

où $`r`$ désigne le rayon spectral sur un espace pondéré à préciser. Sous
irréductibilité, contrôle des états rares et $`\rho_{2,S}(p)<1`$,

```math
\mathbb E_{\mathrm{Palm}}
[H_S(i,j)^2\mid H_R=H]
\le
C_p(\rho_{2,S}(p)+o(1))^H.
```

Ce résultat est conditionnel à la construction de la limite d'environnement.
Il indique néanmoins sans ambiguïté quelle matrice calculer sur les cactus,
les bandes et les troncatures de frontière.

### 5.3 Quenched contre annealed

Un exposant de Lyapunov typique négatif ne suffit pas pour le second moment
annealed. Des hiérarchies rares avec $`|H_S|\simeq1`$ peuvent dominer
$`\mathbb E H_S^2`$. Il faut soit :

- une contraction uniforme de blocs avec une queue quantitative ;
- le rayon spectral du transfert répliqué annealed ;
- une estimation de grandes déviations conjointe du nombre de blocs et de
  leur contraction.

## 6. Du cas critique favorable à toutes les paires

La domination HF du fichier 12 porte sur le LLR d'un unique heat bath au LCA.
Elle ne suffit pas pour un sweep complet. Il faut la remplacer par une
domination du transfert de parité.

Fixons une fenêtre gauche

```math
\mathcal F_{L,\varepsilon}
=
\left\{
d(I_L,J_L)\ge r_L,
\ \beta_c-\varepsilon\le\beta_{I_LJ_L}\le\beta_c
\right\}.
```

Posons

```math
a_{L,\varepsilon}^{S}(p)
=
\mathbb E\left[
H_S(I_L,J_L)^2
\middle|
\mathcal F_{L,\varepsilon}
\right].
```

### Hypothèse HF-S2 — domination favorable du sweep

Pour les paires lointaines fusionnant dans $`(\beta_c,1]`$, on demande

```math
\boxed{
\mathbb E\left[
H_S(I_L,J_L)^2
\middle|
d(I_L,J_L)\ge r_L,
\ \beta_c<\beta_{I_LJ_L}\le1
\right]
\le
a_{L,\varepsilon}^{S}(p)
+\delta_{L,\varepsilon}^{S}(p),
}
```

avec $`\delta_{L,\varepsilon}^{S}(p)\to0`$ dans l'ordre de limites retenu.

Cette hypothèse formalise exactement « la paire critique est le cas le plus
favorable pour garder une corrélation ». Elle est plus forte qu'une
comparaison des moyennes de $`\Lambda_v`$ et distincte de HF-LCA. Elle doit
être démontrée sur le vecteur d'état-frontière complet ou directement sur le
transfert répliqué.

### Théorème 6.1 — réduction globale au transfert critique, statut : établi sous HF-S2

Supposons :

1. les orientations des racines finales sont recolorées indépendamment et
   équitablement ;
2. $`\mathbb P(d(I_L,J_L)<r_L)\to0`$ ;
3. pour tout $`\varepsilon>0`$,

   ```math
   \mathbb P(
   d(I_L,J_L)\ge r_L,
   \ \beta_{I_LJ_L}\le\beta_c-\varepsilon
   )\to0;
   ```

4. HF-S2 est vraie ;
5. dans l'ordre $`L\to\infty`$ puis $`\varepsilon\downarrow0`$,

   ```math
   a_{L,\varepsilon}^{S}(p)
   +\delta_{L,\varepsilon}^{S}(p)
   \longrightarrow0.
   ```

Alors

```math
\mathbb E[H_S(I_L,J_L)^2]\longrightarrow0,
```

et la weak recovery est impossible.

#### Preuve

On partitionne les paires en : proches, fusion strictement sous-critique,
fenêtre critique, fusion postcritique et racines distinctes. Les deux
premières contributions sont majorées par leurs probabilités car
$`H_S^2\le1`$. La fenêtre critique est majorée par
$`a_{L,\varepsilon}^{S}`$. La partie postcritique est contrôlée par HF-S2.
Pour deux racines distinctes, les recolorations indépendantes donnent
$`H_S=0`$. Le théorème 4.1 termine.

## 7. Seuil à viser

La monotonie statistique du BSC donne une simplification stratégique. Si
l'impossibilité est démontrée à un seul $`p_0`$, elle vaut pour tout
$`p\le p_0`$ : l'observation de paramètre plus faible s'obtient en ajoutant
un bruit BSC indépendant à celle de paramètre $`p_0`$.

Les jalons utiles sont :

| cible | $`2p-1`$ | $`\beta_c(p)`$ | $`s_c(p)`$ | statut |
|---:|---:|---:|---:|---|
| $`p_{\mathrm{info}}=0.794659275831`$ | $`0.589318551663`$ | $`0.424567774256`$ | $`0.685399758609`$ | borne rigoureuse à battre |
| $`p_0=0.8`$ | $`0.6`$ | $`0.410716539196`$ | $`0.693582222752`$ | premier certificat recommandé |
| $`p_1=0.81`$ | $`0.62`$ | $`0.386167962329`$ | $`0.708903111615`$ | second jalon |
| $`p_{\mathrm N}^{(0)}=0.835805792367`$ | $`0.671611584734`$ | $`0.330008106302`$ | $`0.748439879301`$ | conjecture, objectif final |

Ici

```math
s_c(p)
=
\frac{p-q_\triangle}{1-q_\triangle}
```

est la probabilité résiduelle qu'une arête fermée au seuil soit satisfaite.

Définissons, lorsqu'une limite de transfert existe,

```math
\boxed{
p_{\mathrm{HT}}
=
\sup\left\{
p:\min_{S\in\{\mathrm{BU},\mathrm{TD}\}}
\rho_{2,S}(p)<1
\text{ et HF-S2 est vraie}
\right\}.
}
```

Une preuve de $`p_{\mathrm{HT}}\ge0.8`$ serait déjà une amélioration
rigoureuse. Une égalité formelle avec Nishimori n'a de valeur que si la loi
Palm et HF-S2 sont dérivées indépendamment ; ajuster une constante de
géométrie après coup ne constitue pas une dérivation.

## 8. Diagnostic géométrique PATH-FAC

La coordonnée $`q`$ permet une simulation sans biais supplémentaire. Donner
à chaque arête un rang uniforme $`U_e`$ produit la percolation
$`\{U_e\le q\}`$. Pour tout $`p`$, le même rang est transformé en temps par

```math
t_p(q)
=
-\frac1{u_p}\log\left(1-\frac qp\right).
```

Une composante critique de taille $`s`$ doit être sélectionnée avec le poids
$`s(s-1)`$ pour obtenir la loi d'une paire uniforme conditionnée à être
connectée. Rejeter ensuite les paires avec $`d(i,j)<\rho L`$ donne exactement
la loi Palm finie recherchée.

### Lemme 8.1 — rôle exact des buckets $`m=2`$, statut : établi dans PATH-FAC

Pour une coupe de taille deux,

```math
\boxed{
\Gamma_2(t;p)=s_p(t)\le p,
\qquad 0\le t\le\beta_c.
}
```

En effet, sous la vraie parité,
$`K=1+\mathrm{Bin}(1,s_p(t))`$. Le cas $`K=1`$ a LLR nul et le cas $`K=2`$
a LLR infini. La fiabilité moyenne vaut donc $`s_p(t)`$, qui décroît de
$`p`$ à $`s_c(p)`$.

Si $`N_{2,L}`$ est le nombre de buckets de taille deux sur les deux bras,

```math
|H_{\mathrm{PATH-FAC}}(i,j)|
\le
p^{N_{2,L}}.
```

Par conséquent,

```math
N_{2,L}\longrightarrow+\infty
\quad\text{en probabilité sous Palm}
```

entraîne la décorrélation PATH-FAC pour tout $`p<1`$ fixé, y compris après
moyenne annealed par bornitude.

### Résultats reproductibles, statut : diagnostic numérique

Le script
[`critical_pair_path_geometry.py`](computations/critical_pair_path_geometry.py)
a été exécuté avec $`200`$ hiérarchies indépendantes, une paire Palm par
hiérarchie, $`\rho=1/4`$ et une graine documentée par la ligne de commande.

```bash
for L in 8 16 32 48 64; do
  python3 research/hierarchical-swendsen-wang/computations/critical_pair_path_geometry.py \
    --side "$L" --repetitions 200 --pairs 1 --distance-fraction 0.25 \
    --p-values 0.8,0.81,0.835805792367 --seed "$((20261000+L))"
done
```

| $`L`$ | $`\mathbb E H_L`$ | $`\mathbb E N_{2,L}`$ | $`\mathbb E[q_p(\beta_{ij})/q_\triangle]`$ | corrélation PATH-FAC à $`p=0.8`$ | à $`p=p_{\mathrm N}^{(0)}`$ |
|---:|---:|---:|---:|---:|---:|
| 8 | 10.41 | 2.14 | 0.750 | 0.278 | 0.389 |
| 16 | 21.13 | 3.55 | 0.876 | 0.0615 | 0.136 |
| 32 | 38.07 | 4.91 | 0.929 | 0.0112 | 0.0384 |
| 48 | 56.18 | 6.37 | 0.949 | 0.00325 | 0.0142 |
| 64 | 68.41 | 7.32 | 0.954 | 0.00528 | 0.0142 |

Les trois signaux sont cohérents avec : localisation du LCA, croissance du
nombre de petits buckets et perte très rapide dans PATH-FAC. Les moyennes de
corrélation aux grandes tailles sont dominées par des événements rares ; la
non-monotonie entre $`L=48`$ et $`64`$ interdit tout ajustement d'exposant sur
ces seules données.

### Contre-audit décisif

Cette perte PATH-FAC pour presque tout $`p`$ ne peut pas être transportée
automatiquement à la vraie dynamique : au-dessus du vrai seuil de recovery,
le théorème 4.1 interdit précisément que le second moment du sweep admissible
disparaisse. Le mode global doit donc être porté par les dépendances entre
nœuds, les ancêtres et les branches latérales.

## 9. Diagnostic du sweep joint exact

Le script
[`joint_hierarchical_sweep.py`](computations/joint_hierarchical_sweep.py)
reconstruit le dendrogramme non marqué, affecte chaque arête physique au LCA
de ses extrémités, puis utilise à chaque proposition les facteurs exacts

```math
\prod_{v\succeq u}
\Lambda_v^{ab}
e^{(1-\beta_v)\Lambda_v^{ab}}.
```

Il visite les feuilles et les nœuds internes, bottom-up ou top-down. Pour
chaque environnement, le carré de $`H_S(i,j)`$ est estimé sans biais de
plug-in à partir de sorties Rademacher indépendantes.

La recoloration équitable des racines n'est pas omise : au nœud interne qui
est une racine finale, les états $`00`$ et $`11`$ ont le même poids, de même
que $`10`$ et $`01`$ ; le bit de flip global est donc uniforme. Une racine
réduite à une feuille a deux poids égaux. Comme la paire échantillonnée est
déjà dans une même composante critique, ce bit global ne change pas sa
relation, mais son inclusion est nécessaire au contre-audit de la dynamique.

Sur les petits tores $`L=4,6,8`$, un seul parcours donne encore des seconds
moments conditionnels typiquement entre $`0.7`$ et $`0.9`$ à $`p=0.8`$, et
encore plus hauts près de Nishimori. Les deux ordres restent proches à la
précision de ce diagnostic. Ces tailles sont trop petites pour une
extrapolation, mais l'écart avec PATH-FAC est de plusieurs ordres de
grandeur.

Un audit apparié à $`L=8`$ conserve les mêmes $`30`$ hiérarchies pour un,
deux et quatre parcours, avec $`100`$ sorties par environnement et la graine
20261804. Il donne :

| parcours à $`D`$ fixé | second moment bottom-up | second moment top-down |
|---:|---:|---:|
| 1 | $`0.869\pm0.041`$ | $`0.870\pm0.036`$ |
| 2 | $`0.845\pm0.044`$ | $`0.842\pm0.039`$ |
| 4 | $`0.819\pm0.049`$ | $`0.846\pm0.041`$ |

Les erreurs affichées sont des erreurs standards entre environnements, pas
des intervalles de confiance uniformes. Elles ne montrent pas de séparation
significative entre les ordres. Surtout, répéter quatre fois le sweep ne fait
pas disparaître le mode persistant à cette taille. Les flux pseudo-aléatoires
de géométrie et de heat bath sont séparés dans le script afin que changer le
nombre de parcours ne change pas l'échantillon de hiérarchies.

Ce résultat oriente le travail : un sweep unique pourrait être trop
« collant » pour battre information-percolation, ou sa contraction pourrait
n'apparaître qu'à une échelle inaccessible à ces tores. Il faut mesurer des
blocs exacts et leur spectre avant de lancer des simulations plus grandes.

## 10. Audit et contre-audit

| Affirmation | Statut | Conséquence |
|---|---|---|
| La paire connectée à $`\beta_c`$ a son LCA près de $`\beta_c`$ | Établi pour deux points distants | le conditionnement critique est légitime |
| $`\mathbb E H_S\to0`$ interdit la weak recovery | Faux | il faut $`\mathbb E H_S^2\to0`$ |
| Le transfert répliqué utilise deux environnements indépendants | Faux | même environnement, deux aléas de sweep |
| Une contraction pairwise $`L^2`$ moyenne contrôle le spectre global | Établi | théorème 4.1 |
| HF-LCA implique HF-S2 | Non établi | le sweep complet demande un nouveau lemme |
| Un nombre divergent de buckets $`m=2`$ tue PATH-FAC | Établi | borne $`p^{N_2}`$ |
| PATH-FAC décrit la vraie dynamique | Réfuté numériquement et non démontré théoriquement | conserver l'état-frontière |
| Bottom-up est toujours meilleur que top-down | Faux en général | comparer les deux rayons spectraux |
| Les données finies prouvent une borne à $`p=0.8`$ | Faux | elles choisissent la bonne cible de preuve |
| Certifier $`\rho_{2,S}(0.8)<1`$ suffit seul | Faux | il faut aussi HF-S2 et le contrôle des états rares |

## 11. Ordre de travail recommandé

1. **Cactus de triangles.** Construire exactement le transfert répliqué pour
   les deux ordres et vérifier l'implémentation par énumération complète.
2. **Blocs de deux à quatre niveaux.** Calculer
   $`\kappa_B^{(\infty)}`$ et identifier les états qui forcent la norme à
   un. Cela dira quelle information de frontière manque.
3. **Jalon $`p=4/5`$.** Chercher un certificat rationnel ou par arithmétique
   d'intervalles avec marge, avant toute optimisation en $`p`$.
4. **Géométrie Palm.** Prouver ou réfuter
   $`N_{2,L}\to\infty`$, puis mesurer la fréquence des blocs réellement
   contractants pour le transfert joint, et pas seulement pour PATH-FAC.
5. **HF-S2.** Coupler les environnements critique et postcritique au niveau
   de l'état répliqué ; une monotonie des seuls temps ou des seuls
   $`\Lambda_v^{00}`$ est insuffisante.
6. **Globalisation.** Appliquer le théorème 6.1, puis utiliser la dégradation
   du BSC pour étendre le résultat à tout $`p\le p_0`$.
7. **Montée du seuil.** Répéter à $`0.81`$, puis approcher
   $`p_{\mathrm N}^{(0)}`$ seulement si les marges de contraction restent
   visibles.

Le verdict actuel est donc double. La géométrie du chemin critique est très
favorable à une perte locale accumulée, mais un sweep joint transporte un
état global que PATH-FAC supprime. La meilleure chance de battre les seuils
connus consiste à certifier la disparition de ce mode global dans le
**transfert répliqué**, d'abord à $`p=4/5`$.

La [représentation du sweep par projections](19_FAVORABLE_SWEEP_PROJECTIONS.md)
donne ensuite une identité $`L^2`$ exacte, prouve l'annulation des racines
distinctes et prouve l'ordre de Blackwell critique/postcritique à taille de
bucket fixée. Le transport à la géométrie et à l'état-frontière reste ouvert.

Le [corridor collapsed](20_COLLAPSED_CORRIDOR_BLACKWELL.md) renforce cette
étape : il est $`L^2`$-optimal parmi les sweeps des mêmes nœuds et tensorise
Blackwell sur tout corridor fixé, même avec un prior corrélé des parités.
La priorité passe ainsi du feedback d'un sweep arbitraire au transfert fini
de l'état de bord sous la loi Palm critique.

Les références primaires utilisées pour la localisation critique, les lois
Palm/IIC et le positionnement par rapport à la synchronisation sur grille
sont recensées dans [l'état de l'art](LITERATURE.md).
