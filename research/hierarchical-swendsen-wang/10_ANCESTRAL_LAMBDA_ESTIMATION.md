# Estimer les $`\Lambda_v`$ au-dessus du LCA critique

Ce fichier reprend le problème dans l'ordre exact des slides 31--33 de la
[présentation du 16 juillet 2026](../../beamer-presentation-reunion-2026-07-16/Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf).
Pour $`u=\mathrm{LCA}(i,j)`$, le facteur du nœud $u$ est simple : un
flip de parité impaire remplace $`\Lambda_u`$ par $`T_u-\Lambda_u`$. En
revanche, le heat bath exact utilise

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
```

La difficulté n'est donc pas de calculer le seul bucket critique $u$, mais
d'estimer simultanément le vecteur à quatre états

```math
\mathcal L_u
:=
\left(
  \bigl(\Lambda_v^{00},\Lambda_v^{01},
        \Lambda_v^{10},\Lambda_v^{11}\bigr)
\right)_{v\succ u}
```

sur toute la chaîne ancestrale. La localisation favorable
$`\beta_u\simeq\beta_c`$ fixe le bas de cette chaîne ; elle ne supprime aucun
de ses ancêtres.

Les résultats nouveaux de ce fichier sont les suivants.

1. **Établi, volume fini.** Conditionnellement au squelette de Kruskal non
   marqué et aux modules des poids, la loi de chaque bucket ancestral est une
   mixture explicite de Poisson-binomiales pondérées. L'arête gagnante n'est
   uniforme que dans le modèle homogène.
2. **Établi, volume fini.** Les moyennes et la matrice de covariance des trois
   poids satisfaits, puis des quatre $`\Lambda_v^{ab}`$, ont des formules
   exactes.
3. **Établi, volume fini.** Dans le cas homogène, les quatre moyennes et
   variances se ferment en fonction de
   $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)`$.
4. **Établi, déterministe.** Une fonctionnelle explicite contrôle l'erreur sur
   $`B_u`$ lorsque des ancêtres sont tronqués.
5. **À prouver.** Sur la grille entière, le verrou restant est la loi du
   squelette ancestral sous la mesure biaisée par une paire lointaine dont le
   LCA tombe dans la fenêtre critique.
6. **Établi sous une domination explicite.** Le
   [théorème de réduction favorable](12_FAVORABLE_HIERARCHICAL_REDUCTION.md)
   montre que l'annulation de la fiabilité critique implique l'impossibilité
   globale dès que l'expérience postcritique est dominée par l'oracle où la
   paire se sépare au seuil.

## 1. Expérience favorable et quantité à calculer

Soit $`G_L`$ une exhaustion finie, et soit $`r_L\to\infty`$. Pour ne pas
confondre « paire uniforme » et « paire lointaine », posons

```math
N_{u,L}^{\mathrm{far}}
:=
\#\left\{
  (x,y)\in C_{u,1}\times C_{u,2}:d(x,y)\ge r_L
\right\}.
```

Sans la contrainte de distance, ce poids vaut exactement
$`|C_{u,1}||C_{u,2}|`$. Avec la contrainte, cette égalité demande un argument
géométrique et ne sera pas supposée.

L'expérience n'est non vide avant la coupe $`1`$ que si le seuil géométrique
est accessible. Sur la grille triangulaire homogène,

```math
\beta_c(p)\le1
\quad\Longleftrightarrow\quad
q_p(1)=2p-1\ge q_c
\quad\Longleftrightarrow\quad
p\ge p_{\mathrm{SW}}:=\frac{1+q_c}{2}.
```

Pour $`p<p_{\mathrm{SW}}`$, il n'existe donc pas de fusion au seuil de
percolation dans le dendrogramme tronqué à $`1`$.

Pour fixer les quantificateurs, prenons d'abord

```math
I_\varepsilon
:=
[\beta_c(p)-\varepsilon,\beta_c(p)+\varepsilon]\cap[0,1].
```

Une fenêtre à l'échelle proche-critique dépendant de $L$ peut ensuite la
remplacer, mais son échelle doit être annoncée explicitement. La loi du
squelette vu depuis la paire favorable est

```math
\mathbb E_{L,\varepsilon}^{\star}[F]
:=
\frac{
\mathbb E\!\left[
  \sum_{u:\,\beta_u\in I_\varepsilon}
  N_{u,L}^{\mathrm{far}}F(u,D,|W|)
\right]
}{
\mathbb E\!\left[
  \sum_{u:\,\beta_u\in I_\varepsilon}
  N_{u,L}^{\mathrm{far}}
\right]
},
```

lorsque le dénominateur est non nul et pour toute $`F`$ mesurable par rapport
aux données non marquées et aux modules des poids. Conditionner sur
$`\beta_u=\beta_c`$ n'est pas une formulation correcte en volume fini : les
temps sont continus. Il faut spécifier une fenêtre, puis l'ordre des limites.

Notons $`\mathscr S_u`$ le squelette ancestral non marqué :

```math
\mathscr S_u
=
\left(
  E_v^{(0)},E_v^{(1)},E_v^{(2)},
  (|W_e|)_{e\in E_v},\beta_v
\right)_{v\succeq u},
```

sans les marques de satisfaction et sans l'identité des arêtes gagnantes. La
fiabilité favorable exacte se désintègre sous la forme

```math
\mathcal G_{L,\varepsilon}^{\mathrm{fav}}(p)
=
\mathbb E_{L,\varepsilon}^{\star}
\left[
  \Gamma_p(\mathscr S_u)
\right],
\qquad
\Gamma_p(\mathscr S_u)
:=
\mathbb E\left[\eta_u\mid\mathscr S_u\right].
```

Le noyau conditionnel qui définit $`\Gamma_p`$ est calculé exactement dans
les sections 3--5. Toute difficulté asymptotique est ainsi reportée, sans être
cachée, dans la loi $`\mathbb P_{L,\varepsilon}^{\star}(d\mathscr S_u)`$.

Notons $`\mathcal K_p(\mathscr S_u,dY)`$ le produit des noyaux de marques de
la section 3. Dès qu'une fonctionnelle dépend des marques — c'est le cas de
$`\ell_v`$ et de $`\mathcal R_u`$ — la mesure pertinente est l'extension

```math
\widehat{\mathbb P}_{L,\varepsilon}^{\star}
(d\mathscr S_u,dY)
:=
\mathbb P_{L,\varepsilon}^{\star}(d\mathscr S_u)
\mathcal K_p(\mathscr S_u,dY).
```

Cette distinction empêche d'attribuer au squelette seul une quantité qui
dépend encore des marques conditionnelles.

On peut définir un **seuil favorable ancestral** à partir de la non-disparition
de $`\mathcal G_{L,\varepsilon}^{\mathrm{fav}}`$, après avoir fixé l'ordre des
limites. Pour une preuve d'impossibilité, son rôle est plus précis qu'un simple
benchmark : si le lemme de domination HF du fichier 12 est établi, toute la
contribution postcritique est majorée par cet oracle critique. Son annulation
devient alors une condition suffisante d'impossibilité. Sans HF, le succès ou
l'échec sous ce seul conditionnement ne détermine pas le seuil global.

## 2. Statistiques minimales d'un ancêtre

Fixons

```math
u:C_u=C_1\mathbin{\dot\cup}C_2
```

et un ancêtre strict $`v\succ u`$. Le fils de $v$ contenant $`C_u`$ est noté
$`P_v`$, et son autre fils $`S_v`$. Le bucket de $v$ se partitionne en

```math
E_v
=
E_v^{(0)}\mathbin{\dot\cup}
E_v^{(1)}\mathbin{\dot\cup}
E_v^{(2)},
```

suivant que l'extrémité située dans $`P_v`$ appartient à
$`P_v\setminus C_u`$, à $`C_1`$ ou à $`C_2`$. Posons

```math
T_{v,r}=\sum_{e\in E_v^{(r)}}|W_e|,
\qquad
\lambda_{v,r}
=\sum_{e\in E_v^{(r)}}
|W_e|\mathbf1_{\{e\text{ satisfaite}\}},
\qquad
X_{v,r}=2\lambda_{v,r}-T_{v,r}.
```

Avec $`\epsilon_a=1-2a`$ et $`\epsilon_b=1-2b`$,

```math
\boxed{
\Lambda_v^{ab}
=
\frac12\left[
T_{v,0}+T_{v,1}+T_{v,2}
+X_{v,0}+\epsilon_aX_{v,1}+\epsilon_bX_{v,2}
\right].
}
```

Cette identité donne le bon objet à estimer. Connaître seulement
$`\Lambda_v^{00}`$ ne suffit pas : il faut au minimum le niveau commun
$`X_{v,0}`$ et les deux déséquilibres retournables
$`X_{v,1},X_{v,2}`$.

Le [certificat de majorité hiérarchique](14_CRITICAL_COMPONENT_BOUNDARY.md)
donne une première cible probabiliste fermée. Sous a priori uniforme, si le
nœud local possède une majorité conforme stricte et si

```math
X_{v,1}\ge0,
\qquad
X_{v,2}\ge0
```

pour chaque ancêtre, alors le heat bath préfère nécessairement la parité
conforme. La preuve utilise la convexité de
$`x\mapsto xe^{(1-\beta_v)x}`$ et la stabilité par produit d'un cône de Walsh
positif. L'estimation exacte reste plus riche : le critère quatre états peut
réussir même lorsqu'un de ces déséquilibres est négatif.

## 3. Course pondérée conditionnelle exacte

### Hypothèse de canal

Conditionnellement aux modules $`w_e:=|W_e|`$, supposons que les marques

```math
Y_e:=\mathbf1_{\{e\text{ satisfaite}\}}
```

soient indépendantes, de paramètres $`\pi_e\in(0,1)`$, puis que

```math
T_e
=
\begin{cases}
\mathrm{Exp}(w_e),&Y_e=1,\\
+\infty,&Y_e=0.
\end{cases}
```

Dans un canal de Nishimori paramétré par le LLR, on a
$`\pi_e=\mathrm{logistic}(w_e)`$. Le cas homogène correspond à
$`w_e=u_p=\log(p/(1-p))`$ et $`\pi_e=p`$.

### Théorème 3.1 — noyau exact d'un bucket non marqué

Fixons le bucket $`E_v`$ et conditionnons son minimum non marqué à valoir
$`\beta_v`$, au sens de la désintégration par rapport à sa densité.
Définissons

```math
s_{v,e}
:=
\mathbb P(Y_e=1\mid T_e>\beta_v,w_e)
=
\frac{\pi_e e^{-w_e\beta_v}}
{1-\pi_e+\pi_e e^{-w_e\beta_v}}.
```

Dans le canal de Nishimori,

```math
s_{v,e}
=
\mathrm{logistic}\!\left(w_e(1-\beta_v)\right).
```

Il existe une arête gagnante latente $`G_v\in E_v`$ telle que

```math
\boxed{
\rho_{v,e}
:=
\mathbb P(G_v=e\mid\mathscr S_u)
=
\frac{w_es_{v,e}}
{\sum_{f\in E_v}w_fs_{v,f}}.
}
```

Conditionnellement à $`G_v=e`$ :

- $`Y_e=1`$ ;
- les $`Y_f`$, $`f\ne e`$, sont indépendantes et suivent
  $`\mathrm{Bernoulli}(s_{v,f})`$ ;
- les buckets distincts sont indépendants conditionnellement au squelette
  complet.

Ainsi, pour $`r\in\{0,1,2\}`$,

```math
\boxed{
\lambda_{v,r}
\ \stackrel d=\
\mathbf1_{\{G_v\in E_v^{(r)}\}}w_{G_v}
+
\sum_{\substack{f\in E_v^{(r)}\\f\ne G_v}}
w_fB_{v,f},
\qquad
B_{v,f}\sim\mathrm{Bernoulli}(s_{v,f}).
}
```

Il s'agit d'une mixture finie explicite de Poisson-binomiales pondérées.

La transformée exponentielle jointe ferme même la loi complète. Si $`r(e)`$
est le groupe de $e$, alors

```math
\boxed{
M_v(t_0,t_1,t_2)
:=
\mathbb E\left[
e^{\sum_{r=0}^2t_r\lambda_{v,r}}
\middle|\mathscr S_u
\right]
=
\sum_{e\in E_v}\rho_{v,e}
e^{t_{r(e)}w_e}
\prod_{f\ne e}
\left(
1-s_{v,f}+s_{v,f}e^{t_{r(f)}w_f}
\right).
}
```

Les buckets étant conditionnellement indépendants, la transformée de toute la
chaîne est $`\prod_{v\succ u}M_v`$. Dans le cas homogène, avec
$`s_v=\mathrm{logistic}(u_p(1-\beta_v))`$, la fonction génératrice des comptes
$`K_{v,r}=\lambda_{v,r}/u_p`$ devient le polynôme

```math
\boxed{
P_v(z_0,z_1,z_2)
=
\sum_{r:\,m_{v,r}>0}
\frac{m_{v,r}}{m_v}
z_r
(1-s_v+s_vz_r)^{m_{v,r}-1}
\prod_{\ell\ne r}
(1-s_v+s_vz_\ell)^{m_{v,\ell}}.
}
```

L'extraction de ses coefficients donne une énumération exacte sans revenir
aux $`2^{m_v}`$ configurations de marques. Pour une chaîne courte, le produit
de ces polynômes et le calcul quatre états donnent donc
$`\Gamma_p(\mathscr S_u)`$ exactement.

### Preuve

La densité pour que $e$ gagne à l'instant $t$ vaut

```math
\pi_ew_ee^{-w_et}
\prod_{f\ne e}
\left(1-\pi_f+\pi_fe^{-w_ft}\right).
```

En factorisant le produit des probabilités de survie de toutes les arêtes, le
facteur dépendant de $e$ devient

```math
\frac{\pi_ew_ee^{-w_et}}
{1-\pi_e+\pi_ee^{-w_et}}
=w_es_e(t).
```

La normalisation donne $`\rho_{v,e}`$. Sachant que $e$ gagne, les autres
horloges sont seulement conditionnées par $`T_f>t`$ ; leurs marques restent
indépendantes, avec les paramètres $`s_f(t)`$. Enfin, les buckets sont des
ensembles d'arêtes disjoints et les contraintes de minimum portent bucket par
bucket. Cela prouve aussi la factorisation le long de la chaîne.

### Contre-audit immédiat

La gagnante est uniforme si tous les $`w_es_{v,e}`$ sont égaux. C'est le cas
homogène, mais pas le cas pondéré général. Remplacer $`\rho_{v,e}`$ par
$`1/|E_v|`$ en présence de poids variables biaise déjà la moyenne de
$`\lambda_{v,r}`$.

## 4. Moyenne et covariance exactes

Posons

```math
\mu_{v,r}
:=
\sum_{e\in E_v^{(r)}}w_es_{v,e},
\qquad
d_{e,r}
:=
\mathbf1_{\{e\in E_v^{(r)}\}}
w_e(1-s_{v,e}),
\qquad
\bar d_{v,r}
:=\sum_{e\in E_v}\rho_{v,e}d_{e,r}.
```

Le théorème précédent donne exactement

```math
\boxed{
\mathbb E[\lambda_{v,r}\mid\mathscr S_u]
=
\mu_{v,r}+\bar d_{v,r}.
}
```

Pour le vecteur $`\boldsymbol\lambda_v=(\lambda_{v,0},\lambda_{v,1},
\lambda_{v,2})`$,

```math
\boxed{
\mathrm{Cov}(\boldsymbol\lambda_v\mid\mathscr S_u)
=
\mathrm{diag}_{r}
\left(
\sum_{e\in E_v^{(r)}}
(1-\rho_{v,e})w_e^2s_{v,e}(1-s_{v,e})
\right)
+
\mathrm{Cov}_{G_v\sim\rho_v}(d_{G_v}).
}
```

Le second terme est le coût exact de la marginalisation de la gagnante. En
particulier, pour $`r\ne s`$,

```math
\mathrm{Cov}(\lambda_{v,r},\lambda_{v,s}\mid\mathscr S_u)
=-\bar d_{v,r}\bar d_{v,s}.
```

Les groupes sont donc indépendants **conditionnellement à la gagnante**, mais
pas après son oubli.

Pour transporter ces moments vers les quatre taux, posons

```math
c^{ab}=(1,1-2a,1-2b),
\qquad
k^{ab}=aT_{v,1}+bT_{v,2}.
```

Alors

```math
\Lambda_v^{ab}=k^{ab}+c^{ab}\boldsymbol\lambda_v,
```

et donc

```math
\boxed{
\begin{aligned}
\mathbb E[\Lambda_v^{ab}\mid\mathscr S_u]
&=k^{ab}+c^{ab}\mathbb E[\boldsymbol\lambda_v\mid\mathscr S_u],\\
\mathrm{Cov}(\Lambda_v^{ab},\Lambda_v^{cd}\mid\mathscr S_u)
&=c^{ab}
\mathrm{Cov}(\boldsymbol\lambda_v\mid\mathscr S_u)
(c^{cd})^{\mathsf T}.
\end{aligned}
}
```

Ce sont des estimateurs conditionnels exacts de tous les $`\Lambda_v`$ au
dessus de $u$, pour des poids hétérogènes.

## 5. Fermeture homogène

Supposons maintenant $`w_e=u_p`$ et $`\pi_e=p`$. Écrivons

```math
m_{v,r}=|E_v^{(r)}|,
\qquad
m_v=\sum_{r=0}^2m_{v,r},
\qquad
q_{v,r}=\frac{m_{v,r}}{m_v},
```

```math
s_v=\mathrm{logistic}(u_p(1-\beta_v)),
\qquad
h_v=2s_v-1,
\qquad
\alpha_v=h_v+\frac{1-h_v}{m_v}.
```

La gagnante est uniforme parmi les $`m_v`$ arêtes. Pour chaque groupe,

```math
\boxed{
\frac{1}{u_p}\mathbb E[X_{v,r}\mid\mathscr S_u]
=m_{v,r}h_v+q_{v,r}(1-h_v).
}
```

Lorsque $`m_{v,r}>0`$, la fraction de déséquilibre satisfait vérifie donc

```math
\boxed{
\mathbb E\left[
\frac{X_{v,r}}{u_pm_{v,r}}
\middle|\mathscr S_u
\right]
=\alpha_v.
}
```

Le terme $`(1-h_v)/m_v`$ est la correction exacte de la gagnante ; il ne doit
pas être jeté dans les petits buckets.

La covariance des déséquilibres est

```math
\boxed{
\frac{\mathrm{Cov}(X_{v,r},X_{v,s}\mid\mathscr S_u)}{4u_p^2}
=
\mathbf1_{\{r=s\}}
(m_{v,r}-q_{v,r})s_v(1-s_v)
+
(1-s_v)^2
\left(
\mathbf1_{\{r=s\}}q_{v,r}-q_{v,r}q_{v,s}
\right).
}
```

Définissons la taille signée vue par l'état $`ab`$,

```math
A_v^{ab}
:=
m_{v,0}+(1-2a)m_{v,1}+(1-2b)m_{v,2}.
```

Alors les quatre moyennes se ferment en une ligne :

```math
\boxed{
\mathbb E[\Lambda_v^{ab}\mid\mathscr S_u]
=
\frac{u_p}{2}
\left(m_v+\alpha_vA_v^{ab}\right).
}
```

Leur variance exacte vaut

```math
\boxed{
\mathrm{Var}(\Lambda_v^{ab}\mid\mathscr S_u)
=
u_p^2\left[
s_v(1-s_v)(m_v-1)
+
(1-s_v)^2
\left(1-\left(\frac{A_v^{ab}}{m_v}\right)^2\right)
\right].
}
```

Plus généralement, si $`c=c^{ab}`$ et $`d=c^{cd}`$,

```math
\frac{\mathrm{Cov}(\Lambda_v^{ab},\Lambda_v^{cd}
\mid\mathscr S_u)}{u_p^2}
=
s_v(1-s_v)\left(1-\frac1{m_v}\right)
\sum_{r=0}^2m_{v,r}c_rd_r
+
(1-s_v)^2
\left[
\sum_{r=0}^2q_{v,r}c_rd_r
-
\left(\sum_{r=0}^2q_{v,r}c_r\right)
\left(\sum_{r=0}^2q_{v,r}d_r\right)
\right].
```

### Conséquence dans le cas favorable critique

Sur l'événement $`\beta_u\simeq\beta_c`$, tout ancêtre strict satisfait

```math
\beta_v>\beta_u,
\qquad
0\le h_v<h_p(\beta_u)\simeq h_p(\beta_c).
```

Pour $`p\ge p_{\mathrm{SW}}`$, de sorte que $`\beta_c(p)\le1`$, le point
critique triangulaire vérifie

```math
h_p(\beta_c)
=
\frac{2p-1-q_c}{1-q_c}
=
\frac{2(p-p_{\mathrm{SW}})}{1-q_c}.
```

Le rapport signal/bruit du déséquilibre $`X_{v,r}`$ est gouverné, à premier
ordre, par

```math
m_{v,r}h_v^2.
```

- si $`m_{v,r}h_v^2\gg1`$, le signe moyen du groupe est concentré ;
- si $`m_{v,r}h_v^2=O(1)`$, les fluctuations sont du même ordre que le
  signal ;
- lorsque $`\beta_v\uparrow1`$, $`h_v\downarrow0`$ et la substitution par la
  moyenne devient particulièrement dangereuse.

Ainsi, localiser $u$ à $`\beta_c`$ ne rend pas homogènes les ancêtres : leurs
paramètres parcourent tout l'intervalle entre $`\beta_c`$ et $`1`$.

## 6. Concentration simultanée conditionnelle

Conditionnons aussi la gagnante à être $`G_v=g`$. Posons

```math
V_{v,r}(g)
=
\sum_{\substack{e\in E_v^{(r)}\\e\ne g}}
w_e^2s_{v,e}(1-s_{v,e}),
\qquad
Q_{v,r}(g)
=
\sum_{\substack{e\in E_v^{(r)}\\e\ne g}}w_e^2,
```

et $`b_{v,r}(g)`$ égal au plus grand poids de cette somme. Si la somme est
vide, $`\lambda_{v,r}`$ est exactement égale à son espérance conditionnelle
et le rayon est nul. Sinon, pour tout $`x>0`$,

```math
\mathbb P\left(
\left|\lambda_{v,r}
-\mathbb E[\lambda_{v,r}\mid\mathscr S_u,G_v=g]\right|
\ge\sqrt{\frac{Q_{v,r}(g)x}{2}}
\middle|\mathscr S_u,G_v=g
\right)
\le2e^{-x},
```

et la borne de Bernstein, plus fine lorsque les $`s_{v,e}`$ sont proches de
$`0`$ ou $`1`$, donne

```math
\mathbb P\left(
\left|\lambda_{v,r}
-\mathbb E[\lambda_{v,r}\mid\mathscr S_u,G_v=g]\right|
\ge
\sqrt{2V_{v,r}(g)x}+\frac{2b_{v,r}(g)x}{3}
\middle|\mathscr S_u,G_v=g
\right)
\le2e^{-x}.
```

Pour une chaîne de $`H`$ ancêtres, prendre
$`x=\log(6H/\delta)`$ donne des intervalles simultanés sur les $`3H`$
groupes avec probabilité au moins $`1-\delta`$. Comme la gagnante n'est pas
observée dans $`\mathscr S_u`$, deux options sont exactes :

1. conserver la mixture en $`G_v`$ dans le calcul ;
2. prendre l'enveloppe des intervalles sur toutes les gagnantes possibles.

Les intervalles des trois groupes se transportent par la formule affine vers
les quatre $`\Lambda_v^{ab}`$. Une borne inférieure strictement positive permet
ensuite de contrôler les logarithmes. Si un intervalle atteint zéro, il faut
revenir aux quatre poids exacts et non diviser par un taux possiblement nul.

## 7. Ce qu'il faut sommer : les contrastes, pas les niveaux communs

Le heat bath est invariant par addition du même log-poids aux quatre états.
Dans cette section, les quatre poids de l'a priori sont supposés strictement
positifs et les quatre taux de chaque ancêtre considéré sont supposés positifs.
Si cette hypothèse échoue, le calcul quatre états reste exact mais la réduction
logarithmique et le certificat ci-dessous ne s'appliquent pas.

Pour un ancêtre, posons

```math
x_v=\Lambda_v^{00},
\qquad
\Delta_{v,1}=T_{v,1}-2\lambda_{v,1}=-X_{v,1},
\qquad
\Delta_{v,2}=T_{v,2}-2\lambda_{v,2}=-X_{v,2}.
```

Les quatre coins sont

```math
x_v,
\quad x_v+\Delta_{v,2},
\quad x_v+\Delta_{v,1},
\quad x_v+\Delta_{v,1}+\Delta_{v,2}.
```

Avec $`\phi_v(z)=\log z+(1-\beta_v)z`$, les champs de Walsh possèdent les
identités intégrales orientées exactes

```math
h_{v,1}
=
-\frac14
\int_0^{\Delta_{v,1}}
\left[
\phi_v'(x_v+s)
+\phi_v'(x_v+\Delta_{v,2}+s)
\right]ds,
```

```math
h_{v,2}
=
-\frac14
\int_0^{\Delta_{v,2}}
\left[
\phi_v'(x_v+t)
+\phi_v'(x_v+\Delta_{v,1}+t)
\right]dt,
```

et le couplage direct vérifie

```math
\boxed{
4J_v
=
\log\frac{x_v(x_v+\Delta_{v,1}+\Delta_{v,2})}
{(x_v+\Delta_{v,1})(x_v+\Delta_{v,2})}
=
-\int_0^{\Delta_{v,1}}
 \int_0^{\Delta_{v,2}}
\frac{dt\,ds}{(x_v+s+t)^2}.
}
```

Le terme linéaire $`(1-\beta_v)\Lambda_v^{ab}`$ contribue aux champs, mais
pas à $`J_v`$.

Soit

```math
\ell_v:=\min_{a,b}\Lambda_v^{ab}>0.
```

Comme $`\phi_v'(z)=(1-\beta_v)+1/z`$,

```math
|h_{v,1}|
\le
\frac{|\Delta_{v,1}|}{2}
\left(1-\beta_v+\frac1{\ell_v}\right),
```

```math
|h_{v,2}|
\le
\frac{|\Delta_{v,2}|}{2}
\left(1-\beta_v+\frac1{\ell_v}\right),
\qquad
|J_v|
\le
\frac{|\Delta_{v,1}\Delta_{v,2}|}{4\ell_v^2}.
```

Pour un ensemble $`I`$ d'ancêtres omis, définissons

```math
\boxed{
\mathcal R_u(I)
:=
\sum_{v\in I}
\left[
(|\Delta_{v,1}|+|\Delta_{v,2}|)
\left(1-\beta_v+\frac1{\ell_v}\right)
+
\frac{|\Delta_{v,1}\Delta_{v,2}|}{2\ell_v^2}
\right].
}
```

La fonction $`\log\cosh`$ étant $`1`$-Lipschitz,

```math
\boxed{
|B_u-B_u^{(-I)}|
\le
\mathcal R_u(I).
}
```

Cette fonctionnelle est le bon certificat de troncature. Elle montre aussi
pourquoi une borne sur les seuls $`\Lambda_v^{00}`$ n'est pas suffisante : les
deux déséquilibres et le plus petit des quatre coins interviennent
explicitement.

## 8. Réduction exacte du verrou asymptotique

Ordonnons les ancêtres

```math
u=v_0\prec v_1\prec\cdots\prec v_{H_L}.
```

Le programme favorable critique serait fermé par les trois résultats
suivants. G1 porte sur $`\mathbb P_{L,\varepsilon}^{\star}`$ ; G2 et G3, qui
dépendent des marques, portent sur
$`\widehat{\mathbb P}_{L,\varepsilon}^{\star}`$.

### G1 — convergence du squelette marqué par les tailles

Pour tout $`K`$ fixé, établir la convergence de

```math
\left(
E_u,(|W_e|)_{e\in E_u},\beta_u;
\left(
\beta_{v_k},
(|W_e|,r_e)_{e\in E_{v_k}}
\right)_{1\le k\le K}
\right),
```

ou, dans le cas homogène, au minimum de

```math
\left(
\beta_u,m_u;
\left(
\beta_{v_k},m_{v_k,0},m_{v_k,1},m_{v_k,2}
\right)_{1\le k\le K}
\right).
```

### G2 — sommabilité de la queue ancestrale

Pour tout $`\zeta>0`$, viser

```math
\lim_{K\to\infty}
\limsup_{\varepsilon\downarrow0}
\limsup_{L\to\infty}
\widehat{\mathbb P}_{L,\varepsilon}^{\star}\!\left(
\mathcal R_u(\{v_k:k>K\})>\zeta
\right)
=0.
```

Une version en espérance, uniformément intégrable, est encore plus forte et
permet de passer directement aux fiabilités moyennes.

### G3 — contrôle des coins proches de zéro

Soit montrer que les événements $`\ell_{v_k}\simeq0`$ ont une contribution
négligeable, soit les conserver exactement dans une limite quatre états. Sans
ce point, le passage à $`\log\Lambda_v^{ab}`$ n'est pas uniforme.

Sous G1--G3, le noyau de marques des sections 3--6 et la borne de la section 7
réduisent la limite de $`\mathcal G_{L,\varepsilon}^{\mathrm{fav}}`$ à un
calcul fini, puis à une limite $`K\to\infty`$. Ce serait une véritable
solution du problème ancestral. Les seules formules locales au nœud $u$ ne
peuvent pas remplacer G1--G3.

Le lemme 6.1 du [fichier 12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md) transporte
quantitativement cette approximation jusqu'à la fiabilité : si
$`B_u^{(K)}`$ conserve les $`K`$ premiers ancêtres, alors

```math
\left|
\tanh^2\left(\frac{\ell_u^{\mathrm{crit}}+B_u}{2}\right)
-
\tanh^2\left(\frac{\ell_u^{\mathrm{crit}}+B_u^{(K)}}{2}\right)
\right|
\le
\min\left(1,\frac{2\mathcal R_u^{(>K)}}{3\sqrt3}\right).
```

Ainsi G1--G3 ferment l'oracle favorable. Pour en déduire une borne globale de
weak recovery sans supposer que tous les temps LCA se concentrent, il reste à
prouver G4 : la domination HF entre la chaîne postcritique et la chaîne
critique.

### Quelle géométrie doit être estimée ?

Le squelette est sélectionné deux fois :

1. Kruskal choisit les buckets et leurs temps par une règle de minimum ;
2. la paire lointaine critique biaise les nœuds par
   $`N_{u,L}^{\mathrm{far}}`$.

La loi pertinente n'est donc ni celle d'un bucket déterministe, ni celle d'un
nœud uniforme du dendrogramme. Pour la grille, il faut contrôler le bord des
deux composantes critiques $`C_1,C_2`$, puis leurs incidences avec chaque
frère ancestral $`S_v`$. Les trois nombres $`m_{v,r}`$ sont précisément ces
incidences.

## 9. Voies de résolution

### 9.1 Cactus et arbres de blocs — exact

L'état de transfert doit mémoriser les deux composantes étiquetées
$`C_1,C_2`$, la connectivité frontière, le temps du dernier minimum et les
trois tailles de groupe. La mixture de gagnantes de la section 3 s'intègre à
chaque transition. On obtient alors exactement la loi jointe des
$`(X_{v,0},X_{v,1},X_{v,2})`$ et le message $`B_u`$.

### 9.2 Bandes de largeur fixée — certifié

Une matrice de transfert finie traite exactement la géométrie. Les poids
continus et les log-poids peuvent être propagés par arithmétique d'intervalles.
La fonctionnelle $`\mathcal R_u`$ fournit un critère d'arrêt vérifiable et une
erreur explicite sur $`B_u`$, puis sur $`\eta_u`$.

### 9.3 Grille entière — estimations de bord

Il faut relier

```math
m_{v,1}=|E(C_1,S_v)|,
\qquad
m_{v,2}=|E(C_2,S_v)|,
\qquad
m_{v,0}=|E(P_v\setminus C_u,S_v)|
```

à la géométrie multi-bras des composantes sélectionnées par le pivot critique.
Les objectifs minimaux sont des bornes de moments uniformes pour ces trois
quantités, les écarts $`1-\beta_v`$ et $`\ell_v^{-1}`$. Une estimation du seul
volume des composantes ne contrôle pas leurs interfaces et ne suffit donc pas.

### 9.4 Calcul Monte-Carlo à variance réduite — diagnostic, pas preuve

On peut simuler uniquement le squelette géométrique, puis intégrer exactement
les marques par le noyau conditionnel. Cette Rao--Blackwellisation mesure
directement $`\Gamma_p(\mathscr S_u)`$ et évite de resimuler les horloges
marquées. Les sorties doivent enregistrer chaque terme de $`\mathcal R_u`$ et
la contribution de chaque ancêtre ; une moyenne numérique seule ne constitue
pas une preuve de sommabilité.

## 10. Contre-audits obligatoires

1. **Moyenne dans un logarithme.** En général,
   $`\mathbb E\log\Lambda\ne\log\mathbb E\Lambda`$ et
   $`\mathbb E\log\Lambda\le\log\mathbb E\Lambda`$. Un taux nul rend même le
   premier membre égal à $`-\infty`$ sur l'événement correspondant.
2. **Gagnante pondérée.** L'uniformité de la gagnante est fausse hors du cas
   homogène ; le bon poids est $`w_es_{v,e}`$.
3. **Dépendance des groupes.** Les groupes sont indépendants sachant la
   gagnante, mais anticorrélés après sa marginalisation.
4. **Dépendance des ancêtres.** Les marques se factorisent sachant le squelette
   complet ; les géométries $`(m_{v,r},\beta_v)`$ restent fortement corrélées
   le long de la chaîne.
5. **Quatre états nécessaires.** Estimer $`\Lambda_v^{00}`$ seul ne contrôle
   ni les compléments, ni $`J_v`$, ni un coin proche de zéro.
6. **Criticité non héritée.** $`\beta_u\simeq\beta_c`$ n'implique pas
   $`\beta_v\simeq\beta_c`$ pour $`v\succ u`$.
7. **Mauvais biais de paire.** Un nœud uniforme n'a pas la loi du LCA d'une
   paire ; il faut le poids $`N_{u,L}^{\mathrm{far}}`$, ou
   $`|C_{u,1}||C_{u,2}|`$ lorsque toutes les paires sont incluses.
8. **Oracle favorable et portée logique.** Le succès conditionnel des fusions
   les plus favorables n'est pas une condition suffisante de weak recovery.
   En revanche, leur impossibilité donne une impossibilité globale dès que la
   domination HF du fichier 12 est démontrée ; ce lemme ne doit pas rester
   implicite.

## 11. Implémentation et vérification

Le module
[computations/ancestral_lambda_estimation.py](computations/ancestral_lambda_estimation.py)
implémente la course pondérée, les moments exacts, le transport quatre états et
le certificat $`\mathcal R_u`$. Les tests
[computations/test_ancestral_lambda_estimation.py](computations/test_ancestral_lambda_estimation.py)
contre-auditent les formules par une énumération indépendante de toutes les
marques sur de petits buckets, vérifient la réduction homogène et testent les
bornes de contraste sur des réalisations aléatoires reproductibles.
