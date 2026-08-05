# Le problème central : seuil de fusion d'une paire lointaine et chaîne ancestrale des $`\Lambda_v`$

Cette note est la **préface mathématique** du dossier ancestral et la
formulation canonique du problème que poursuit tout le programme
hiérarchique. Elle définit le problème, fixe une table de notations unique,
rassemble les énoncés exacts déjà démontrés dans la trilogie
[08](08_ANCESTRAL_LAMBDA_CHAIN.md) /
[10](10_ANCESTRAL_LAMBDA_ESTIMATION.md) /
[14](14_CRITICAL_COMPONENT_BOUNDARY.md), ajoute trois identités
structurelles compactes, puis isole le verrou unique restant et les routes
de résolution. Chaque énoncé porte son statut : **établi**,
**établi sous références standard**, **conditionnel**, **réfuté** ou
**ouvert**. Cette version intègre un contre-audit adversarial complet
(cinq lentilles indépendantes) ; les énoncés qui y ont survécu sont
signalés, les autres ont été restreints.

## 1. Énoncé du problème

Le cadre est celui de la [mesure jointe exacte](../01_MATHEMATICAL_FRAMEWORK.md) :
conditionnellement à une réplique postérieure $`\sigma`$, chaque arête
satisfaite reçoit une horloge $`\xi_e\sim\mathrm{Exp}(|W_e|)`$, les autres
$`\xi_e=+\infty`$ ; la filtration $`\Pi_t`$ des composantes du graphe
$`\{\xi_e\le t\}`$ définit le dendrogramme de partitions **non marqué** $D$,
et

```math
\nu_O(\sigma\mid D)
\propto
\mu_0(\sigma)
\prod_{v\in D}
\Lambda_v(\sigma)\,e^{(1-\beta_v)\Lambda_v(\sigma)}.
```

Soient $`i_L,j_L`$ deux sommets avec $`d(i_L,j_L)\to\infty`$ et

```math
\beta_{ij}
:=
\inf\{t:\ i\leftrightarrow j\ \text{dans }\Pi_t\},
\qquad
u=\mathrm{LCA}_D(i,j),
\qquad
\beta_u=\beta_{ij}.
```

Le problème central se décompose en deux volets.

> **P1 — seuil de fusion.** La bonne variable pour la weak recovery
> pairwise est le niveau $`\beta_{ij}`$ où les clusters de $i$ et $j$
> fusionnent. Fixer les conventions de conditionnement licites autour de
> $`\beta_c(p)`$, démontrer la localisation conditionnelle de
> $`\beta_{ij}`$ au seuil, et délimiter ce que cette restriction
> favorable peut — et ne peut pas — donner pour le score global.
>
> **P2 — chaîne ancestrale.** Sur l'événement de fusion, la persistance
> de la relation $`\sigma_i\sigma_j`$ sous le heat bath du nœud $u$ est
> gouvernée par $`L_u=B_u+\ell_u`$, où le message ancestral $`B_u`$
> agrège les quatre taux $`\Lambda_v^{ab}`$ de **tous** les ancêtres
> $`v\succ u`$. Estimer ces $`\Lambda_v^{ab}`$ — c'est-à-dire la loi du
> squelette groupé $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succ u}`$ sous
> la loi de la paire — est **le** problème d'estimation du programme.

La portée est fixée par la
[réduction pairwise](../03_HIERARCHICAL_WEAK_RECOVERY.md) : pour interdire
la weak recovery, il suffit d'annuler asymptotiquement le second moment
pairwise. La borne LCA $`Q_n\le H_n^{\mathrm{LCA}}`$
([03 §4](../03_HIERARCHICAL_WEAK_RECOVERY.md)) le relie à
$`\eta_u=\tanh^2(L_u/2)`$ ; elle est **établie conditionnellement à A1**,
la formalisation complète de $`\nu_O`$ et de ses conventions de censure
sur l'exhaustion finie
([07](../../diagnostics/07_CRITICAL_BAND_CRITERION.md),
[02](../02_CHAPTER_11_BASELINE.md)).

## 2. Table de notations

Cette table est la référence du dossier ancestral ; elle lève les
collisions de symboles entre les notes historiques.

| symbole | définition | référence |
|---|---|---|
| $`u_p=\log\frac p{1-p}`$ | taux des horloges (GSBM homogène, $`p>1/2`$) | [04 §1](../04_TRIANGULAR_GSBM.md) |
| $`q_p(t)=p(1-e^{-u_pt})`$ | paramètre de percolation au temps $t$ | [04](../04_TRIANGULAR_GSBM.md) |
| $`q_c=2\sin(\pi/18)`$, $`\beta_c(p)=q_p^{-1}(q_c)`$ | seuil et temps critiques (réseau triangulaire) | [04](../04_TRIANGULAR_GSBM.md) |
| $`p_{\mathrm{SW}}=\frac{1+q_c}2=0.673648\ldots`$ | domaine d'existence de $`\beta_c\le1`$ | [09](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md) |
| $`s_p(t)=\mathrm{logistic}(u_p(1-t))`$ | qualité d'une marque fermée au temps $t$ | [14 §3](14_CRITICAL_COMPONENT_BOUNDARY.md) |
| $`h_p(t)=2s_p(t)-1=\tanh\frac{u_p(1-t)}2`$ | biais de marque ; $`s_c,h_c`$ à $`t=\beta_c`$ | [14 §3](14_CRITICAL_COMPONENT_BOUNDARY.md) |
| $`E_v^{(r)},\ T_{v,r},\ \lambda_{v,r}`$ | groupes d'arêtes d'une coupe ancestrale, poids totaux et satisfaits | [08 §3](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`m_{v,r}=|E_v^{(r)}|`$, $`m_v=\sum_r m_{v,r}`$ | tailles des groupes | [08 §6](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`X_{v,r}=2\lambda_{v,r}-T_{v,r}`$ | déséquilibre signé du groupe $r$ | [08 §3](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`\Delta_{v,r}=T_{v,r}-2\lambda_{v,r}=-X_{v,r}`$ | déséquilibre orienté (convention de [10 §7](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) | [10 §7](10_ANCESTRAL_LAMBDA_ESTIMATION.md) |
| $`\Lambda_v^{ab}`$ | taux de la coupe $`E_v`$ dans l'état flippé $`\sigma^{ab}`$ | [08 §3](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`\lambda_v^{\min}=\min_{ab}\Lambda_v^{ab}`$ | plus petit des quatre taux (noté $`\ell_v`$ dans 10 §7) | [10 §7](10_ANCESTRAL_LAMBDA_ESTIMATION.md) |
| $`\phi_v(x)=\log x+(1-\beta_v)x`$ | log-facteur d'un nœud | [08 §4](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`\ell_u=\log\frac{\Lambda_u}{T_u-\Lambda_u}+(1-\beta_u)(2\Lambda_u-T_u)`$ | log-odds **local** du nœud $u$ | [01 §4](../01_MATHEMATICAL_FRAMEWORK.md) (décomposition), [16, lemme 4.2](../../archive/oracles/16_FLIP_PROBABILITIES_DESCENDANT_PATH.md) (symbole) |
| $`B_u`$, $`L_u=B_u+\ell_u`$ | message ancestral, log-odds de parité | [01 §4](../01_MATHEMATICAL_FRAMEWORK.md) |
| $`\eta_u=\tanh^2(L_u/2)`$ | fiabilité pairwise | [03](../03_HIERARCHICAL_WEAK_RECOVERY.md) |
| $`h_1,h_2,J`$ | coordonnées de Walsh de $`(\Phi_u^{ab})`$ | [08 §5](08_ANCESTRAL_LAMBDA_CHAIN.md) |
| $`H_u`$ | **nombre d'ancêtres stricts** de $u$ (noté $h$ dans 08, $H$ dans 10) | ici |
| $`\mathscr S_u`$ | squelette ancestral non marqué (noté $`\mathscr D`$ dans 08 §6 et 14 §6) | [10 §1](10_ANCESTRAL_LAMBDA_ESTIMATION.md) |
| $`R_u(I)`$ | certificat de troncature (noté $`\mathcal R_u(I)`$ dans 10) | [10 §7](10_ANCESTRAL_LAMBDA_ESTIMATION.md) |
| $`N_{u,L}^{\mathrm{far}}`$ | nombre de paires lointaines séparées par $u$ | [10 §1](10_ANCESTRAL_LAMBDA_ESTIMATION.md) |

Symboles historiques à ne pas confondre : $`h_p(t)`$ (biais de marque),
$`h_1,h_2`$ (champs de Walsh), $`H_u`$ (profondeur ancestrale) ; et
$`\ell_u`$ (log-odds local du nœud) contre $`\lambda_v^{\min}`$ (plus
petit des quatre taux d'un ancêtre, $`\ell_v`$ dans 10 §7). La constante
$`(2+q_c)/3=0.782432\ldots`$ est notée uniformément
$`p_{\partial,\mathrm{late}}`$.

## 3. Pourquoi seuls les ancêtres comptent

Trois faits déterministes, valables à poids hétérogènes, délimitent
exactement ce qui entre dans le heat bath de $u$. Les lemmes (b) et (c)
ont été vérifiés par simulation exhaustive (300 dendrogrammes, tous les
nœuds internes, zéro violation) lors du contre-audit.

**(a) Annulation des descendants et des arêtes internes — établi.**
Pour tout descendant strict $w$ de $u$ et toute arête interne à $`C_1`$ ou
$`C_2`$, la satisfaction est invariante sous les flips $(a,b)$ ; les
facteurs correspondants se simplifient dans les quatre poids
([14, lemme 2.1](14_CRITICAL_COMPONENT_BOUNDARY.md),
[16, lemme 6.1](../../archive/oracles/16_FLIP_PROBABILITIES_DESCENDANT_PATH.md)).
Le heat bath exact de $u$ ne contient que $u$, ses ancêtres et l'a priori.

**(b) Partition des arêtes par leur LCA — établi.** Chaque arête
potentielle $`e=\{x,y\}`$ dont les extrémités appartiennent à la même
racine finale de $`\Pi_1`$ appartient à exactement une coupe du
dendrogramme : $`e\in E_{v(e)}`$ avec $`v(e)=\mathrm{LCA}_D(x,y)`$. Une
arête entre deux racines finales distinctes n'appartient à aucune coupe.

*Preuve.* $`E_v`$ est l'ensemble des arêtes joignant les deux enfants de
$v$ ; $`e\in E_v`$ si et seulement si les enfants de $v$ séparent $x$ de
$y$ et $v$ les réunit, c'est-à-dire $`v=\mathrm{LCA}_D(x,y)`$, qui est
unique. $\square$

**(c) Comptabilité de frontière — établi.** Fixons $`u:C_u=C_1\dot\cup C_2`$
et $`r\in\{1,2\}`$. Au temps $`\beta_u`$, toute arête quittant $`C_u`$ est
fermée. Les groupes ancestraux $`E_v^{(r)}`$, $`v\succ u`$, forment une
partition exacte de l'ensemble

```math
\partial^{\to}C_r
:=
\bigl\{
e=\{x,y\}:\ x\in C_r,\ y\notin C_u,\
y\ \text{dans la même racine finale que}\ C_u
\bigr\},
```

chaque arête étant absorbée au niveau $`\beta_{v(e)}\in(\beta_u,1]`$ de
l'unique ancêtre $`v(e)`$ où la composante de son extrémité extérieure
rejoint celle de $`C_u`$. En particulier

```math
\boxed{
\sum_{v\succ u}m_{v,r}
=
|\partial^{\to}C_r|,
\qquad r\in\{1,2\}.
}
```

*Preuve.* Par (b), $`e\in E_{v(e)}`$ avec $`v(e)=\mathrm{LCA}_D(x,y)`$.
Comme $`\beta_{v(e)}>\beta_u`$, la composante de $x$ juste avant
$`\beta_{v(e)}`$ contient $`C_u`$ tout entier : l'enfant de $`v(e)`$
contenant $x$ est donc $`P_{v(e)}`$ et $y$ est dans $`S_{v(e)}`$, d'où
$`e\in E_{v(e)}^{(r)}`$ puisque l'extrémité intérieure est dans $`C_r`$.
Réciproquement toute arête de $`E_v^{(r)}`$ est de cette forme. Les
arêtes de frontière dont l'extrémité extérieure vit dans une autre racine
finale (censure) n'apparaissent dans aucun facteur. $\square$

La conséquence est une reformulation géométrique du volet P2 : estimer les
$`\Lambda_v^{ab}`$ revient à estimer **comment la frontière fermée de
$`C_1`$ et $`C_2`$ se répartit entre les niveaux de fusion successifs**
$`\beta_v\in(\beta_u,1]`$ — chaque arête non gagnante absorbée au niveau
$`\beta_v`$ portant une marque de qualité $`s_p(\beta_v)`$, la gagnante
étant conforme ([théorème de course, 10 Th. 3.1](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) —
plus le contrôle des tailles $`m_{v,0}`$ du niveau commun.

## 4. P1 — le seuil de fusion et sa localisation critique

**Hypothèse permanente de cette section :** $`p\ge p_{\mathrm{SW}}`$, de
sorte que $`\beta_c(p)\le1`$ et que $`\Pi_{\beta_c}`$ appartient au
dendrogramme censuré à $`1`$
([10 §1](10_ANCESTRAL_LAMBDA_ESTIMATION.md),
[14 §1](14_CRITICAL_COMPONENT_BOUNDARY.md)).

### 4.1 Loi du seuil et identité déterministe — établi

Pour $`t\le1`$,

```math
\mathbb P(\beta_{ij}\le t)
=
\tau_{ij}(q_p(t)),
\qquad\text{et}\qquad
\{i\leftrightarrow j\ \text{dans }\Pi_{\beta_c}\}
=
\{\beta_{ij}\le\beta_c\}.
```

### 4.2 Localisation par la gauche — esquissé sous intrants standard

Pour une **paire déterministe** $`(i_L,j_L)`$ avec
$`d_L(i_L,j_L)\ge\rho L`$, sur le tore triangulaire ou une exhaustion à
conditions de bord contrôlées (choix à fixer,
[14 §10.1](14_CRITICAL_COMPONENT_BOUNDARY.md)) : pour tout $`\delta>0`$,

```math
\mathbb P\bigl(\beta_{ij}\le\beta_c-\delta
\bigm|\ \beta_{ij}\le\beta_c\bigr)
=
\frac{\tau_{ij}(q_p(\beta_c-\delta))}{\tau_{ij}(q_c)}
\longrightarrow0,
```

par décroissance exponentielle sous-critique au numérateur et minoration
polynomiale de type box-crossing au dénominateur
([14, prop. 5.1](14_CRITICAL_COMPONENT_BOUNDARY.md) — esquisse de
preuve ; intrants : sharpness de Duminil-Copin–Tassion, box-crossing
isoradial de Grimmett–Manolescu). Pour des extrémités uniformes
conditionnées par la distance, le membre de droite devient un quotient
d'espérances, et la conclusion $`\to0`$ subsiste par uniformité de la
borne exponentielle. Donc, **sous ce conditionnement**,
$`\beta_{ij}\to\beta_c`$ en probabilité : la restriction « fusion très
proche du seuil » est automatique dès que l'on conditionne la paire
lointaine à être connectée au niveau critique.

Deux limites essentielles, à ne jamais omettre :

- **la masse de l'événement s'annule** :
  $`\rho^c_{L,\rho}=\mathbb P(d\ge\rho L,\ \beta_{ij}\le\beta_c)\to0`$
  ([14, (4.4)/(8.1)](14_CRITICAL_COMPONENT_BOUNDARY.md)). La localisation
  est donc **purement conditionnelle** : aucune borne sur $`Q_n`$ n'en
  découle sans réinsérer le facteur de masse
  ([14 §8](14_CRITICAL_COMPONENT_BOUNDARY.md),
  [09, cor. 2.2](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md)) ;
- l'échelle **quantitative** de la fenêtre (exposants near-critical)
  n'est pas requise pour la localisation qualitative, mais elle l'est
  pour transporter les identités critiques exactes de §5.5 (voir la mise
  en garde qui y figure) ; elle reste ouverte pour la percolation par
  liens.

### 4.3 Conventions de conditionnement — quatre expériences distinctes

À volume fini les temps sont continus : $`\mathbb P(\beta_{ij}=\beta_c)=0`$.
Quatre formulations licites, à ne pas mélanger :

1. **Conditionnement gauche** $`\{\beta_{ij}\le\beta_c\}`$
   ([15 §1](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)) :
   probabilité strictement positive à tout $L$, aucune fenêtre auxiliaire
   nécessaire ; la localisation 4.2 y est un théorème (esquissé). C'est
   la convention de l'expérience favorable P1 et des hypothèses CUT/ANC
   ([15 §5](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)).
   *(Une fenêtre $`[\beta_c-\varepsilon,\beta_c]`$ avec l'ordre
   $`L\to\infty`$ puis $`\varepsilon\downarrow0`$ définit la même limite
   conditionnelle : par 4.2,
   $`\mathbb P(\beta_{ij}<\beta_c-\varepsilon\mid\text{gauche})\to0`$,
   d'où une équivalence en variation totale — valable pour les
   fonctionnelles **bornées** ($`\eta_u`$, probabilités), pas pour des
   espérances non uniformément intégrables telles que $`B_u`$ ou
   $`R_u`$.)*
2. **Palm de flux** : pondération des nœuds par
   $`N_{u,L}^{\mathrm{far}}`$ et fenêtre
   $`I_\varepsilon=[\beta_c-\varepsilon,\beta_c+\varepsilon]`$, avec
   l'ordre $`L\to\infty`$ puis $`\varepsilon\downarrow0`$
   ([10 §1](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) ; conventions de flux
   pré-saut $`m\,N_\rho`$ contre nœud réalisé $`N_\rho`$ dans
   [25](../25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md).
   **Attention** : à $`\varepsilon`$ fixé, la fenêtre bilatérale est
   asymptotiquement dominée par sa moitié **surcritique**
   ($`\mathbb P(\beta_{ij}\le\beta_c)\to0`$ polynomialement, contre une
   masse $`\theta(q_p(\beta_c+\varepsilon))^2>0`$ à droite) : avec cet
   ordre de limites, la convention 2 coïncide avec la bande 3, pas avec
   le conditionnement gauche 1. Une expérience réellement critique exige
   une fenêtre $`\varepsilon_L`$ d'échelle proche-critique, à annoncer.
3. **Bande surcritique** $`(\beta_c,\beta_c+\delta]`$ : sert au bilan de
   masse
   $`\mathbb E[\eta_u\mathbf 1_{\{\beta_c<\beta_{ij}\le\beta_c+\delta\}}]
   \le S_n(\beta_c+\delta)-S_n(\beta_c)`$
   ([09 §2](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md)).
4. **Désintégration régulière au niveau $\beta$** (par bucket) : la
   convention sous laquelle sont établis le théorème de course
   ([10, th. 3.1](10_ANCESTRAL_LAMBDA_ESTIMATION.md), §5.3 ci-dessous)
   et le calcul local critique de
   [09 §1](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md) repris en
   §5.5. Ce n'est pas une loi de paire mais une désintégration par
   niveau de fusion.

Enfin, deux réductions exactes délimitent le rôle de l'événement
favorable : si $`\beta_{ij}>1`$ (racines distinctes), la persistance d'un
sweep complet est **exactement nulle**
([19, lemme 2.1 et cor. 2.2](../19_FAVORABLE_SWEEP_PROJECTIONS.md)) ; et
l'oracle favorable **ne domine pas** les fusions postcritiques en
général — la domination HF est réfutée en multiport
([29](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md)). Le cas
favorable est donc un benchmark conditionnel exact, pas un raccourci de
preuve globale.

## 5. P2 — structure exacte de la chaîne au-dessus de $u$

**Hypothèse permanente de cette section : a priori i.i.d. uniforme.**
Pour un $`\mu_0`$ général, le coefficient de Walsh $`xy`$ de
$`\log\mu_0(\sigma^{ab})`$ s'ajoute à $J$ (et ses coefficients $`x,y`$ à
$`h_1,h_2`$) ; le certificat de majorité 5.4(3) ne s'applique alors plus
([14, th. 7.1, encadré Portée](14_CRITICAL_COMPONENT_BOUNDARY.md)).

### 5.1 Les quatre taux et la réduction de Walsh — établi

Pour chaque ancêtre $`v\succ u`$
([08 §3–§5](08_ANCESTRAL_LAMBDA_CHAIN.md)) :

```math
\Lambda_v^{ab}
=
\tfrac12\bigl[T_v+X_{v,0}+(-1)^aX_{v,1}+(-1)^bX_{v,2}\bigr],
\qquad
T_v=T_{v,0}+T_{v,1}+T_{v,2},
```

et, lorsque les quatre taux sont strictement positifs, avec
$`\Phi_u^{ab}=\log\mu_0(\sigma^{ab})+\sum_{v\succ u}\phi_v(\Lambda_v^{ab})
=C_u+h_1x+h_2y+Jxy`$ :

```math
\boxed{
B_u
=
2J+\log\cosh(h_1+h_2)-\log\cosh(h_1-h_2),
}
\qquad
J=\sum_{v\succ u}J_v,\quad
J_v=\tfrac14\log\frac{\Lambda_v^{00}\Lambda_v^{11}}
{\Lambda_v^{01}\Lambda_v^{10}}.
```

Le terme linéaire $`(1-\beta_v)\Lambda_v^{ab}`$, affine en $(a,b)$, ne
contribue jamais à $J$ : le couplage direct des deux flips provient du
seul préfacteur de course $`\log\Lambda_v^{ab}`$.

### 5.2 Trois identités structurelles compactes — établi

Les deux premières sont implicites dans la forme affine de 5.1 ; la
troisième est nouvelle sous cette forme fermée. Elles sont démontrées
ci-dessous et vérifiées numériquement, ainsi qu'en arithmétique
rationnelle exacte, par
[computations/ancestral_walsh_identities.py](../../computations/ancestral_walsh_identities.py).

**(i) Sommes appariées.** Pour chaque ancêtre,

```math
\Lambda_v^{00}+\Lambda_v^{11}
=
\Lambda_v^{01}+\Lambda_v^{10}
=
T_v+X_{v,0}.
```

**(ii) Produit apparié et signe du couplage.**

```math
\boxed{
\Lambda_v^{00}\Lambda_v^{11}-\Lambda_v^{01}\Lambda_v^{10}
=
-\,X_{v,1}X_{v,2},
}
\qquad\text{donc}\qquad
\mathrm{sgn}(J_v)=-\,\mathrm{sgn}(X_{v,1}X_{v,2}).
```

*Preuve.* $`4\Lambda_v^{00}\Lambda_v^{11}=(T_v+X_{v,0})^2-(X_{v,1}+X_{v,2})^2`$
et $`4\Lambda_v^{01}\Lambda_v^{10}=(T_v+X_{v,0})^2-(X_{v,1}-X_{v,2})^2`$ ;
la différence vaut $`-4X_{v,1}X_{v,2}`$. $\square$

Deux groupes alignés dans le même sens créent donc un couplage direct
$`J_v<0`$ défavorable à la parité paire. Corollaire de signe (établi,
c'est l'étape « cône de Walsh » de la preuve de
[14, th. 7.1](14_CRITICAL_COMPONENT_BOUNDARY.md), qui n'utilise pas
l'hypothèse locale (7.7)) : **si $`X_{v,1},X_{v,2}\ge0`$ pour tout
ancêtre $`v\succ u`$, alors $`h_1,h_2\ge0`$,
$`\log\cosh(h_1+h_2)-\log\cosh(h_1-h_2)\ge0`$ et $`B_u\ge0`$** — les
quatre coefficients de Walsh de chaque facteur sont positifs et le cône
est stable par produit. Hors de cet événement d'alignement simultané,
aucun contrôle du signe de $`B_u`$ n'est démontré : un seul ancêtre
aligné dans une chaîne non alignée peut laisser $`B_u<0`$.

**(iii) Ancêtre terminal isolé muet.** Si la chaîne se réduit à un
**unique** ancêtre de niveau $`\beta_v=1`$ (a priori uniforme), alors
$`B_u=0`$ : $`\phi_v=\log`$ et, par (i),
$`\mathrm{LSE}(\log\Lambda_v^{00},\log\Lambda_v^{11})
=\log(T_v+X_{v,0})
=\mathrm{LSE}(\log\Lambda_v^{01},\log\Lambda_v^{10})`$.
**La muteté n'est pas additive** : au sein d'une chaîne, un ancêtre à
$`\beta_v=1`$ garde $`J_v\ne0`$ et $`h_{v,1},h_{v,2}\ne0`$ par la seule
courbure du $`\log`$ (identité (ii)), et peut même renverser le signe du
message — contre-exemple exact du contre-audit :
$`v_1=(\beta=0.25,\,T=6,\,X_0=0,\,X_1=-0.5,\,X_2=-3)`$ seul donne
$`B_u=+0.44`$, la chaîne $`\{v_1,v_2\}`$ avec
$`v_2=(\beta=1,\,T=6,\,X_0=0,\,X_1=3,\,X_2=2.5)`$ donne $`B_u=-1.77`$.
Seule l'inclinaison exponentielle $`(1-\beta_v)`$ s'éteint à la censure.

### 5.3 Loi exacte des marques sachant le squelette — établi

**Cas homogène** ([08 §6](08_ANCESTRAL_LAMBDA_CHAIN.md)) :
conditionnellement au squelette non marqué $`\mathscr S_u`$ (au sens de
la désintégration par niveau, convention 4.3(4)), les buckets sont
indépendants, une catégorie gagnante $`G_v`$ est tirée avec
$`\mathbb P(G_v=r)=m_{v,r}/m_v`$, puis

```math
K_{v,r}
\stackrel{d}{=}
\mathbf 1_{\{G_v=r\}}
+
\mathrm{Bin}\bigl(m_{v,r}-\mathbf 1_{\{G_v=r\}},\,s_p(\beta_v)\bigr),
\qquad
K_u\stackrel{d}{=}1+\mathrm{Bin}(m_u-1,s_p(\beta_u)).
```

**Cas pondéré** ([10, th. 3.1](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) : la
gagnante est tirée avec $`\rho_{v,e}\propto w_es_{v,e}`$ et les comptes
sont des Poisson-binomiales pondérées — la loi homogène ci-dessus ne s'y
réduit pas, et remplacer $`\rho_{v,e}`$ par l'uniforme biaise déjà la
moyenne. Transformée génératrice jointe :
[10 §3](10_ANCESTRAL_LAMBDA_ESTIMATION.md) ; moyennes, covariances
(anticorrélation des groupes après oubli de la gagnante) et concentration
simultanée : [10 §4–§6](10_ANCESTRAL_LAMBDA_ESTIMATION.md). Dans le cas
homogène, **tout le biais restant porte sur la géométrie**
$`(m_{v,r},\beta_v)_{v\succ u}`$ ; dans le cas pondéré, les modules
$`|W_e|`$ y entrent aussi.

### 5.4 Contrôles certifiés du message — établi

1. **Troncature** : sur l'événement
   $`\{\lambda_v^{\min}>0\ \forall v\in I\}`$ (hypothèse permanente de
   [10 §7](10_ANCESTRAL_LAMBDA_ESTIMATION.md) ; sa défaillance est
   exactement le verrou G3),
   $`|B_u-B_u^{(-I)}|\le R_u(I)`$ avec
   $`R_u(I)=\sum_{v\in I}\bigl[(|\Delta_{v,1}|+|\Delta_{v,2}|)
   (1-\beta_v+1/\lambda_v^{\min})+|\Delta_{v,1}\Delta_{v,2}|/(2(\lambda_v^{\min})^2)\bigr]`$ ;
   transport à la fiabilité par la constante $`2/(3\sqrt3)`$.
2. **Sandwich des comptes** : tout ancêtre strict vérifie
   $`1/2\le s_p(\beta_v)<s_p(\beta_u)`$ **sans condition**
   ($`\beta_v>\beta_u`$, $`s_p`$ décroissante,
   [14, (6.8)](14_CRITICAL_COMPONENT_BOUNDARY.md)) ; le sandwich se
   resserre en $`s_c=s_p(\beta_c)`$ — couplage
   $`\mathrm{Bin}(\cdot,1/2)\le K\le\mathrm{Bin}(\cdot,s_c)`$ —
   **seulement si** $`\beta_u\ge\beta_c`$
   ([09 §7, prop. 7.1](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md)).
   Sous la localisation gauche 4.2, un ancêtre peut vivre dans
   $`(\beta_u,\beta_c]`$ avec $`h_p(\beta_v)\ge h_c`$ ; le sandwich vaut
   alors avec le paramètre $`s_p(\beta_u)`$
   ([08, préambule](08_ANCESTRAL_LAMBDA_CHAIN.md),
   [14 §6](14_CRITICAL_COMPONENT_BOUNDARY.md)).
3. **Certificat de majorité hiérarchique** : **sous a priori i.i.d.
   uniforme** (hypothèse essentielle), si $`2\Lambda_u-T_u\ge0`$ et
   $`X_{v,1},X_{v,2}\ge0`$ pour tout $`v\succ u`$, alors
   $`q_u^{00}+q_u^{11}\ge q_u^{10}+q_u^{01}`$, strictement si la majorité
   locale est stricte ([14, th. 7.1](14_CRITICAL_COMPONENT_BOUNDARY.md)) ;
   probabilité conditionnelle **exacte** de l'événement du certificat :
   un produit de $`2H_u`$ queues binomiales
   ([14, cor. 7.2](14_CRITICAL_COMPONENT_BOUNDARY.md)), qui tend vers
   zéro dès que $`H_u\to\infty`$ hors du régime (7.17)
   $`\min_{v,r}m_{v,r}h_p(\beta_v)^2\gg\log H_u`$.

### 5.5 Réponse au nœud critique — établi dans l'oracle local

Dans l'oracle $`B_u=0`$ et sous la désintégration au niveau $`\beta_u`$
(convention 4.3(4)), avec
$`\ell_{m,K}(\beta_u)=\log\frac{K}{m-K}+u_p(1-\beta_u)(2K-m)`$ et
$`K=1+\mathrm{Bin}(m-1,s_p(\beta_u))`$ : la probabilité de parité
conforme moyenne vaut $`(1+\Gamma_m(\beta_u))/2`$
([15, prop. 4.1](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)). Au
niveau critique $`\beta_u=\beta_c`$ : valeur $`1/2+1/(2m)`$ au bord
$`p=p_{\mathrm{SW}}`$ ([15, cor. 4.2](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)),
et, **pour tout $`p>p_{\mathrm{SW}}`$ fixé**, déficit

```math
1-\bar P_m^c\sim C_{m\bmod2}(p)\,m^{-1/2}e^{-mI_c(p)},
\qquad I_c=\log\cosh(a_c/2)
```

([15, th. 4.3](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)) — au
bord $`p_{\mathrm{SW}}`$, $`I_c=0`$ et la série $`C_r`$ diverge :
l'équivalent y est faux, la valeur exacte étant $`1/2+1/(2m)`$.
**Mise en garde de transport** : sous les conventions 4.3(1)–(2),
$`\beta_u\ne\beta_c`$ ; substituer $`(s_c,a_c)`$ aux valeurs réalisées
$`(s_p(\beta_u),u_p(1-\beta_u))`$ commet une erreur extensive
$`\asymp u_p\,h\,m\,|\beta_c-\beta_u|`$, négligeable seulement si
$`|\beta_c-\beta_u|=o(1/m_u)`$ — l'échelle proche-critique que la
localisation qualitative 4.2 ne fournit pas. La réponse hiérarchique
complète est
$`\bar P^{\mathrm{hier}}=\tfrac12[1+\mathbb E^\star\tanh^2((\ell_{M_L,K_L}(\beta_u)+B_L)/2)]`$
([15 §7](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)) : tout le
reste du problème est la loi jointe de $`(\beta_u,M_L,B_L)`$ sous la
convention choisie.

## 6. Le verrou unique et ce qu'il n'est pas

### 6.1 Formulation du verrou — ouvert

Établir ([10 §8](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) :

- **G1** : la convergence (ou un contrôle uniforme) de la loi du squelette
  groupé $`(\beta_u,m_u;(\beta_v,m_{v,0},m_{v,1},m_{v,2})_{v\succ u})`$
  sous la Palm de flux 4.3(2) — de façon équivalente, par la comptabilité
  de frontière 3(c), la loi de la répartition de
  $`\partial^{\to}C_1,\partial^{\to}C_2`$ entre les niveaux de fusion,
  plus celle des niveaux communs $`m_{v,0}`$ ;
- **G2** : la sommabilité de la queue $`R_u`$ des ancêtres profonds
  (sous la mesure étendue aux marques,
  [10 §8](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) ;
- **G3** : le contrôle des coins $`\lambda_v^{\min}\approx0`$ ;
- **CUT / ANC** ([15 §5](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md),
  sous le conditionnement gauche 4.3(1)) : $`M_L\to\infty`$ et
  $`B_L/M_L\to0`$ pour la version limite de 5.5.

C'est un problème de **percolation near-critical géométrique** : les
marques sont fermées (5.3), les fonctionnelles sont certifiées (5.4), il
ne manque que la loi de la géométrie sous le biais de la paire — chaque
famille d'hypothèses étant attachée à sa convention de 4.3.

### 6.2 Ce que le verrou n'est pas — réfuté

- Remplacer les $`\Lambda_v`$ réalisés par leur version criticalisée :
  **réfuté** par le contre-exemple multiport exact
  ([29 §2](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md)).
- Poser $`B_u=0`$ comme approximation conservatrice : ce n'est **pas une
  marginalisation** ([08, préambule](08_ANCESTRAL_LAMBDA_CHAIN.md) ;
  [15 §6, contre-audit 3](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md)),
  et l'oracle local retombe sur $`p_{\mathrm{SW}}`$
  ([09 §5, cor. 5.2](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md)).
- Factoriser les décisions de heat bath le long du corridor : réfuté par
  le contre-audit abstrait minimal de
  [16 §7](../../archive/oracles/16_FLIP_PROBABILITIES_DESCENDANT_PATH.md)
  ($`\mathbb E\prod(-1)^{A_r}\ne\prod\mathbb E(-1)^{A_r}`$) ; l'oracle
  factorisé PATH-FAC de [16 §8](../../archive/oracles/16_FLIP_PROBABILITIES_DESCENDANT_PATH.md)
  reste licite comme benchmark mais n'est ni majorant ni minorant.
- Supposer $`|B_u|\le B_0`$ uniformément : le no-go du potentiel non borné
  ([28 §4](../../archive/roadmaps/28_FIRST_CORRIDOR_P0805_RESULTS.md))
  l'interdit sans preuve de screening.

## 7. Conséquence stratégique : où la décorrélation peut-elle naître ?

**7.1 Signe du message ancestral : ce qui est établi, ce qui est
heuristique.** Sur l'événement du certificat 5.4(3) — majorités
groupées $`X_{v,1},X_{v,2}\ge0`$ pour **tous** les ancêtres, a priori
uniforme — on a $`B_u\ge0`$ (corollaire de 5.2). La probabilité de cet
événement est exactement le produit des $`2H_u`$ queues binomiales de
[14, cor. 7.2](14_CRITICAL_COMPONENT_BOUNDARY.md) : elle tend vers
zéro pour une chaîne longue à petits groupes, et **rien n'est démontré
hors de cet événement** — l'annulation ancestrale reste une
possibilité ouverte (hypothèse ANC,
[15 §5–§6](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md) ;
[08 §10](08_ANCESTRAL_LAMBDA_CHAIN.md)). La forme quantitative exacte
est le **test de non-contraction** (lemme 6.3 de la
[note 12 archivée](../../archive/roadmaps/12_FAVORABLE_HIERARCHICAL_REDUCTION.md),
établi) :

```math
\eta_u
\ge
\tanh^2\left(
\frac{\bigl||\ell_u|-|B_u|\bigr|}{2}
\right),
\qquad\text{et}\qquad
\Gamma^{\mathrm{fav}}\to0
\ \Longrightarrow\
B_u+\ell_u\to0
\ \text{en probabilité,}
```

via $`\mathbb P^\star(|L_u|>x)\le\Gamma^{\mathrm{fav}}/\tanh^2(x/2)`$.
Si la coupe critique est grande ($`\ell_u`$ divergent), une
contraction favorable exigerait la compensation exacte
$`B_u\simeq-\ell_u`$ à la même échelle — pas seulement
$`B_u=O(1)`$. **Heuristique de travail (non démontrée)** : le biais
conforme des marques rend cette compensation improbable, et
l'annulation ancestrale n'est pas la source attendue de
décorrélation ; le rôle de P2 dans une borne supérieure est de
certifier G2/G3/ANC.

**7.2 L'atténuation vient des bras descendants.** Sur le cactus ($h$ =
nombre de triangles, pas la profondeur ancestrale $`H_u`$), le LCA
seul donne $`A_h^{\mathrm{LCA\ only}}=\kappa_{\mathrm{flux}}`$,
indépendant de la distance
([22, th. 5.1](../../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md)),
tandis que le corridor complet donne
$`A_h^{\mathrm{full}}=\kappa_{\mathrm{flux}}\kappa_{\mathrm{conn}}^{h-1}\to0`$
([22, th. 5.2](../../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md) ;
c'est le $`A_h^{\mathrm{LCA}}`$ de
[21, (4.3)](../../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md)).

**7.3 Lien avec la cible répliquée — deux événements complémentaires.**
Avec les indicatrices de
[41, (1.2)](../../active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md) :
le conditionnement gauche est exactement l'événement **same-block**,
$`\{\beta_{ij}\le\beta_c\}=\{A_D(i)=A_D(j)\}=\{s_D=1\}`$ — le terme
diagonal que [41, (3.3)](../../active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md)
borne par $`\sqrt{S_L^c}\to0`$ et **élimine**. La cible
$`\mathcal D_L^\times`$ vit sur l'événement complémentaire
**cross-block** $`g_D(1-s_D)=1`$, c'est-à-dire
$`\{i,j\in R_\star\}\cap\{\beta_c<\beta_{ij}\le1\}`$ : toute la bande
postcritique jusqu'à la censure, **sans** localisation près de
$`\beta_c`$ (fausse pour la connectivité brute,
[08 §1](08_ANCESTRAL_LAMBDA_CHAIN.md) ; masse contrôlée seulement par
$`S_n(\beta_c+\delta)-S_n(\beta_c)`$, 4.3(3)). L'expérience favorable
P1 et la cible répliquée sont donc **complémentaires**, pas
identiques : sur la bande cross-block, le sandwich 5.4(2) s'applique
avec $`\beta_u>\beta_c`$ et toute la machinerie P2 (qui est
indifférente au niveau) reste l'outillage pertinent. L'état récursif
minimal de [41 §7](../../active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md)
contient les données dont dérivent les quatre poids de P2 (rangs
réels, buckets, incidences), mais y ajoute les ports du séparateur et
les lois de bord, et interdit toute réduction scalaire non démontrée.

Routes de travail, par ordre de fermeture réaliste :

| route | objet | statut |
|---|---|---|
| cactus : coefficients collapsed $`\kappa_{\mathrm{flux}},\kappa_{\mathrm{conn}}`$ | contraction exacte du corridor | **fermée** ([21, th. 4.1](../../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md), [22, th. 5.1–5.2](../../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md)) |
| cactus : loi complète de $`(h_1,h_2,J)`$ et de $`B_u`$ | propagation exacte de la chaîne ancestrale | ouverte, outillée ([08 §9](08_ANCESTRAL_LAMBDA_CHAIN.md), [10 §9.1](10_ANCESTRAL_LAMBDA_ESTIMATION.md)) |
| bandes de largeur fixée | matrice de transfert certifiée, état = connectivité + $`(h_1,h_2,J)`$ par intervalles | ouverte, outillée par 5.4 |
| grille : G1 near-critical | répartition de $`\partial^{\to}C_r`$ sur les niveaux $`(\beta_v)`$ sous la Palm 4.3(2) | ouverte (verrou principal) |
| grille : cible répliquée | $`\mathcal D_L^\times\to0`$ via TRI1-o/TRI2, sur la bande cross-block | ouverte ([41](../../active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md)) |

## 8. Statut synthétique

| énoncé | statut |
|---|---|
| loi jointe, quatre poids, $`L_u=B_u+\ell_u`$ | établi ([01](../01_MATHEMATICAL_FRAMEWORK.md)) |
| $`Q_n\le H_n^{\mathrm{LCA}}`$ | établi **conditionnellement à A1** ([03](../03_HIERARCHICAL_WEAK_RECOVERY.md), [07](../../diagnostics/07_CRITICAL_BAND_CRITERION.md)) |
| seuls $u$ et ses ancêtres entrent dans le heat bath | établi (§3a) |
| partition des arêtes par leur LCA, comptabilité de frontière | établi (§3b–c) |
| localisation $`\beta_{ij}\to\beta_c`$ conditionnelle ($`p\ge p_{\mathrm{SW}}`$) | esquissé sous intrants standard (RSW/sharpness ; tore/bord à fixer) |
| masse de l'événement favorable $`\to0`$ | établi ([14 (8.1)](14_CRITICAL_COMPONENT_BOUNDARY.md)) |
| quatre taux, Walsh, $`B_u=2J+\log\cosh(h_1{+}h_2)-\log\cosh(h_1{-}h_2)`$ (a priori uniforme) | établi (§5.1) |
| sommes appariées, produit apparié $`=-X_1X_2`$ ; corollaire de signe | établi (§5.2) |
| ancêtre terminal **isolé** muet ; non additif en chaîne | établi (§5.2iii) |
| course conditionnelle des marques (homogène + pondérée), moments, concentration | établi (§5.3) |
| troncature $`R_u`$ (sur $`\lambda_v^{\min}>0`$), sandwich, certificat de majorité (a priori uniforme) | établi (§5.4) |
| oracle local critique : moyenne, bord, équivalent aigu ($`p>p_{\mathrm{SW}}`$) | établi (§5.5) |
| test de non-contraction ($`\eta_u\ge\tanh^2(||\ell_u|-|B_u||/2)`$) | établi (§7.1) |
| loi du squelette groupé sous Palm (G1), queue (G2), coins (G3), CUT/ANC | **ouvert** |
| « l'annulation ancestrale n'est pas la source de décorrélation » | **heuristique** (§7.1) |
| criticalisation uniforme, $`B_u=0`$, PATH-FAC, $`|B_u|\le B_0`$ gratuit | **réfuté / interdit** (§6.2) |
| nouveau seuil de weak recovery par cette voie | **à prouver** |

Aucun seuil nouveau n'est revendiqué. La borne rigoureuse du dépôt reste
$`p_{\mathrm{WR}}\ge0.809439`$ par la
[voie non hiérarchique](../../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).
