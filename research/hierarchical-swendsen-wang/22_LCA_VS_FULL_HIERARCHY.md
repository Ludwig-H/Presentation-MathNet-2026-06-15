# Jusqu'où descendre dans la hiérarchie ?

Cette note compare rigoureusement les deux dynamiques proposées pour une
paire $`i,j`$ dont le LCA $`u`$ fusionne au niveau critique
$`\beta_u=\beta_c(p)`$ :

1. rééchantillonner seulement les orientations des deux enfants de $`u`$ ;
2. poursuivre les rééchantillonnages sur les deux bras jusqu'aux feuilles
   $`i,j`$.

En volume fini, l'égalité de niveaux est comprise au sens de la densité de
Palm du flux de fusions, ou comme limite d'une fenêtre autour de
$`\beta_c`$. L'événement ponctuel a probabilité nulle.

La réponse dépend du sens du mot « mieux ».

- Pour construire l'obstruction de weak recovery la plus forte, il faut
  intégrer **tout le corridor descendant**, idéalement en un heat bath
  collapsed conjoint.
- Pour fabriquer le cas le plus favorable à la conservation de la vérité,
  le LCA seul est plus favorable, car il laisse toutes les relations
  descendantes figées.
- Si l'on impose un seul sweep séquentiel, bottom-up possède une comparaison
  $`L^2`$ avec le LCA seul ; top-down n'en possède pas sans structure
  supplémentaire.

Sur le cactus du fichier 21, ces trois affirmations deviennent des identités
fermées. Le LCA seul conserve une persistance indépendante de la distance,
alors que le corridor complet gagne un facteur strictement contractant par
triangle descendant.

Comme dans le fichier 21, « conserver la vérité » signifie ici préserver la
relation de la réplique de référence dans le couplage annealed de Nishimori.
Ce n'est pas une probabilité quenched de succès d'un estimateur à $`O`$ fixé.

## 1. Trois opérations qu'il ne faut pas confondre

Fixons $`O,D`$ et écrivons

```math
\pi_D(d\sigma)=\nu_O(d\sigma\mid D),
\qquad
f_{ij}(\sigma)=\sigma_i\sigma_j.
\tag{1.1}
```

Soit $`u=\mathrm{LCA}_D(i,j)`$, avec enfants $`C_i`$ et $`C_j`$ contenant
respectivement $`i`$ et $`j`$. Notons $`\mathcal G_v`$ la tribu conservée
par le heat bath du nœud $`v`$ et

```math
P_vg=\mathbb E_{\pi_D}[g\mid\mathcal G_v].
\tag{1.2}
```

Soit $`\mathcal C_{ij}^{\downarrow}`$ l'union de $`u`$ et des deux bras qui
descendent de $`C_i,C_j`$ jusqu'aux feuilles $`i,j`$. Posons

```math
\mathcal A_{ij}^{\downarrow}
=
\bigcap_{v\in\mathcal C_{ij}^{\downarrow}}\mathcal G_v,
\qquad
P_{\downarrow}g
=
\mathbb E_{\pi_D}[g\mid\mathcal A_{ij}^{\downarrow}].
\tag{1.3}
```

Les trois opérations pertinentes sont alors :

| opération | états rééchantillonnés | opérateur |
|---|---|---|
| flip commun des deux enfants | seulement $`(0,0)`$ et $`(1,1)`$ | ne change jamais $`f_{ij}`$ |
| LCA quatre états | $`(0,0),(1,0),(0,1),(1,1)`$ au nœud $`u`$ | $`P_u`$ |
| corridor descendant collapsed | toutes les orientations des deux bras en un bloc | $`P_{\downarrow}`$ |

Le premier choix n'est pas un mécanisme de décorrélation de paire. Dans la
suite, « LCA seul » signifie toujours le heat bath exact à quatre états.

## 2. Théorème général de profondeur

Comme le bloc descendant rééchantillonne toutes les variables du LCA et des
deux bras, on a

```math
\mathcal A_{ij}^{\downarrow}\subseteq\mathcal G_u.
\tag{2.1}
```

### Théorème 2.1 — le corridor complet domine le LCA seul, statut : établi

Pour tout $`g\in L^2(\pi_D)`$,

```math
\boxed{
\|P_{\downarrow}g\|_2^2
\le
\|P_ug\|_2^2.
}
\tag{2.2}
```

Plus précisément,

```math
\boxed{
\|P_ug\|_2^2
=
\|P_{\downarrow}g\|_2^2
+
\|(P_u-P_{\downarrow})g\|_2^2.
}
\tag{2.3}
```

#### Preuve

L'inclusion (2.1) donne les identités de projections conditionnelles

```math
P_{\downarrow}P_u=P_{\downarrow},
\qquad
P_uP_{\downarrow}=P_{\downarrow}.
\tag{2.4}
```

Ainsi $`P_u-P_{\downarrow}`$ est la projection orthogonale sur
$`\mathrm{Ran}(P_u)\cap\mathrm{Ran}(P_{\downarrow})^\perp`$, et

```math
P_ug
=
P_{\downarrow}g+(P_u-P_{\downarrow})g
```

est une décomposition orthogonale. Pythagore donne (2.3), puis (2.2).

### Critère de stricte amélioration

L'inégalité est stricte pour $`g=f_{ij}`$ si et seulement si

```math
(P_u-P_{\downarrow})f_{ij}\ne0.
\tag{2.5}
```

Autrement dit, après le seul heat bath du LCA, il reste de l'information sur
la relation de paire dans les orientations descendantes. C'est précisément
l'information que le corridor complet rééchantillonne.

### Conséquence pour la weak recovery

Les seconds moments conditionnels associés aux deux couplages sont

```math
A_u(i,j;O,D)=\|P_uf_{ij}\|_2^2,
\qquad
A_{\downarrow}(i,j;O,D)=\|P_{\downarrow}f_{ij}\|_2^2.
\tag{2.6}
```

Le théorème donne point par point

```math
A_{\downarrow}(i,j;O,D)\le A_u(i,j;O,D).
\tag{2.7}
```

Pour prouver une impossibilité, le membre de gauche est donc le meilleur
certificat. Pour éprouver le scénario le plus favorable à la récupération,
le membre de droite est le test le plus sévère.

## 3. Pourquoi le flip commun ne sert pas

Si $`a,b\in\{0,1\}`$ indiquent les flips des deux enfants de $`u`$, alors

```math
\frac{f_{ij}(\sigma^{ab})}{f_{ij}(\sigma)}
=
(-1)^{a+b}.
\tag{3.1}
```

Les états pairs $`(0,0)`$ et $`(1,1)`$ conservent la relation ; les états
impairs $`(1,0)`$ et $`(0,1)`$ l'inversent. Restreindre la dynamique aux deux
états pairs donne donc exactement

```math
H(i,j)=1.
\tag{3.2}
```

Il faut utiliser les quatre poids

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
\tag{3.3}
```

La probabilité de retournement de la relation au LCA est

```math
\mathbb P_u(\text{impair}\mid O,D,\sigma)
=
\frac{q_u^{10}+q_u^{01}}
{q_u^{00}+q_u^{10}+q_u^{01}+q_u^{11}}.
\tag{3.4}
```

Le fait que $`\beta_u=\beta_c`$ améliore la qualité du bucket de $`u`$, mais
ne produit aucune atténuation multiplicative automatique par le seul nombre
de niveaux sous $`u`$. Sur la grille, le heat bath du LCA peut encore dépendre
des tailles et des messages ancestraux ; sur le cactus libre, son coefficient
est exactement indépendant de la distance. L'accumulation contrôlable de
distance vient des deux bras descendants.

## 4. Un sweep séquentiel n'est pas le bloc collapsed

Pour un ordre physique $`S=(v_1,\ldots,v_M)`$, l'opérateur de Markov agit sur
les fonctions selon

```math
K_S=P_{v_1}\cdots P_{v_M}.
\tag{4.1}
```

Les projections agissent donc sur $`f_{ij}`$ dans l'ordre inverse du sweep
physique.

### 4.1 Bottom-up : une comparaison garantie

Si le sweep visite les descendants de bas en haut et termine physiquement au
LCA, alors

```math
K_{\rm BU}=Q_{\downarrow}P_u,
\tag{4.2}
```

où $`Q_{\downarrow}`$ est un produit de projections descendantes. Comme tout
produit de contractions est une contraction,

```math
\boxed{
\|K_{\rm BU}f_{ij}\|_2
\le
\|P_uf_{ij}\|_2.
}
\tag{4.3}
```

Ainsi, si l'on impose un seul parcours et si l'on veut une amélioration
au sens large, c'est-à-dire une non-augmentation garantie par rapport au LCA
seul, l'ordre bottom-up est le choix canonique. L'inégalité peut être une
égalité.

### 4.2 Top-down : aucune comparaison universelle en un passage

Pour un sweep qui commence physiquement au LCA,

```math
K_{\rm TD}=P_uQ_{\downarrow}.
\tag{4.4}
```

Les projections ne commutent généralement pas. La seule contractivité ne
permet donc pas de comparer $`\|P_uQ_{\downarrow}f\|_2`$ à
$`\|P_uf\|_2`$.

### Contre-exemple 4.1 — non-commutation, statut : établi exactement

Sur $`\Omega=\{1,2,3,4\}`$ uniforme, prenons les partitions

```math
\mathcal P_A=\{\{1,2\},\{3\},\{4\}\},
\qquad
\mathcal P_B=\{\{1\},\{2,3\},\{4\}\},
```

et les espérances conditionnelles $`P_A,P_B`$ correspondantes. Pour

```math
f=(1,-1,0,0),
```

on a

```math
P_Af=0,
\qquad
P_AP_Bf=\left(\frac14,\frac14,-\frac12,0\right),
```

donc

```math
\|P_Af\|_2^2=0,
\qquad
\|P_AP_Bf\|_2^2=\frac3{32}.
\tag{4.5}
```

Ce contre-exemple ne dit pas que top-down est mauvais sur tout dendrogramme.
Il interdit seulement de déduire sa supériorité de la géométrie abstraite
des projections. Une comparaison top-down exige la loi réelle des messages
ancestraux et descendants.

### 4.3 Le collapsed reste l'enveloppe optimale

Le théorème 2.1 du fichier 20 donne, pour tout sweep $`S`$ utilisant les
mêmes nœuds,

```math
\boxed{
\|P_{\downarrow}f_{ij}\|_2^2
\le
\|K_Sf_{ij}\|_2^2.
}
\tag{4.6}
```

En volume fini, répéter cycliquement tous les heat baths du corridor fait
converger le produit de projections vers la projection sur l'intersection de
leurs espaces :

```math
K_S^m f
\longrightarrow
P_{\downarrow}f
\qquad(m\to\infty).
\tag{4.7}
```

Le heat bath collapsed est donc à la fois l'optimum $`L^2`$ et la limite
d'un balayage local répété. Il évite de choisir artificiellement entre un
seul passage top-down et bottom-up.

## 5. Calcul exact sur le cactus sous LCA critique

Reprenons la chaîne de $`h`$ triangles du fichier 21 et conditionnons, au
sens de la densité de Palm,

```math
q_{\rm LCA}(a_0,a_h)=q.
\tag{5.1}
```

Il existe presque sûrement un unique triangle pivotal au rang $`q`$. Les
$`h-1`$ autres triangles sont seulement conditionnés à avoir connecté leurs
articulations avant $`q`$.

### Théorème 5.1 — LCA seul, statut : établi

Si l'on rééchantillonne uniquement les deux enfants du LCA global, alors

```math
\boxed{
A_h^{\rm LCA\ only}(p,q)
=
\kappa_{\rm flux}(p,q).
}
\tag{5.2}
```

Cette quantité ne dépend pas de $`h`$.

#### Preuve

Le heat bath du LCA ne rééchantillonne que la relation du triangle pivotal.
Son coefficient sous le flux de fusions vaut $`\kappa_{\rm flux}`$ par le
théorème 4.1 du fichier 21. Toutes les relations sur les bras descendants
restent égales à celles de la réplique de référence et contribuent le facteur
$`1`$.

### Théorème 5.2 — tous les descendants, statut : établi

Pour le corridor descendant collapsed complet,

```math
\boxed{
A_h^{\rm full}(p,q)
=
\kappa_{\rm flux}(p,q)
\kappa_{\rm conn}(p,q)^{h-1}.
}
\tag{5.3}
```

Ainsi le gain exact de profondeur vaut

```math
\boxed{
\frac{A_h^{\rm full}(p,q)}
{A_h^{\rm LCA\ only}(p,q)}
=
\kappa_{\rm conn}(p,q)^{h-1}.
}
\tag{5.4}
```

Pour tout $`p<1`$ admissible et $`q>0`$, ce rapport tend
exponentiellement vers zéro.

### Valeurs à $`p=0.8`$ et $`q=q_\triangle`$

On a

```math
\kappa_{\rm flux}=0.791530736866\ldots,
\qquad
\kappa_{\rm conn}=0.886752566857\ldots.
\tag{5.5}
```

La conformité Nishimori après le LCA seul est donc, pour tout $`h`$,

```math
P_h^{\rm LCA\ only}
=
\frac{1+\kappa_{\rm flux}}2
=0.895765368433\ldots.
\tag{5.6}
```

| $`h`$ | second moment LCA seul | rapport full/LCA | second moment full | conformité full |
|---:|---:|---:|---:|---:|
| 2 | $`0.791530736866`$ | $`0.886752566857`$ | $`0.701891912662`$ | $`0.850945956331`$ |
| 5 | $`0.791530736866`$ | $`0.618315049484`$ | $`0.489415366734`$ | $`0.744707683367`$ |
| 10 | $`0.791530736866`$ | $`0.339017477840`$ | $`0.268342754045`$ | $`0.634171377023`$ |
| 20 | $`0.791530736866`$ | $`0.101917000003`$ | $`0.080670438112`$ | $`0.540335219056`$ |
| 40 | $`0.791530736866`$ | $`0.009210765320`$ | $`0.007290603861`$ | $`0.503645301931`$ |

Le LCA critique fournit un bucket de bonne qualité, mais il n'accumule aucune
distance. Toute la perte asymptotique vient des $`h-1`$ relations
descendantes non pivotales.

## 6. Les ancêtres stricts au-dessus du LCA

Un flip d'un ancêtre strict de $`u`$ retourne simultanément $`i`$ et $`j`$ ;
il n'apparaît donc pas dans l'identité déterministe de parité. Il influence
néanmoins les probabilités des décisions par les facteurs
$`\Lambda_v(\sigma^{ab})`$ de (3.3).

Il faut distinguer deux opérations.

1. **Évaluer** tous les facteurs ancestraux : c'est obligatoire même pour le
   LCA seul et pour le corridor descendant.
2. **Rééchantillonner** aussi leurs orientations dans un bloc collapsed :
   cela réduit encore la tribu conservée et ne peut qu'améliorer le
   certificat $`L^2`$.

Si $`P_{\rm root}`$ désigne ce bloc étendu jusqu'à la racine, alors

```math
\|P_{\rm root}f_{ij}\|_2^2
\le
\|P_{\downarrow}f_{ij}\|_2^2
\le
\|P_uf_{ij}\|_2^2.
\tag{6.1}
```

Étendre le bloc à tout l'arbre finit cependant par reformuler directement
le problème de la postérieure complète. Le corridor descendant est le
meilleur compromis actuel entre contraction et calculabilité.

## 7. Conséquence pour la stratégie $`p=0.8`$

Deux optimisations indépendantes interviennent.

1. **Géométrie favorable.** Parmi les fusions postcritiques, placer le LCA à
   $`\beta_c`$ maximise la persistance sur le cactus.
2. **Couplage contractant.** À cette géométrie fixée, rééchantillonner tout le
   corridor minimise la persistance parmi les programmes utilisant ces
   coordonnées.

Il n'y a pas de contradiction : on majore les géométries défavorables par la
géométrie critique la plus informative, puis on choisit dans cette géométrie
le couplage postérieur le plus utile au théorème d'obstruction.

L'ordre (2.7) ne suffit toutefois pas à prouver une limite nulle. Les fusions
descendantes ont des niveaux plus précoces que $`\beta_c`$ et peuvent être
presque déterministes. Dans un corridor factorisé de coefficients
$`\kappa_r`$, la condition exacte est

```math
\prod_r\kappa_r\longrightarrow0
\quad\Longleftrightarrow\quad
\sum_r-\log\kappa_r\longrightarrow+\infty.
\tag{7.1}
```

Sur la grille, cette somme doit être remplacée par une contraction de blocs
ou un rayon spectral du transfert avec état de frontière. « Descendre » est
donc le bon ordre de dynamique, pas encore un théorème de perte.

La cible correcte sur la grille devient

```math
\mathbb E\left[
\|P_{\downarrow}f_{I_LJ_L}\|_2^2
\middle|
\mathrm{Palm}(\beta_{I_LJ_L}=\beta_c)
\right]
\longrightarrow0.
\tag{7.2}
```

Le LCA seul ne peut pas produire cette limite par un simple argument de
distance : sur le cactus, son second moment reste exactement
$`\kappa_{\rm flux}`$. Le prochain certificat de bande doit donc représenter
les deux bras complets et non seulement le bloc pivotal.

## 8. Audit et contre-audit

| affirmation | statut | raison ou limite |
|---|---|---|
| Retourner simultanément les deux enfants change la relation $`i,j`$ | Faux | les états $`(0,0)`$ et $`(1,1)`$ sont pairs |
| Le heat bath du LCA utilise quatre états | Établi | seuls $`(1,0)`$ et $`(0,1)`$ inversent la relation |
| Le corridor collapsed est au moins aussi contractant que le LCA seul | Établi point par point | projections conditionnelles imbriquées |
| Descendre donne une amélioration stricte dans tous les environnements | Faux | égalité possible dans (2.5) |
| Un sweep bottom-up complet est au plus persistant que le LCA seul | Établi en $`L^2`$ | le LCA est le dernier update physique |
| Un sweep top-down complet est toujours meilleur que le LCA seul | Non démontré et faux pour des projections générales | contre-exemple (4.5) |
| Répéter les sweeps converge vers le collapsed | Établi en volume fini | produit cyclique de projections |
| Sur le cactus, le LCA seul perd la corrélation avec la distance | Faux | coefficient constant $`\kappa_{\rm flux}`$ |
| Sur le cactus, le corridor complet perd la corrélation | Établi | facteur $`\kappa_{\rm conn}^{h-1}`$ |
| Descendre suffit toujours à faire tendre le second moment vers zéro | Faux | il faut une atténuation cumulée divergente |
| La même factorisation vaut sur la grille | Ouvert | cycles chevauchants et état de frontière |
| La conformité calculée est un succès quenched d'estimation | Faux | conformité de la réplique de référence dans le couplage Nishimori |

## Conclusion opérationnelle

Pour utiliser au mieux la dynamique hiérarchique dans la preuve de weak
recovery :

1. conditionner le LCA de $`i,j`$ à fusionner à $`\beta_c`$ ;
2. conserver les quatre états exacts au LCA et tous les
   $`\Lambda_v`$ ancestraux ;
3. descendre sur les deux bras jusqu'aux feuilles ;
4. utiliser le heat bath collapsed de ce corridor pour la preuve ;
5. utiliser un sweep bottom-up en un passage, puis top-down comme
   contre-audit séquentiel ;
6. certifier sur une bande de largeur deux que le transfert complet conserve
   un rayon spectral strictement inférieur à un.

Le LCA critique est le meilleur **point de départ**. Il n'est pas une
dynamique suffisante : la distance ne devient exploitable qu'en utilisant
toute la hiérarchie descendante.
