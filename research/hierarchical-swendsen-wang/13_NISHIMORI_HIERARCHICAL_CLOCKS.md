# Conjecture de Nishimori et horloges hiérarchiques

> **Statut dans le programme.** Cette note est une calibration entropique
> auxiliaire. Elle fournit une identité exacte de face, mais n'est pas une
> route autonome vers la weak recovery. La preuve prioritaire passe par le
> [corridor collapsed](20_COLLAPSED_CORRIDOR_BLACKWELL.md).

Ce fichier répond à une question précise : la dynamique hiérarchique peut-elle
retrouver le nombre triangulaire

```math
p_{\mathrm N}^{(0)}=0.8358057923\ldots
```

de Nishimori--Ohzeki, puis transformer ce calcul en une borne de weak
recovery ? La réponse comporte trois niveaux qu'il faut absolument séparer.

1. **Établi exactement.** L'équation publiée par Nishimori--Ohzeki est
   équivalente à une identité d'entropie conditionnelle à quatre états. Elle
   admet une formulation sans répliques par une course de quatre horloges
   exponentielles et possède une unique racine dans $`(1/2,1)`$.
2. **Compatible exactement avec la dynamique hiérarchique.** Au LCA critique,
   les quatre poids $`q_u^{ab}`$ forment eux aussi une course exponentielle.
   Chacun dépend de tous les $`\Lambda_v^{ab}`$ pour $`v\succeq u`$. Une
   élimination collapsed des orientations fournit la conditionnelle exacte à
   $D$ fixé.
3. **À prouver.** L'égalité autoduale du premier niveau n'identifie pas, à
   elle seule, le seuil de weak recovery. Il faut contrôler, sous une même
   loi, les trois groupes de liens de chaque ancêtre, les quatre
   $`\Lambda_v^{ab}`$, puis établir le lemme favorable HF du fichier 12.

Ainsi, la dynamique réobtient exactement la **constante conjecturée au niveau
face** et fournit un bilan hiérarchique précis à tester. Elle ne donne pas
encore une preuve que cette constante est le vrai seuil.

## 1. Le canal de bruit d'une face triangulaire

Sur une grille triangulaire finie, écrivons

```math
O_{xy}=\Sigma_x\Sigma_y Z_{xy},
\qquad
\mathbb P(Z_{xy}=+1)=p,
```

avec les $`Z_e`$ indépendants. Pour une face triangulaire $`f`$ de bord
$`e_1,e_2,e_3`$, son syndrome observé est

```math
S_f
:=
\prod_{e\in\partial f}O_e
=
Z_{e_1}Z_{e_2}Z_{e_3}.
```

Les spins disparaissent parce que chaque spin de la face apparaît deux fois.
Posons $`q=1-p`$ et

```math
r_p
:=
\mathbb P(S_f=+1)
=p^3+3pq^2
=\frac{1+(2p-1)^3}{2}.
```

Conditionnellement à $`S_f`$, le triplet de bruit prend quatre valeurs :

| syndrome | mot de bruit | probabilité conditionnelle |
|---|---|---:|
| $`S_f=+1`$ | $`(+,+,+)`$ | $`p^3/r_p`$ |
| $`S_f=+1`$ | chacune des trois permutations de $`(+,-,-)`$ | $`pq^2/r_p`$ |
| $`S_f=-1`$ | $`(-,-,-)`$ | $`q^3/(1-r_p)`$ |
| $`S_f=-1`$ | chacune des trois permutations de $`(-,+,+)`$ | $`qp^2/(1-r_p)`$ |

Notons

```math
h_2(x)=-x\log_2x-(1-x)\log_2(1-x).
```

### Proposition 1.1 — identité d'entropie de face

Pour tout $`p\in[0,1]`$,

```math
\boxed{
H(Z_{e_1},Z_{e_2},Z_{e_3}\mid S_f)
=
3h_2(p)-h_2(r_p).
}
```

### Preuve

L'indépendance donne

```math
H(Z_{e_1},Z_{e_2},Z_{e_3})=3h_2(p).
```

Comme $`S_f`$ est une fonction déterministe du triplet et a pour loi
$`\mathrm{Bernoulli}(r_p)`$, la règle de chaîne donne

```math
H(Z_{e_1},Z_{e_2},Z_{e_3})
=H(S_f)+H(Z_{e_1},Z_{e_2},Z_{e_3}\mid S_f),
```

d'où la formule.

## 2. Réduction exacte de l'équation de Nishimori--Ohzeki

L'équation (28) de Nishimori--Ohzeki s'écrit

```math
\begin{aligned}
&2p^2(3-2p)\log p
+2q^2(1+2p)\log q+\log2\\
&\quad=
p(4p^2-6p+3)\log(4p^2-6p+3)\\
&\qquad\quad
+q(4p^2-2p+1)\log(4p^2-2p+1).
\end{aligned}
```

Or

```math
r_p=p(4p^2-6p+3),
\qquad
1-r_p=q(4p^2-2p+1),
```

et

```math
2p^2(3-2p)+r_p=3p,
\qquad
2q^2(1+2p)+(1-r_p)=3q.
```

En faisant passer le membre droit à gauche, on obtient exactement

```math
\log2-3H_{\mathrm e}(p)+H_{\mathrm e}(r_p)=0,
```

où $`H_{\mathrm e}`$ désigne l'entropie de Shannon en nats, appliquée ici à
une loi binaire. Par la proposition 1.1,
cette équation est donc

```math
\boxed{
H(Z_{e_1},Z_{e_2},Z_{e_3}\mid S_f)=1\ \text{bit}.
}
\qquad\text{(NO)}
```

Cette simplification n'est ni une approximation numérique ni une limite de
répliques : c'est une identité algébrique exacte entre l'équation publiée et
une entropie de canal fini.

### Proposition 2.1 — existence et unicité de la racine supérieure

L'équation (NO) possède une unique solution dans $`(1/2,1)`$, donnée par

```math
0.835805792366
<p_{\mathrm N}^{(0)}
<0.835805792368,
\qquad
p_{\mathrm N}^{(0)}=0.835805792367\ldots.
```

### Preuve

Posons $`x=2p-1`$ et

```math
G(p)=3H_{\mathrm e}(p)-H_{\mathrm e}(r_p).
```

Pour $`p\in(1/2,1)`$,

```math
G'(p)
=-6\left(
\mathop{\mathrm{artanh}}(x)
-x^2\mathop{\mathrm{artanh}}(x^3)
\right)<0.
```

En effet,
$`0<x^2<1`$ et
$`0<\mathop{\mathrm{artanh}}(x^3)<\mathop{\mathrm{artanh}}(x)`$.
Enfin,

```math
G(1/2)=2\log2,
\qquad
G(1)=0.
```

La continuité et la stricte décroissance donnent une unique solution de
$`G(p)=\log2`$. La valeur affichée est obtenue par dichotomie dans le script
associé. Pour contre-audit, une évaluation décimale à 80 chiffres donne pour
$`G(p)-\log2`$ respectivement
$`4.0296\times10^{-12}`$ et $`-4.0416\times10^{-12}`$ aux deux bornes
décimales ci-dessus ; la monotonie fournit alors l'encadrement annoncé.

Au point racine, les valeurs de contrôle sont

```math
r_p=0.651469272865\ldots,
```

```math
H(Z\mid S_f=+1)=0.645278785002\ldots\ \text{bit},
```

```math
H(Z\mid S_f=-1)=1.663040455296\ldots\ \text{bit},
```

et leur moyenne pondérée vaut exactement $`1`$ bit à la précision du calcul.

## 3. Formulation exacte par horloges exponentielles

### Lemme 3.1 — course de normalisation

Soit $`\mathcal X`$ fini et soient $`w_x\ge0`$, avec
$`W=\sum_xw_x>0`$. Tirons indépendamment

```math
T_x\sim\mathrm{Exp}(w_x)
```

et posons $`T_\star=\min_xT_x`$ et
$`X_\star=\mathop{\mathrm{argmin}}_xT_x`$. Alors

```math
\mathbb P(X_\star=x)=\frac{w_x}{W},
\qquad
T_\star\sim\mathrm{Exp}(W).
```

En particulier, si $`\pi_x=w_x/W`$,

```math
\mathbb E[-\log\pi_{X_\star}]
=-\sum_x\pi_x\log\pi_x
=H(\pi),
```

et

```math
\mathbb E[-\log T_\star]
=\gamma+\log W,
```

où $`\gamma`$ est la constante d'Euler.

### Preuve

La probabilité de survie commune jusqu'au temps $`t`$ est

```math
\prod_xe^{-w_xt}=e^{-Wt}.
```

La densité que l'horloge $`x`$ gagne à l'instant $`t`$ vaut
$`w_xe^{-Wt}`$. Son intégrale est $`w_x/W`$. La dernière identité suit de
$`WT_\star\sim\mathrm{Exp}(1)`$.

Appliquons ce lemme, conditionnellement à $`S_f=s`$, aux quatre taux

```math
w_z=\mathbb P(Z=z\mid S_f=s).
```

Si $`Z_\star`$ est le gagnant de cette course, alors (NO) devient

```math
\boxed{
\mathbb E[-\log_2\mathbb P(Z_\star\mid S_f)]=1.
}
```

C'est une représentation exacte de l'équation conjecturée par des horloges
exponentielles, sans continuation analytique en un nombre de répliques.
L'identité
$`\mathbb E[-\log T_\star]=\gamma+\log W`$ représente de plus exactement
les `log-sum-exp` et les log-partitions quenched qui apparaissent dans une
condition d'autodualité.

> **Contre-audit.** Ces horloges de normalisation ne sont pas encore les
> horloges d'arêtes de Kruskal. Elles donnent une réalisation exponentielle
> du même calcul fini. Leur identification avec une étape de la dynamique
> hiérarchique exige la construction collapsed ci-dessous.

### Proposition 3.2 — la face comme heat bath hiérarchique collapsed

Soit $`\triangle`$ un triangle à trois sommets et notons
$`\overline\Sigma_\triangle`$ la configuration de ses spins modulo flip
global. Pour une observation $`O_\triangle`$ fixée, l'application

```math
[\sigma]
\longmapsto
z(\sigma)
:=
(O_{xy}\sigma_x\sigma_y)_{xy\in E_\triangle}
```

est une bijection entre les quatre configurations de spins modulo flip et le
coset de syndrome

```math
\mathcal C_{S_\triangle}
=
\left\{z\in\{\pm1\}^3:\prod_ez_e=S_\triangle\right\}.
```

De plus,

```math
\boxed{
\mathbb P(\overline\Sigma_\triangle=[\sigma]\mid O_\triangle)
=
\mathbb P(Z=z(\sigma)\mid S_\triangle).
}
\qquad\text{(F)}
```

Par conséquent,

```math
\mathbb E\left[
H(\overline\Sigma_\triangle\mid O_\triangle)
\right]
=
H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3).
```

Un dendrogramme binaire à trois feuilles possède deux nœuds internes, donc
deux bits d'orientation relative. Le heat bath collapsed qui marginalise
$D$ et rafraîchit **conjointement** ces deux bits est exactement la course à
quatre états de la section 3. Ainsi (NO) est bien une identité sur une mise à
jour hiérarchique complète d'une face.

### Preuve

Pour tout $`[\sigma]`$, le mot $`z(\sigma)`$ a le syndrome observé. Réciproquement,
un mot de ce coset fixe les trois produits $`\sigma_x\sigma_y`$ et donc la
classe de spins. Par Bayes et l'a priori uniforme,

```math
\mathbb P([\sigma]\mid O_\triangle)
\propto
\prod_{e\in E_\triangle}
p^{\mathbf1_{\{z_e(\sigma)=+1\}}}
(1-p)^{\mathbf1_{\{z_e(\sigma)=-1\}}}.
```

La constante de normalisation est $`\mathbb P(S_\triangle)`$, ce qui donne
(F). Le dernier énoncé utilise les deux coordonnées relatives de l'arbre binaire à trois feuilles et le lemme 3.1.

> **Contre-audit.** La course de quatre états d'une face n'est pas le heat
> bath d'un seul LCA. Un nœud interne ne porte qu'un bit relatif ; les quatre
> choix $`(a,b)`$ de ses deux flips contiennent deux copies liées par le flip
> global. Il faut rafraîchir les deux nœuds de la face et marginaliser $D$ pour
> obtenir les quatre classes distinctes de (F).


## 4. Ce que l'identité de face ne démontre pas

Le résultat précédent réobtient exactement l'équation conjecturée au niveau
d'une face. Trois opérations différentes ne doivent toutefois pas être
confondues.

1. **Heat bath collapsed de la face.** On marginalise le dendrogramme et l'on
   rafraîchit conjointement les deux bits relatifs. On retrouve la postérieure
   à quatre états et l'équation (NO).
2. **Heat bath d'un seul LCA.** On conditionne sur les autres orientations et
   sur $`D`$. Les quatre choix de flips ne représentent que deux parités,
   chacune dupliquée par le flip global.
3. **Oracle de fusion critique.** On conditionne en plus par le fait que la
   paire appartient au même arbre et fusionne au voisinage de
   $`\beta_c`$. Cette sélection révèle de l'information sur les spins.

La mesure jointe vérifie exactement

```math
\nu(d\sigma,dD\mid O)
=
\mu(d\sigma\mid O)P(dD\mid\sigma,O).
```

Si l'on intègre complètement $`D`$, on retrouve $`\mu`$ et les facteurs
hiérarchiques disparaissent. Si l'on conserve $`D`$, les
$`\Lambda_v^{ab}`$ interviennent, mais l'expérience devient plus informative.
Il n'existe donc pas de défaut autodual hiérarchique canonique obtenu en
supprimant simplement les facteurs au-delà d'une profondeur $`K`$.

Le [calcul des frontières critiques](14_CRITICAL_COMPONENT_BOUNDARY.md)
donne le test hiérarchique concret qui remplace une troncature autoduale :
calculer les majorités dans les deux groupes affectés de chaque ancêtre, puis
évaluer exactement

```math
q_u^{00}+q_u^{11}
>
q_u^{10}+q_u^{01}.
```

Ce critère, et non l'ancienne notation provisoire $`\widehat\Psi_K`$, est
maintenant calculé sur le cactus dans le fichier 21. Il doit ensuite être
transporté aux bandes triangulaires.

## 5. Audit critique de la constante

L'équation

```math
3h_2(p)-h_2\left(\frac{1+(2p-1)^3}{2}\right)=1
```

et sa racine supérieure sont établies exactement comme objets finis. Leur
identification au point multicritique reste l'ansatz principal de
Nishimori--Ohzeki.

La méthode améliorée par amas d'Ohzeki déplace déjà la première approximation
triangulaire vers $`0.835985\ldots`$. Une hiérarchie complète peut donc corriger
$`p_{\mathrm N}^{(0)}`$ au lieu de le reproduire à chaque profondeur.

Le fichier 14 fournit un contre-audit directement adapté au cas favorable.
Les deux seuils de majorité **de frontière** sont

```math
p_{\mathrm{SW}}=0.673648\ldots,
\qquad
p_{\partial,\mathrm{late}}=0.782432\ldots.
```

Ils sont tous deux strictement inférieurs à la baseline
$`p_{\mathrm{info}}=0.794659\ldots`$ et, a fortiori, à
$`p_{\mathrm N}^{(0)}`$. La condition de majorité est donc déjà très
fortement satisfaite au point de Nishimori ; elle n'y change ni de signe ni
de régime. Plus précisément,

```math
h_c(p_{\mathrm N}^{(0)})
=
0.4968797586\ldots,
\qquad
s_c(p_{\mathrm N}^{(0)})
=
0.7484398793\ldots.
```

Réobtenir la constante conjecturée exige donc une identité sur les
**amplitudes** des quatre poids collapsed, pas un simple vote conforme.

## 6. Statut et prochain calcul

| Énoncé | Statut |
|---|---|
| Équation (28) $`\Longleftrightarrow H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1`$ bit | **Établi exactement** |
| Unicité de la racine dans $`(1/2,1)`$ | **Établie** |
| Course exponentielle à quatre états d'une face collapsed | **Établie** |
| Identification de cette racine au seuil de weak recovery | **Conjecture** |
| Troncature autoduale brute des facteurs ancestraux | **Abandonnée : non canonique** |
| Décomposition exacte de la majorité postcritique | **Établie dans le fichier 14** |
| Canal collapsed sur une chaîne de cactus triangulaires | **Établi exactement dans le fichier 21** |
| Identification de l'équation de Nishimori par ce canal cactus | **Non obtenue : le coefficient reste strictement contractant pour $`(1+q_\triangle)/2<p<1`$** |

L'ordre de travail est maintenant :

1. conserver le calcul de face comme calibration exacte ;
2. utiliser le cactus du fichier 21 comme cas-test exact : il perd toute
   persistance pour chaque $`(1+q_\triangle)/2<p<1`$ et ne sélectionne donc
   pas la constante de Nishimori ;
3. construire le transfert collapsed d'une bande triangulaire de largeur
   deux, où les cycles se chevauchent ;
4. chercher seulement ensuite une identité autoduale sur ce transfert ou
   une limite de bandes.

Le calcul de face reproductible est dans
[`computations/nishimori_hierarchical_entropy.py`](computations/nishimori_hierarchical_entropy.py),
avec son
[`contre-audit à huit états`](computations/test_nishimori_hierarchical_entropy.py).

## Sources primaires

- [Nishimori--Ohzeki, 2006](https://arxiv.org/abs/cond-mat/0601356) :
  dualité triangulaire face-vers-face, condition principale conjecturale et
  équation (28).
- [Ohzeki, 2009](https://arxiv.org/abs/0811.0464) : amélioration par amas et
  flot de renormalisation ; première approximation triangulaire
  $`0.835985\ldots`$.
