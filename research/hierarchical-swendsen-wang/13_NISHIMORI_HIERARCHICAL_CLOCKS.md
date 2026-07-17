# Conjecture de Nishimori et horloges hiérarchiques

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
   élimination collapsed des orientations fournit une décomposition
   entropique exacte du dendrogramme conditionné.
3. **À prouver.** L'égalité autoduale du premier niveau n'identifie pas, à
   elle seule, le seuil de weak recovery. Il faut encore relier une limite de
   blocs hiérarchiques autoduals à la fiabilité des paires lointaines dont le
   LCA est critique, contrôler l'information révélée par $D$, puis établir le
   lemme favorable HF du fichier 12.

Ainsi, la dynamique réobtient exactement la **constante conjecturée au niveau
face**, et propose une hiérarchie naturelle de corrections. Elle ne donne pas
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
\tag{NO}
```

Cette simplification n'est ni une approximation numérique ni une limite de
répliques : c'est une identité algébrique exacte entre l'équation publiée et
une entropie de canal fini.

### Proposition 2.1 — existence et unicité de la racine supérieure

L'équation (NO) possède une unique solution dans $`(1/2,1)`$, donnée par

```math
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
$`G(p)=\log2`$. La valeur affichée est obtenue par dichotomie certifiée dans
le script associé.

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
\tag{F}
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
(F). Le dernier énoncé utilise les deux bits du lemme 5.1 et le lemme 3.1.

> **Contre-audit.** La course de quatre états d'une face n'est pas le heat
> bath d'un seul LCA. Un nœud interne ne porte qu'un bit relatif ; les quatre
> choix $`(a,b)`$ de ses deux flips contiennent deux copies liées par le flip
> global. Il faut rafraîchir les deux nœuds de la face et marginaliser $D$ pour
> obtenir les quatre classes distinctes de (F).

## 4. La course quatre états au-dessus du LCA critique

Fixons le dendrogramme de partitions non marqué $`D`$ et un LCA
$`u:C_u=C_1\mathbin{\dot\cup}C_2`$. Pour
$`a,b\in\{0,1\}`$, retournons $`C_1`$ si $`a=1`$ et $`C_2`$ si $`b=1`$.
Les poids exacts de la dynamique sont

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
```

Posons

```math
\pi_u^{ab}
=
\frac{q_u^{ab}}{\sum_{c,d}q_u^{cd}}.
```

### Proposition 4.1 — horloge hiérarchique locale exacte

Conditionnellement à toutes les orientations extérieures au mouvement de
$u$, le heat bath des deux fils est le gagnant de quatre horloges
indépendantes de taux $`q_u^{ab}`$. Son entropie vaut

```math
H_u^{(4)}
=-\sum_{a,b}\pi_u^{ab}\log\pi_u^{ab}.
```

La parité relative $`R_u=a\mathbin{\mathsf{xor}}b`$ a pour log-rapport

```math
L_u
=
\log
\frac{q_u^{00}+q_u^{11}}
{q_u^{01}+q_u^{10}},
```

et son entropie est

```math
H(R_u\mid D,O,\text{orientations extérieures})
=h_2\left(\frac{1}{1+e^{-L_u}}\right).
```

### Justification

Le premier énoncé est le lemme 3.1 appliqué aux poids du heat bath. Le second
résulte de la sommation des deux orientations absolues dans chaque classe de
parité. Surtout, chaque poids contient

```math
\prod_{v\succ u}
\Lambda_v^{ab}e^{(1-\beta_v)\Lambda_v^{ab}}.
```

Il est donc impossible de retrouver la bonne course en ne gardant que
$`\Lambda_u`$ ou en posant le message ancestral $`B_u`$ à zéro. Les formules
affines exactes du fichier 08 calculent les quatre taux de chaque ancêtre.

## 5. Entropie collapsed du dendrogramme

Considérons d'abord un arbre binaire plein $`D_C`$ dont les feuilles sont les
sommets d'une composante $`C`$. Fixons un ordre des enfants et une feuille de
référence dans chaque sous-arbre. À chaque nœud interne $u$, associons le bit
$`R_u`$ qui vaut $`0`$ si les spins des feuilles de référence de ses deux fils
sont égaux et $`1`$ sinon. C'est le sens précis de l'orientation relative des
deux fils.

### Lemme 5.1 — coordonnées relatives de l'arbre

L'application

```math
\{\pm1\}^{C}/\{\pm1\}
\longrightarrow
\{0,1\}^{|C|-1},
\qquad
\sigma\longmapsto(R_u)_{u\in D_C^{\mathrm{int}}},
```

est une bijection.

### Preuve

Un arbre binaire plein à $`|C|`$ feuilles possède $`|C|-1`$ nœuds internes.
Fixons le spin d'une feuille de référence. En parcourant l'arbre, chaque bit
relatif détermine l'orientation d'un fils à partir de celle de l'autre ; tous
les spins sont donc déterminés. Inversement, une configuration fixe chaque
orientation relative. Changer tous les spins ne change aucun bit.

Ordonnons les nœuds internes $`u_1,\ldots,u_{|C|-1}`$ et définissons le
log-rapport **collapsed**

```math
L_k^{\mathrm{coll}}
=
\log
\frac{
\mathbb P(R_{u_k}=0\mid O,D,R_{u_1},\ldots,R_{u_{k-1}})
}{
\mathbb P(R_{u_k}=1\mid O,D,R_{u_1},\ldots,R_{u_{k-1}})
}.
```

La règle de chaîne donne l'identité exacte

```math
\boxed{
H(\Sigma_C/\{\pm1\}\mid O,D)
=
\sum_{k=1}^{|C|-1}
\mathbb E\left[
h_2\left(\frac{1}{1+e^{-L_k^{\mathrm{coll}}}}\right)
\right].
}
\tag{CH}
```

Le calcul de $`L_k^{\mathrm{coll}}`$ somme les orientations non révélées par
`log-sum-exp`. Le lemme 3.1 permet de réaliser chacune de ces sommes par une
course exponentielle. Pour un flip relatif au nœud $`u_k`$, les facteurs
strictement sous $`u_k`$ sont invariants, tandis que le facteur de $`u_k`$ et
tous ceux de ses ancêtres changent. Les quatre
$`\Lambda_v^{ab}`$ pour $`v\succeq u_k`$ sont donc exactement les données
nécessaires à (CH).

Cette identité justifie mathématiquement une dynamique hiérarchique
**collapsed**, dans laquelle on élimine successivement les orientations le
long de $`u\leadsto\mathrm{racine}`$ au lieu de les figer.

### Contre-audit : le dendrogramme révèle de l'information

L'entropie pertinente pour la weak recovery conditionne sur $`O`$, pas sur
un dendrogramme auxiliaire conservé. Toujours exactement,

```math
\boxed{
H(\Sigma/\{\pm1\}\mid O)
=
H(\Sigma/\{\pm1\}\mid O,D)
+I(\Sigma;D\mid O).
}
\tag{D}
```

Le terme $`I(\Sigma;D\mid O)`$ est positif en général, car les taux des
horloges dépendent de la configuration. L'omettre ferait artificiellement
paraître la posterior plus informative. Deux voies seulement sont sûres :

1. borner explicitement cette fuite d'information ;
2. analyser le noyau marginal qui rééchantillonne $D$ puis l'oublie.

## 6. Pourquoi le bilan d'un bit ne prouve pas la weak recovery

Sur un tore triangulaire à $`n`$ sommets, les $`n`$ faces montantes ont des
bords disjoints et partitionnent les $`3n`$ arêtes. Leurs syndromes sont donc
indépendants et la proposition 1.1 donne

```math
\boxed{
H(Z_E\mid(S_f)_{f\in F_\triangle})
=
n\bigl(3h_2(p)-h_2(r_p)\bigr).
}
```

Au point $`p_{\mathrm N}^{(0)}`$, ce membre vaut $`n`$ bits. Un arbre
hiérarchique couvrant $`n`$ sommets porte $`n-1`$ bits d'orientation relative.
Il existe donc un appariement entropique asymptotique très suggestif : un bit
résiduel de bruit de face par bit de spin relatif.

Cet appariement n'est pas un critère de reconstruction. Sur un arbre, et déjà
sur un chemin, chaque observation d'arête transmet une quantité d'information
strictement positive ; l'information mutuelle totale est extensive pour tout
$`p>1/2`$, alors que la corrélation entre deux sommets à distance croissante
peut tendre vers zéro. Une entropie totale ne dit pas à quelles échelles
l'information est transportée.

Sur la grille triangulaire, trois difficultés supplémentaires subsistent :

- les syndromes des faces descendantes couplent les blocs de faces montantes ;
- la contrainte de cycle et les interactions effectives créées par
  marginalisation ne sont pas décrites par une face isolée ;
- la weak recovery exige une masse non négligeable de paires lointaines
  informatives, donc un contrôle du LCA critique et non du seul budget global.

Le nombre $`p_{\mathrm N}^{(0)}`$ est ainsi une calibration autoduale exacte,
pas encore une conséquence de la dynamique sur les paires lointaines.

## 7. Hiérarchie autoduale guidée par les $`\Lambda_v`$

La bonne extension consiste à faire croître simultanément le bloc de dualité
et la profondeur ancestrale au-dessus d'un LCA critique. Il faut ici
distinguer le bloc de face, qui est exact, de son plongement dans la loi d'une
paire critique, qui reste à construire.

### Défaut exact au niveau face

En nats, posons

```math
\Psi_0^{\mathrm{face}}(p)
:=
\log2-H_{\mathrm e}(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3).
```

Les sections 2 et 3 établissent

```math
\Psi_0^{\mathrm{face}}(p)=0
\quad\Longleftrightarrow\quad
p=p_{\mathrm N}^{(0)}
```

sur $`(1/2,1)`$. Par la proposition 3.2, c'est le défaut entropique du heat
bath conjoint des deux bits relatifs d'une face après marginalisation de
$D$. C'est le sens exact, et limité, dans lequel le niveau zéro de la
dynamique retrouve la conjecture.

### Candidat ancestral conditionné par la paire critique

Fixons $`u=v_0\prec v_1\prec\cdots`$ et, pour $`K\ge0`$, conservons le nœud
$u$, ses $K$ premiers ancêtres et tous les nœuds internes du bloc triangulaire
retenu. Notons cet ensemble $`\mathcal A_K`$. Après fixation d'une condition
de bord, notons $`\mathcal X_K`$ l'ensemble fini des orientations compatibles
et $`\mathcal D_K^{\mathrm c}`$ l'ensemble des dendrogrammes internes dont le
LCA distingué appartient à la fenêtre critique. Pour $`D_K`$ fixé,
définissons le poids positif primal

```math
W_K^{\mathrm{pr}}(x,D_K)
=
\mu_{0,K}(x)
\prod_{v\in\mathcal A_K}
\Lambda_v(x)e^{(1-\beta_v)\Lambda_v(x)}.
```

Si $`\rho_K`$ désigne la mesure de référence sur les temps et les squelettes,
le poids où le dendrogramme interne est réellement oublié est

```math
\overline W_K^{\mathrm{pr}}(x)
=
\int_{\mathcal D_K^{\mathrm c}}
W_K^{\mathrm{pr}}(x,D_K)\,\rho_K(dD_K).
```

Notons $`\mathcal X_K^{\mathrm{pr}}`$ le sous-ensemble imposé par la condition
principale. La quantité collapsed principale est

```math
Z_K^{\mathrm{pr}}
=\sum_{x\in\mathcal X_K^{\mathrm{pr}}}\overline W_K^{\mathrm{pr}}(x).
```

Une transformation triangulaire face-vers-face, appliquée avant de relâcher
la condition principale, définit le bloc dual et $`Z_K^{\mathrm{du}}`$. Soit
$`a_K`$ la constante déterministe de normalisation imposée par cette
transformation. Le candidat de travail est

```math
\widehat\Psi_K(p)
:=
\mathbb E_p\left[
\log Z_K^{\mathrm{pr}}-
\log Z_K^{\mathrm{du}}-a_K
\ \middle|\
u\text{ est un LCA critique favorable}
\right].
\tag{SD-K}
```

La condition « LCA critique » doit être une fenêtre en volume fini, puis une
mesure de Palm dans la limite. L'intégration précédente est indispensable :
une espérance de $`\log Z(D_K)`$ ne peut pas remplacer le logarithme du poids
marginalisé sans terme correctif. Par le lemme 3.1, chaque log-partition de
(SD-K) se représente exactement par le minimum d'une famille d'horloges
exponentielles. Cette formulation conserve tous les termes non linéaires
$`\log\Lambda_v^{ab}`$ sans employer de répliques.

### Niveau $`K=0`$

Pour une **face isolée**, après marginalisation complète de son dendrogramme,
la proposition 3.2 et la condition principale de Nishimori--Ohzeki donnent
exactement

```math
\Psi_0^{\mathrm{face}}(p)=0
\quad\Longleftrightarrow\quad
3h_2(p)-h_2(r_p)=1,
```

donc redonne $`p_{\mathrm N}^{(0)}=0.835805792367\ldots`$.

En revanche, l'égalité

```math
\widehat\Psi_0=\Psi_0^{\mathrm{face}}
```

sous la loi biaisée d'un LCA critique n'est **pas** automatique : le bucket
de $u$, sa condition de bord et la fuite par $D$ diffèrent d'une face isolée.
Cette identification fait partie de NH1.

### Niveaux $`K\ge1`$

Les niveaux suivants ne doivent pas être obtenus en remplaçant
$`\Lambda_v`$ par sa moyenne. Il faut :

1. garder les quatre valeurs $`\Lambda_v^{00},\Lambda_v^{01},
   \Lambda_v^{10},\Lambda_v^{11}`$ pour chaque ancêtre ;
2. effectuer les sommes collapsed avant l'espérance quenched ;
3. conserver les interactions multi-corps produites par ces sommes ;
4. dualiser un poids de bloc positif avec condition de bord, et non une
   transformée de Fourier naïve d'une loi normalisée qui pourrait changer de
   signe.

Cette hiérarchie est cohérente avec la méthode de blocs améliorée d'Ohzeki :
la condition principale sur une face est son niveau zéro, tandis que des
sommes partielles sur des amas croissants suivent davantage le flot de
renormalisation. La première correction publiée sur la grille triangulaire
donne environ $`0.835985`$, et non $`0.8358058`$. Cette différence est un
signal utile : une hiérarchie complète peut **corriger** la conjecture
originale au lieu de la confirmer exactement.

## 8. Lemme de pont nécessaire vers la weak recovery

Définissons en parallèle la fiabilité collapsed de profondeur $K$ de la paire
critique,

```math
\Gamma_K^{\mathrm{fav}}(p)
=
\mathbb E_p^\star\left[
\tanh^2\left(\frac{L_{u,K}^{\mathrm{coll}}}{2}\right)
\right],
```

où $`L_{u,K}^{\mathrm{coll}}`$ marginalise exactement les orientations des
$K$ premiers ancêtres à partir de leurs quatre taux. La troncature du fichier
12 contrôle le passage à $`K=\infty`$ dès que la queue ancestrale est
sommable.

Pour transformer (SD-K) en seuil de weak recovery, il faut démontrer les
quatre énoncés suivants.

### NH1 — bloc critique canonique

Sous le biais d'une paire lointaine du même arbre dont le LCA est à
$`\beta_c`$, les $K$ premiers buckets ancestraux convergent vers un bloc
aléatoire primal auquel la dualité triangulaire s'applique avec une condition
de bord contrôlée. Au niveau minimal, cette construction doit en outre
identifier $`\widehat\Psi_0`$ à $`\Psi_0^{\mathrm{face}}`$, ou quantifier
explicitement leur défaut.

### NH2 — limite autoduale

Les $`\widehat\Psi_K`$ convergent localement uniformément vers une fonction
$`\widehat\Psi_\infty`$ ayant un unique zéro $`p_{\mathrm H}`$. Les limites
$`K\to\infty`$, taille du volume $`L\to\infty`$ et moyenne quenched doivent
être justifiées dans cet ordre ou par une domination commune.

### NH3 — pont dualité--fiabilité

Le signe de $`\widehat\Psi_\infty(p)`$ doit impliquer une contraction quantitative de
$`\Gamma_\infty^{\mathrm{fav}}(p)`$, après oubli de $D$. Une égalité de libres
énergies sans cette implication ne contrôle aucune corrélation à deux points.

### NH4 — retour non oracle

Le lemme HF du fichier 12 et le contrôle de $`I(\Sigma;D\mid O)`$ transportent
la contraction de l'expérience critique favorable vers toutes les paires
postcritiques réelles. Alors seulement le théorème 3.1 du fichier 12 donne
l'absence de weak recovery.

Si NH1--NH4 sont établis et si $`p_{\mathrm H}=p_{\mathrm N}^{(0)}`$, la
dynamique hiérarchique fournit une nouvelle dérivation du seuil conjecturé.
Si les zéros $`p_K`$ se déplacent avec $K$, elle fournit plutôt une suite
d'approximants contrôlés du véritable seuil.

## 9. Verdict et ordre de travail

Le résultat immédiatement exploitable est le suivant.

| Énoncé | Statut |
|---|---|
| Équation (28) $`\Longleftrightarrow H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1`$ | **Établi exactement** |
| Unicité et valeur $`0.835805792367\ldots`$ | **Établi pour l'équation conjecturée** |
| Réalisation à quatre horloges de l'entropie de face | **Établi exactement** |
| Course locale $`q_u^{ab}`$ avec tous les $`\Lambda_v^{ab}`$ | **Établi, volume fini** |
| Décomposition collapsed (CH) conditionnée par $`D`$ | **Établi, volume fini** |
| Identité $`\Psi_0^{\mathrm{face}}=0`$ pour le heat bath collapsed d'une face isolée | **Établi exactement** |
| Identification de $`\widehat\Psi_0`$ sous le biais LCA à $`\Psi_0^{\mathrm{face}}`$ | **À prouver : NH1** |
| Construction canonique de $`\widehat\Psi_K`$ depuis un LCA critique | **À formaliser : NH1** |
| Convergence et autodualité à profondeur infinie | **À prouver : NH2** |
| Passage de l'autodualité à la fiabilité d'une paire lointaine | **À prouver : NH3** |
| Seuil de weak recovery égal à $`0.8358058\ldots`$ | **Conjecture** |

L'ordre de calcul recommandé est :

1. vérifier $`\Psi_0^{\mathrm{face}}`$ par l'énumération à huit mots déjà fournie ;
2. construire $`\widehat\Psi_0`$, puis $`\widehat\Psi_1`$, sur un cactus de deux niveaux avec les quatre
   $`\Lambda_v^{ab}`$ exacts ;
3. comparer son zéro à $`0.8358058`$ et à $`0.835985`$ ;
4. calculer simultanément $`\Gamma_1^{\mathrm{fav}}`$ afin de tester NH3 ;
5. passer aux bandes triangulaires par matrices de transfert et arithmétique
   d'intervalles ;
6. n'extrapoler au plan qu'après contrôle de la queue ancestrale et du terme
   de fuite (D).

Le calcul reproductible est dans
[`computations/nishimori_hierarchical_entropy.py`](computations/nishimori_hierarchical_entropy.py),
avec un contre-audit indépendant par énumération des huit mots dans
[`computations/test_nishimori_hierarchical_entropy.py`](computations/test_nishimori_hierarchical_entropy.py).

## Sources primaires

- [Nishimori--Ohzeki, 2006](https://arxiv.org/abs/cond-mat/0601356) : dualité triangulaire face-vers-face, condition principale conjecturale et équation (28).
- [Ohzeki, 2009](https://arxiv.org/abs/0811.0464) : amélioration par blocs et flot de renormalisation ; valeur triangulaire de première approximation $`0.835985`$.
