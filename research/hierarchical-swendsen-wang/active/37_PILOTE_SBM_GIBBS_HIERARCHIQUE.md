# Pilote obligatoire : tester la stratégie hiérarchique sur le SBM

**Statut : calibration exacte du broadcast $`\mathrm{PGW}(d)`$ et no-go de
la coupe partagée ; transfert au SBM fini, preuve dynamique et seuils
almost/exact par la hiérarchie encore ouverts.**

**Question testée :** la stratégie

```math
\text{dendrogramme}
\ \longrightarrow\
\text{coupe à }\beta_c
\ \longrightarrow\
\text{Gibbs des sous-arbres}
```

retrouve-t-elle le seuil exact de weak recovery du SBM symétrique à deux
communautés ?

**Réponse courte :** le Gibbs exact sur l'arbre entier fournit une
formulation cohérente et falsifiable, mais il faut l'étudier dans le secteur
à deux répliques.
Conditionner une unique réplique par tout son dendrogramme fige les parités
sur l'arbre local du SBM et ne retrouve pas Kesten--Stigum. En revanche,
deux Gibbs postérieurs exacts, munis de deux hiérarchies indépendantes et
entièrement marginalisées, ont un Jacobien local égal à $\theta^2$. Leur
loi marginale est celle du broadcast classique, dont la densité d'évolution
a le seuil $`d\theta^2=1`$.

Cette dernière identité vaut pour **toute** coupe exactement marginalisée.
Elle calibre le bookkeeping répliqué, mais ne démontre pas que
$`\beta_c`$, un sweep hiérarchique ou la dynamique sur le SBM fini produit
le seuil. Sur le graphe fini, la contrainte de balance ou les non-arêtes
recouplent les arbres finaux.

La coupe à $`\beta_c`$ reste obligatoire **dans la stratégie testée** : elle
décompose chaque arbre en blocs critiques et fixe l'ordre d'élimination.
Elle n'est pas l'origine de Kesten--Stigum. Aucun facteur postcritique n'est
supprimé, remplacé par un canal plus bruité ou traité indépendamment. Tous
sont conservés dans le Gibbs joint de l'arbre final.

Ce résultat valide le **bookkeeping répliqué** sur l'arbre de broadcast. Il
ne prouve pas encore qu'un pas du noyau hiérarchique, ni son transfert au
graphe SBM fini, possède à lui seul le seuil exact.

Cette note est un test go/no-go de la stratégie, pas une nouvelle preuve du
théorème SBM.

## 1. Le benchmark exact

Considérons le SBM équilibré à deux classes :

```math
\mathbb P(A_{ij}=1\mid X_iX_j=+1)
=
\frac an,
\qquad
\mathbb P(A_{ij}=1\mid X_iX_j=-1)
=
\frac bn,
\qquad
a>b>0.
\qquad\text{(1.1)}
```

Posons

```math
d
=
\frac{a+b}{2},
\qquad
\theta
=
\frac{a-b}{a+b},
\qquad
\lambda
=
d\theta^2
=
\frac{(a-b)^2}{2(a+b)}.
\qquad\text{(1.2)}
```

Le seuil exact de weak recovery est

```math
\boxed{\lambda=1.}
\qquad\text{(1.3)}
```

La recovery est impossible pour $`\lambda\le1`$ et possible pour
$`\lambda>1`$. Le voisinage local marqué d'un sommet uniforme converge vers
un $`\mathrm{PGW}(d)`$ sur lequel chaque arête transmet un bit par un canal
binaire symétrique de corrélation $\theta$.

Le pilote travaille d'abord sur cette limite locale de broadcast. Le passage
au graphe fini n'est pas une simple identification arête par arête. Deux
formulations exactes mettent en évidence le verrou.

Avec des labels i.i.d., posons

```math
h_1
=
\frac12\log\frac ab,
\qquad
h_0
=
\frac12
\log\frac{1-a/n}{1-b/n}.
\qquad\text{(1.4)}
```

Pour une matrice d'adjacence $A$, le posterior exact est

```math
\mu_A(x)
\propto
\exp\left[
h_0\sum_{i<j}x_ix_j
+
(h_1-h_0)
\sum_{\{i,j\}:A_{ij}=1}x_ix_j
\right].
\qquad\text{(1.5)}
```

Le premier terme vaut

```math
\frac{h_0}{2}
\left[
\left(\sum_i x_i\right)^2-n
\right]
\qquad\text{(1.6)}
```

et couple toutes les racines du dendrogramme des arêtes présentes.

Dans le planted bisection exactement équilibré, (1.6) est constant, mais la
contrainte $`\sum_i x_i=0`$ couple à son tour les orientations des racines.
Le taux arête-seule fini exact est alors

```math
u_n
=
\log
\frac{
a(1-b/n)
}{
b(1-a/n)
},
\qquad\text{(1.7)}
```

et non exactement $`\log(a/b)`$. Ainsi, on peut éliminer les non-arêtes ou
avoir un a priori produit, mais pas obtenir gratuitement les deux à la fois.
La présente note ne prétend pas avoir contrôlé ce port global ni effectué le
transfert broadcast--graphe.

Le seuil global est le théorème de
[Mossel--Neeman--Sly](https://arxiv.org/abs/1311.4115), obtenu aussi par
[Massoulié](https://arxiv.org/abs/1311.3085). Le pilote ci-dessous demande si
la dynamique proposée retrouve ce nombre sans l'insérer dans ses
définitions.

Ce benchmark sépare trois seuils qu'une stratégie peut confondre :

| objet | paramètre critique |
|---|---:|
| composante géante du graphe | $`d=1`$ |
| composante géante Swendsen--Wang du broadcast-tree | $`d\theta=1`$ |
| weak recovery / Kesten--Stigum | $`d\theta^2=1`$ |

Une méthode qui renvoie $`d=1`$ ou $`d\theta=1`$ ne voit pas le bon secteur
d'information.

## 2. Transposition littérale du dendrogramme

Conditionnellement à l'existence d'une arête du SBM, la relation des labels
est satisfaite avec probabilité

```math
p
=
\frac a{a+b}
=
\frac{1+\theta}{2},
\qquad
u
=
\log\frac p{1-p}
=
\log\frac ab.
\qquad\text{(2.1)}
```

Donnons une horloge $`\mathrm{Exp}(u)`$ à chaque arête satisfaite d'une
réplique, et aucune horloge finie à une arête non satisfaite. La probabilité
annealed qu'une arête présente soit ouverte avant $\beta$ vaut

```math
q_p(\beta)
=
p(1-e^{-u\beta}),
\qquad
q_p(1)=\theta.
\qquad\text{(2.2)}
```

La coupe géométrique critique du graphe local vérifie

```math
dq_p(\beta_c^{\mathrm{geom}})
=
1,
\qquad
\beta_c^{\mathrm{geom}}
=
-\frac1u
\log\left(1-\frac1{dp}\right).
\qquad\text{(2.3)}
```

La solution finie de (2.3) existe si $`dp>1`$ ; elle appartient à l'horizon
$`[0,1]`$ du dendrogramme si et seulement si $`d\theta\ge1`$. Au seuil de
Kesten--Stigum, $`d\theta^2=1`$ implique
$`d\theta=1/\theta>1`$ dès que $`\theta<1`$ : la forêt finale possède donc
déjà une composante géante. Le benchmark reproduit bien la situation du GSBM
triangulaire.

## 3. No-go exact : figer le dendrogramme complet

Le SBM sparse est localement arborescent. Une fusion typique de son
dendrogramme est donc portée par un bucket contenant une seule arête. Pour
une telle fusion $u$,

```math
\Lambda_u(\sigma)
=
u\,
\mathbf1_{\{
\text{l'arête de }u\text{ est satisfaite par }\sigma
\}}.
\qquad\text{(3.1)}
```

Dans la conditionnelle à dendrogramme fixé, le facteur

```math
F_u(\sigma)
=
\Lambda_u(\sigma)
\exp\left[
(1-\beta_u)\Lambda_u(\sigma)
\right]
\qquad\text{(3.2)}
```

est nul si le flip relatif des deux enfants rend l'unique arête non
satisfaite. La fusion impose donc leur parité exactement.

### Proposition 3.1 — échec du dendrogramme figé

Sur une composante arborescente, conditionner par tout le dendrogramme final
fixe toutes les parités relatives le long de ses arêtes de fusion. Le Gibbs
de cette racine ne conserve qu'un flip global uniforme.

Par conséquent, couper cette racine à $`\beta_c^{\mathrm{geom}}`$ puis tirer
un Gibbs **joint** de ses sous-arbres ne crée aucune perte entre deux
sous-arbres reliés par les fusions postcritiques : ces fusions les recollent
par des contraintes dures.

La persistance de paire est alors celle de la percolation Swendsen--Wang
finale, de degré moyen $`d\theta`$. Cette version ne peut prouver
l'impossibilité que dans le régime

```math
d\theta\le1,
\qquad\text{(3.3)}
```

strictement plus petit que le vrai régime $`d\theta^2\le1`$.

Par exemple, $`d=3`$ et $`\theta=1/2`$ donnent

```math
d\theta
=
1.5>1,
\qquad
d\theta^2
=
0.75<1.
\qquad\text{(3.4)}
```

Le dendrogramme figé voit une géante parfaitement corrélée alors que la weak
recovery est impossible. C'est un contre-test exact de la formulation
littérale sur le broadcast libre, pas une objection à l'idée de Gibbs par
sous-arbres. Sur le SBM fini équilibré, le port global de balance doit être
ajouté avant de transposer ce no-go.

## 4. Gibbs exact sur tout l'arbre : formulation sans perte

Pour remplacer la vérité par la sortie d'un noyau, le chapitre 11 impose que
ce noyau laisse la postérieure invariante. Son théorème de percolation exige
en plus un recoloriage indépendant et uniforme des clusters. La construction
proposée sort précisément de cette seconde hypothèse et utilise le Gibbs
exact de l'arbre entier.

Pour une observation $O$, notons $`\mu_O`$ la postérieure et
$`R_O(dD\mid\sigma)`$ la loi conditionnelle du dendrogramme complet, horloges
postcritiques comprises. La mesure augmentée est

```math
\nu_O(d\sigma,dD)
=
\mu_O(d\sigma)R_O(dD\mid\sigma).
\qquad\text{(4.1)}
```

À $D$ fixé, soit $`G_D`$ un noyau qui préserve exactement
$`\nu_O(\cdot\mid D)`$. Dans le broadcast ou le modèle edge-only avec a
priori produit, la version sans perte la plus directe est le heat bath
global : elle tire exactement le Gibbs joint de chaque arbre final, avec des
hasards indépendants entre racines. À l'intérieur d'une même racine, toutes
les fusions et tous les facteurs ancêtres sont conservés.

Cette factorisation ne vaut pas telle quelle sur le SBM fini. Si $z_R$ est
l'orientation d'une racine et $m_R$ sa magnétisation interne, le planted
bisection porte au minimum le facteur global

```math
\mathbf1_{\{
\sum_R z_Rm_R=0
\}},
\qquad\text{(4.2a)}
```

tandis que le modèle i.i.d. porte le potentiel mean-field de (1.6). Le
prochain transfert doit traiter ce facteur comme un port global, ou montrer
qu'il est asymptotiquement négligeable pour l'observable étudiée ; aucune de
ces deux étapes n'est encore établie.

Le noyau marginal

```math
K_O(\sigma,d\sigma')
=
\int
R_O(dD\mid\sigma)
G_D(\sigma,d\sigma')
\qquad\text{(4.2)}
```

laisse $`\mu_O`$ invariante, puis $D$ est oublié. Si $`G_D`$ est réversible
pour chaque conditionnelle, (4.2) est réversible par augmentation de données.
Il n'y a ici ni canal dégradé ni suppression des fusions supérieures.

Cette propriété ne suffit toutefois pas à prouver le mélange ou la
non-reconstruction. Sur l'arbre local du SBM, le Gibbs global à $D$ fixé
conserve toutes les parités internes et ne tire qu'un flip par racine finale.
Le noyau (4.2) y coïncide donc avec le couplage Swendsen--Wang usuel. Son
invariance n'explique pas encore Kesten--Stigum.

### 4.1 La coupe physique exacte

Pour une arête orientée parent--enfant, écrivons

```math
P_\theta(t\mid s)
=
\frac{1+\theta st}{2},
\qquad
s,t\in\{-1,+1\}.
\qquad\text{(4.3)}
```

À un niveau $\beta$, introduisons le bit de coupe $B_e$ par

```math
\mathbb P(B_e=1\mid s,t)
=
\left(1-e^{-u\beta}\right)
\mathbf1_{\{s=t\}}.
\qquad\text{(4.4)}
```

Le facteur joint d'arête est

```math
\Psi_{\theta,\beta}(s,t,b)
=
P_\theta(t\mid s)
\mathbb P(B_e=b\mid s,t),
\qquad
\sum_{b=0}^1
\Psi_{\theta,\beta}(s,t,b)
=
P_\theta(t\mid s).
\qquad\text{(4.5)}
```

À $`\beta=\beta_c^{\mathrm{geom}}`$,

```math
\mathbb P(B_e=1)
=
q_p(\beta_c^{\mathrm{geom}})
=
\frac1d.
\qquad\text{(4.6)}
```

La coupe demandée découpe donc bien le $`\mathrm{PGW}(d)`$ en composantes
critiques. Conditionnellement à $`B_e=1`$, l'arête impose $`s=t`$.
Conditionnellement à $`B_e=0`$, son canal résiduel a pour corrélation

```math
\theta_{\beta}^{\mathrm{res}}
=
\frac{
p e^{-u\beta}-(1-p)
}{
p e^{-u\beta}+(1-p)
}.
\qquad\text{(4.7)}
```

À la coupe critique, l'identité des espérances donne aussi

```math
\theta_c^{\mathrm{res}}
=
\frac{\theta-1/d}{1-1/d}.
\qquad\text{(4.8)}
```

Révéler la coupe augmente strictement le second moment transmis :

```math
\begin{aligned}
\eta_{\mathrm{cut}}
&=
\mathbb E\left[
\mathbb E[st\mid B_e]^2
\right]
\\
&=
\frac1d
+
\left(1-\frac1d\right)
\left(\theta_c^{\mathrm{res}}\right)^2
\\
&=
\theta^2
+
\frac{
(1/d)(1-\theta)^2
}{
1-1/d
}
>
\theta^2.
\end{aligned}
\qquad\text{(4.9)}
```

Cette distinction ne dépend pas du choix critique. Pour
$`q=q_p(\beta)`$, posons

```math
\pi_1=q,
\qquad
\pi_0=1-q,
\qquad
c_1=1,
\qquad
c_0=\frac{\theta-q}{1-q}.
\qquad\text{(4.9a)}
```

La marginalisation d'une coupe donne

```math
\sum_b\pi_bc_b=\theta.
\qquad\text{(4.9b)}
```

Par conséquent, deux coupes indépendantes donnent, pour tout $\beta$,

```math
\sum_{b_1,b_2}
\pi_{b_1}\pi_{b_2}c_{b_1}c_{b_2}
=
\theta^2,
\qquad\text{(4.9c)}
```

alors qu'une coupe partagée donne

```math
\sum_b\pi_bc_b^2
=
\theta^2
+
\frac{q(1-\theta)^2}{1-q}.
\qquad\text{(4.9d)}
```

Ainsi $`\beta_c`$ organise la géométrie, mais n'est pas sélectionné par le
facteur d'information $`\theta^2`$. À criticité exacte, la susceptibilité des
tailles de blocs diverge en outre ; une coupe légèrement sous-critique doit
rester dans les comparaisons numériques de variance.

Par exemple, pour $`d=3`$ et $`\theta=1/2`$,

```math
d\theta^2
=
0.75,
\qquad
d\eta_{\mathrm{cut}}
=
1.125.
\qquad\text{(4.10)}
```

Ainsi, une preuve qui conditionne les deux objets comparés par la même coupe
critique crée artificiellement un régime supercritique sous le vrai seuil.
La coupe doit être utilisée dans le calcul exact, puis marginalisée dans
chaque réplique.

### 4.2 Lift canonique à deux répliques

La quantité de weak recovery est quadratique. Pour deux répliques
postérieures indépendantes conditionnellement à $O$, munies de deux
dendrogrammes indépendants conditionnellement à leurs spins, la loi exacte
est

```math
\widehat\nu_O
=
\nu_O(d\sigma^{(1)},dD^{(1)})
\otimes
\nu_O(d\sigma^{(2)},dD^{(2)}).
\qquad\text{(4.11)}
```

L'indépendance des deux hiérarchies est essentielle. Imposer
$`D^{(1)}=D^{(2)}`$ ou seulement
$`B^{(1)}=B^{(2)}`$ révèle une variable auxiliaire commune et produit
l'inflation (4.9). Il ne suffit pas non plus de tirer deux dendrogrammes
conditionnellement à la même vérité : ils resteraient corrélés par cette
vérité. Il faut deux tirages i.i.d. de la mesure augmentée conditionnellement
à la seule observation.

Pour une perturbation d'overlap

```math
L(t_1,t_2)
=
1+\varepsilon t_1t_2+O(\varepsilon^2),
\qquad\text{(4.12)}
```

l'opérateur marginal d'une arête vérifie

```math
\begin{aligned}
&\sum_{\substack{
t_1,t_2\\
b_1,b_2
}}
\Psi_{\theta,\beta}(s_1,t_1,b_1)
\Psi_{\theta,\beta}(s_2,t_2,b_2)
L(t_1,t_2)
\\
&\qquad
=
1+\varepsilon\theta^2s_1s_2
+O(\varepsilon^2).
\end{aligned}
\qquad\text{(4.13)}
```

Posons $`Y_v=\sigma_v^{(1)}\sigma_v^{(2)}`$. Sur une arête sans champ
extérieur, la dérivée d'un message d'une réplique vaut $\theta$. La dérivée
du message overlap exact vaut donc

```math
J_{\mathrm{overlap}}
=
\theta\theta
=
\theta^2.
\qquad\text{(4.14)}
```

Comme le nombre moyen de descendants est $d$, la linéarisation du Gibbs
répliqué entier est

```math
\mathcal L_{\mathrm{SBM}}
=
d\theta^2.
\qquad\text{(4.15)}
```

Le carré est ici dérivé du produit exact de deux Jacobiennes Gibbs. Il n'est
pas inséré en remplaçant la coupe physique par une percolation
d'information.

Les deux copies de (4.11) sont déjà à l'équilibre. L'identité (4.15) ne
contrôle donc ni un sweep du noyau $K_O$, ni son gap spectral, ni son temps de
mélange. Une preuve réellement dynamique doit majorer, pour un nombre de
pas explicite,

```math
\left\|
K_O^t f_{ij}
\right\|_{L^2(\mu_O)}^2
\qquad\text{(4.15a)}
```

ou un fonctionnel spectral équivalent.

### 4.3 La coupe à $\beta_c$ comme séparateur exact

Pour un arbre fini $T$ et une condition au bord $L_{\partial T}$,

```math
Z_T(s_\rho)
=
\sum_{\sigma_{T\setminus\{\rho\}}}
\left[
\prod_{e=(v,w)\in T}
P_\theta(\sigma_w\mid\sigma_v)
\right]
L_{\partial T}(\sigma).
\qquad\text{(4.16)}
```

En insérant (4.5) à $`\beta_c`$,

```math
Z_T(s_\rho)
=
\sum_B
\sum_{\sigma_{T\setminus\{\rho\}}}
\left[
\prod_{e=(v,w)\in T}
\Psi_{\theta,\beta_c}
(\sigma_v,\sigma_w,B_e)
\right]
L_{\partial T}(\sigma).
\qquad\text{(4.17)}
```

On somme d'abord les spins à l'intérieur des composantes critiques, puis
les séparateurs, toutes les fusions postcritiques et enfin les bits de coupe.
L'associativité des sommes finies garantit que ce Gibbs hiérarchique est
exactement le Gibbs de l'arbre entier.

La même élimination est faite séparément dans $`D^{(1)}`$ et
$`D^{(2)}`$. Aucun sous-arbre critique d'une même racine n'est déclaré
indépendant après marginalisation de ses séparateurs.

Si $m_{v\to u}$ est la magnétisation transmise par le sous-arbre enraciné en
$v$, la récursion exacte d'une réplique est

```math
m_{u\to\mathrm{par}(u)}
=
\tanh\left[
\sum_{v\text{ enfant de }u}
\mathrm{atanh}(\theta m_{v\to u})
\right].
\qquad\text{(4.18)}
```

La coupe critique n'est donc pas une approximation de (4.18). Elle en est
un ordre d'élimination adapté à la géométrie que l'on veut ensuite transférer
au GSBM.

## 5. Ce que la calibration certifie exactement

Soit $`M_t`$ la magnétisation postérieure de la racine d'un
$`\mathrm{PGW}(d)`$ tronqué à profondeur $t$, lorsque les spins de la
frontière sont révélés. Posons

```math
q_t
=
\mathbb E[M_t^2].
\qquad\text{(5.1)}
```

$`q_t`$ est à la fois l'énergie $`L^2`$ d'un message Gibbs, une information
$`\chi^2`$ et un overlap à deux répliques. Sous le contrôle de moments
nécessaire pour uniformiser le reste, la linéarisation de (4.18) en zéro
donne

```math
q_{t+1}
=
d\theta^2q_t
+
O(q_t^2).
\qquad\text{(5.2)}
```

Le carré n'est pas une convention : une branche transmet $\theta$ dans le
secteur spin et $\theta^2$ dans le secteur overlap. Le nombre moyen de
branches vaut $d$. Le sandwich suivant, et non le seul terme
$`O(q_t^2)`$, porte la conclusion rigoureuse.

La linéarisation seule ne suffit pas. Ici, on dispose du sandwich global

```math
\boxed{
\ell_t(\lambda)
\le
q_t
\le
r_t(\lambda),
}
\qquad\text{(5.3)}
```

avec

```math
\ell_t(\lambda)
=
\left(
\sum_{s=0}^{t}\lambda^{-s}
\right)^{-1},
\qquad
r_0=1,
\qquad
r_{t+1}
=
1-e^{-\lambda r_t}.
\qquad\text{(5.4)}
```

La borne inférieure vient du score linéaire
$`Z_t=\sum_{|v|=t}\sigma_v`$ et de Cauchy--Schwarz. La borne supérieure est
la probabilité qu'une percolation auxiliaire de transmission
$`\chi^2`$, de nombre moyen d'enfants $\lambda$, atteigne la profondeur $t$.
Elle certifie ici le benchmark d'équilibre ; elle n'est pas identifiée à la
coupe physique et ne constitue pas une preuve nouvelle issue du seul
dendrogramme. Cette majoration est la spécialisation arborescente des
bornes de
[Evans--Kenyon--Peres--Schulman](https://users.cms.caltech.edu/~schulman/Papers/evansKPS00Ising.pdf)
et de l'[information-percolation
$`\chi^2`$](https://arxiv.org/abs/1806.03227).

Les trois régimes du broadcast classique sont alors exacts.

```math
\begin{array}{c|c|c}
\lambda<1
&
r_t\longrightarrow0\text{ exponentiellement}
&
q_t\longrightarrow0
\\
\lambda=1
&
\dfrac1{t+1}\le q_t\le r_t\sim\dfrac2t
&
q_t\longrightarrow0
\\
\lambda>1
&
\displaystyle
\liminf_{t\to\infty}q_t
\ge
\dfrac{\lambda-1}{\lambda}
&
\text{reconstruction}
\end{array}
\qquad\text{(5.5)}
```

Conclusion exacte de la calibration :

```math
\boxed{
\text{Gibbs postérieur exact sur le broadcast-tree}
\quad+\quad
\text{lift à deux répliques}
\quad\Longrightarrow\quad
d\theta^2=1.
}
\qquad\text{(5.6)}
```

La coupe à $`\beta_c`$ fournit une décomposition exacte de ce Gibbs, et
(4.9) montre pourquoi elle doit être marginalisée séparément dans les deux
répliques. Ce qui n'est pas encore prouvé est que la contraction d'un nombre
fixé de pas de (4.2), analysée uniquement par cette coupe, possède elle-même
la frontière de Kesten--Stigum.

## 6. Quel rôle reste-t-il à $\beta_c$ ?

Le benchmark ne dit pas d'abandonner la coupe physique. Il précise son rôle.

- $`\beta_c^{\mathrm{geom}}`$ choisit des blocs de taille critique et organise
  un sweep multiscalaire.
- Il ne détermine pas seul la transition d'information.
- La transition statistique est gouvernée par l'opérateur à deux répliques
  de la Gibbs marginale exacte.
- Sur le broadcast-tree du SBM, sa linéarisation spatiale est le scalaire
  $`d\theta^2`$, indépendamment de l'ordre exact d'élimination des
  sous-arbres.

Autrement dit,

```math
\boxed{
\text{coupe géométrique pour organiser}
\quad+\quad
\text{transfert }L^2\text{ pour décider}.
}
\qquad\text{(6.1)}
```

Une « percolation d'information » de rétention $\theta^2$ est un certificat
externe du second terme. Elle ne doit pas être confondue avec la forêt
géométrique de rétention $`q_p(\beta)`$, ni présentée comme la partie nouvelle
de la stratégie.

## 7. Deux benchmarks supplémentaires : almost exact et exact recovery

La coupe à $`\beta_c`$ ne dépend pas de l'objectif statistique. Ce qui change
est le fonctionnel transporté par le Gibbs hiérarchique.

Considérons maintenant

```math
\mathbb P(A_{ij}=1\mid X_iX_j=+1)
=
\frac{a_n}{n},
\qquad
\mathbb P(A_{ij}=1\mid X_iX_j=-1)
=
\frac{b_n}{n}.
\qquad\text{(7.1)}
```

Posons encore $`d_n=(a_n+b_n)/2`$,
$`p_n=a_n/(a_n+b_n)`$ et $`u_n=\log(a_n/b_n)`$. La coupe imposée par la
stratégie est

```math
d_n p_n
\left(
1-e^{-u_n\beta_{c,n}}
\right)
=
1.
\qquad\text{(7.2)}
```

Elle produit toujours des composantes critiques de degré moyen un, même
quand $`d_n\to\infty`$. Le même sum-product groupé par ces composantes reste
exact après marginalisation des bits de coupe et conservation de tout le
Gibbs postcritique.

### 7.1 Almost exact recovery

L'énergie quadratique ne distingue pas une erreur constante d'une erreur
qui tend vers zéro. Il faut propager le log-likelihood $L_v$ et son
fonctionnel de Hellinger

```math
H_v
=
\mathbb E
\left[
e^{-L_v/2}
\mid X_v=+1
\right].
\qquad\text{(7.3)}
```

Pour un prior i.i.d., révélons les autres labels. Dans le planted bisection,
révéler tous les autres labels dévoilerait $X_v$ par la contrainte de compte :
il faut cacher au moins une paire opposée, ou un sous-ensemble résiduel dont
la taille diverge. Sous l'approximation équilibrée avec $m_n$ labels révélés
dans chaque groupe, le coefficient de Bhattacharyya exact du test local vaut

```math
H_{v,n}
=
\left[
\frac{\sqrt{a_nb_n}}n
+
\sqrt{
\left(1-\frac{a_n}n\right)
\left(1-\frac{b_n}n\right)
}
\right]^{2m_n}.
\qquad\text{(7.4)}
```

Lorsque $`2m_n=n+O(1)`$ et $`a_n,b_n=O(\log n)`$, son exposant de Chernoff
est

```math
C_n
=
-\log H_{v,n}
=
\frac{
\left(
\sqrt{a_n}-\sqrt{b_n}
\right)^2
}{2}
+
o(1).
\qquad\text{(7.5)}
```

Si $`P_{e,n}^{\mathrm{oracle}}`$ est l'erreur de Bayes locale, les inégalités
élémentaires entre variation totale et affinité donnent

```math
\frac{H_{v,n}^2}{4}
\le
P_{e,n}^{\mathrm{oracle}}
\le
\frac{H_{v,n}}2.
\qquad\text{(7.6)}
```

Le test oracle a donc une erreur tendant vers zéro exactement lorsque
$`C_n\to\infty`$ dans ce régime. Le benchmark almost exact demandé à la
hiérarchie est

```math
C_n\longrightarrow\infty
\quad\Longrightarrow\quad
\mathbb P(L_vX_v\le0)\longrightarrow0.
\qquad\text{(7.7)}
```

À degré borné, $C_n$ reste borné et une proportion positive de sommets est
isolée. La dynamique doit alors conclure que l'almost exact recovery est
impossible ; toute conclusion contraire invaliderait la méthode.

Cette cible est cohérente avec les résultats d'accuracy optimale en degré
divergent d'[Abbe--Sandon](https://arxiv.org/abs/1506.03729). Elle reste ici
un benchmark nécessaire. Pour la suffisance globale, il faut encore
construire une initialisation presque exacte, puis un raffinement
leave-one-out ou par séparation aléatoire des arêtes, et contrôler la
contamination créée par les labels initialement faux. Aucune concentration
du Gibbs hiérarchique à cette échelle n'est encore démontrée.

### 7.2 Exact recovery

Dans le régime logarithmique

```math
a_n=A\log n,
\qquad
b_n=B\log n,
\qquad
A>B>0,
\qquad\text{(7.8)}
```

(7.5) devient

```math
C_n
=
\frac{
\left(
\sqrt A-\sqrt B
\right)^2
}{2}
\log n
+
o(\log n).
\qquad\text{(7.9)}
```

Pour contrôler simultanément les $n$ sommets, l'exposant doit dépasser
$`\log n`$. La frontière de premier ordre à retrouver est donc

```math
\boxed{
\left(
\sqrt A-\sqrt B
\right)^2
=
2.
}
\qquad\text{(7.10)}
```

Au-dessus de cette frontière, une initialisation presque exacte suivie d'un
raffinement local doit fournir une erreur $`o(1/n)`$ uniformément en $v$.
Le regroupement Gibbs à $`\beta_{c,n}`$ ne donne pas encore ce résultat. En
dessous, (7.9) suggère l'impossibilité, mais l'affinité seule ne la prouve
pas : il faut l'asymptotique précise de l'erreur locale

```math
P_{e,n}^{\mathrm{oracle}}
=
n^{-C_n/\log n+o(1)}
\qquad\text{(7.11)}
```

et un second moment du nombre de sommets ambigus. Le traitement du cas
d'égalité exige des corrections de second ordre et n'est pas inclus dans ce
premier pilote.

La frontière stricte est celle du théorème d'exact recovery
d'[Abbe--Bandeira--Hall](https://arxiv.org/abs/1405.3267).

### 7.3 Le lift hiérarchique adapté à Hellinger

Almost exact et exact recovery ne demandent pas simplement deux répliques
ordinaires. Si $`W_\pm(o,D)`$ est le poids non normalisé de l'observation
locale $o$ et de son dendrogramme sous les deux hypothèses sur $X_v$, alors
l'affinité exacte marginale est

```math
H_v
=
\sum_o
\sqrt{
\left[
\int W_+(o,D_+)\,dD_+
\right]
\left[
\int W_-(o,D_-)\,dD_-
\right]
}.
\qquad\text{(7.12)}
```

Chaque dendrogramme doit être sommé **à l'intérieur** de sa propre fonction
de partition avant la racine carrée. La quantité à dendrogramme partagé,

```math
\int
\sum_o
\sqrt{
W_+(o,D)W_-(o,D)
}\,dD,
\qquad\text{(7.13)}
```

correspond à une expérience augmentée plus informative et n'est pas
$`H_v`$. La cible algorithmique est donc un opérateur de
Bhattacharyya--Chernoff sur deux éliminations séparées, pas l'opérateur
quadratique d'overlap réutilisé sans modification.

### 7.4 Même hiérarchie, trois fonctionnels

| objectif | quantité propagée | seuil à retrouver |
|---|---|---|
| weak recovery | $`\mathbb E[M^2]`$ | $`d\theta^2=1`$ |
| almost exact | $`\mathbb E[e^{-L/2}]`$ | $`C_n\to\infty`$ |
| exact | $`n\,\mathbb P(LX\le0)`$ | $`(\sqrt A-\sqrt B)^2=2`$ |

Ainsi, la généralisation ne demande pas trois géométries incompatibles, mais
elle demande des lifts probabilistes différents. Pour une quantité
quadratique, les deux hiérarchies répliquées restent indépendantes
conditionnellement à l'observation. Pour Hellinger, les deux sommes de
(7.12) sont marginalisées séparément avant de prendre la moyenne
géométrique.

## 8. Transfert intelligent vers le GSBM triangulaire

La cible triangulaire doit maintenant reproduire la même architecture, sans
postuler un mélange global.

### 8.1 Objet à construire

Tirer deux répliques postérieures et leurs deux dendrogrammes complets
$`D^{(1)},D^{(2)}`$, indépendamment conditionnellement à l'observation.
Couper chacun à $`\beta_c(p)`$, sans supprimer sa partie postcritique. Pour
chaque paire de types d'interfaces $s$, définir le transfert linéarisé du
Gibbs répliqué exact

```math
\mathcal L_{p}(s,ds').
\qquad\text{(8.1)}
```

Cet opérateur doit combiner :

- le nombre de descendants de type $s'$ créés par une paire de cellules ;
- le produit exact des deux Jacobiennes Gibbs ;
- les deux géométries critiques, leurs rangs réels et leurs buckets
  multiports ;
- tous les facteurs postcritiques de chaque arbre entier ;
- le biais de paire $`|A||B|`$ dans les composantes géantes.

Sur le pilote SBM, il doit se réduire identiquement à

```math
\mathcal L_{\mathrm{SBM}}
=
d\theta^2.
\qquad\text{(8.2)}
```

### 8.2 Critère falsifiable

La première question n'est plus « le Gibbs géant mélange-t-il ? », mais

```math
\rho(\mathcal L_p)<1\ ?
\qquad\text{(8.3)}
```

Si la réponse est oui avec une marge certifiée, on cherche ensuite une
enveloppe non linéaire analogue à $r_{t+1}=1-e^{-\lambda r_t}$. Si la réponse
est non dès la linéarisation, cette famille de blocs ne peut pas améliorer la
borne et doit être abandonnée.

### 8.3 Condition exacte de non-perte

Pour le GSBM aussi, un bucket postcritique de taille un est une contrainte
dure lorsque le dendrogramme complet est fixé. Une chaîne macroscopique de
tels buckets crée un canal parfaitement persistant dans une conditionnelle
quenched. Le calcul répliqué doit donc sommer exactement chaque dendrogramme
dans sa propre copie de la mesure jointe. Il ne doit ni partager la coupe
entre les répliques, ni contracter une fusion postcritique, ni tirer
artificiellement les sous-arbres critiques indépendamment.

À $`D^{(r)}`$ fixé, les racines finales de la réplique $r$ se factorisent et
leurs Gibbs peuvent être tirés avec des hasards indépendants. Après
marginalisation de $`D^{(r)}`$, on obtient un mélange de produits : cette
indépendance ne doit pas être prolongée au-delà de son conditionnement
valide.

## 9. Portes go/no-go révisées

| porte | test | verdict requis |
|---|---|---|
| SBM0 | figer tout le dendrogramme du broadcast libre | doit reproduire le no-go $`d\theta=1`$ |
| SBM1 | Gibbs exact du broadcast, groupé à $`\beta_c`$ | doit redonner la récursion (4.18) sans perte |
| SBM2 | deux hiérarchies indépendantes, jamais une coupe partagée | doit donner exactement $`d\theta^2`$ |
| SBM3 | fermer le régime non linéaire | doit reproduire le sandwich (5.3) |
| SBM-F | ajouter le port global de balance ou des non-arêtes sur le graphe fini | doit comparer quantitativement l'overlap fini au broadcast |
| SBM4 | propager le fonctionnel de Hellinger (7.12) | doit retrouver (7.5) puis les frontières almost/exact |
| TRI0 | définir $`\mathcal L_p`$ sur la paire de quotients critiques géants | aucun mélange global supposé avant cette définition |
| TRI1 | estimer puis certifier $`\rho(\mathcal L_{0.81})`$ | continuer seulement si la marge est strictement négative |
| TRI2 | construire une enveloppe non linéaire | annoncer une borne seulement après fermeture pairwise |

## 10. Diagnostic reproductible

Le module
[`sbm_broadcast_density_evolution.py`](../computations/sbm_broadcast_density_evolution.py)
simule la récursion (4.18) et compare $`\widehat q_t`$ aux deux certificats de
(5.3). Par exemple :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_broadcast_density_evolution.py \
  --degree 3 \
  --lambdas 0.8 0.95 1.0 1.05 1.2 \
  --depth 30 \
  --particles 50000 \
  --batches 8 \
  --seed 20260726
```

Le module
[`sbm_critical_cut_replica_diagnostic.py`](../computations/sbm_critical_cut_replica_diagnostic.py)
calcule en outre exactement l'inflation (4.9). Par exemple :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_critical_cut_replica_diagnostic.py \
  --degree 3 \
  --theta 0.5
```

Le Monte-Carlo illustre la densité d'évolution ; les bornes
$`\ell_t,r_t`$ portent la conclusion analytique sur le broadcast. Le second
diagnostic est une identité d'arête exacte, pas une extrapolation numérique.

## 11. Verdict

La stratégie initiale survit comme calibration du broadcast sous la forme
sans perte suivante :

```math
\boxed{
\text{Gibbs exact de chaque arbre entier}
\ \longrightarrow\
\text{couper à }\beta_c\text{ pour éliminer}
\ \longrightarrow\
\text{conserver et sommer tout le postcritique}
\ \longrightarrow\
\text{répliquer indépendamment}
\ \longrightarrow\
\text{analyser l'overlap}.
}
\qquad\text{(11.1)}
```

La coupe critique reste pertinente pour choisir les échelles. Dans le
broadcast edge-only, l'indépendance reste pertinente entre branches
**conditionnellement à leurs séparateurs**, entre racines finales à
dendrogramme fixé, et entre les deux hiérarchies conditionnellement à
l'observation. Sur le SBM fini, le port global de (1.6) ou (4.2a) interdit
encore la factorisation par racines. Ce qui doit être retiré dans tous les
cas est une coupe commune révélée aux deux répliques : (4.9) montre
exactement qu'elle donne trop d'information.

Le pilote certifie le Jacobien $`d\theta^2`$ et le seuil d'équilibre du
broadcast grâce à des bornes classiques, pas grâce à la coupe. Les preuves
nouvelles encore manquantes sont : contrôler la dynamique pour un nombre de
sweeps explicite, intégrer le port global du SBM fini, puis transférer un
calcul réellement hiérarchique au quotient géant triangulaire. Les seuils
almost/exact restent eux aussi des benchmarks, avec le lift de Hellinger
(7.12).
