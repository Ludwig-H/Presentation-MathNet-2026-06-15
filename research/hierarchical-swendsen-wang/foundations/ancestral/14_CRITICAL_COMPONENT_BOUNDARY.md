# Frontières critiques, composante macroscopique et chaîne des $`\Lambda_v`$

Cette note repart de l'objet exact des slides 31--33. Si
$`u:C=C_1\mathbin{\dot\cup}C_2`$ est un nœud de Kruskal, alors

```math
E_u
=
\bigl\{\{x,y\}\in E:x\in C_1, y\in C_2\bigr\},
\qquad
\Lambda_u(\sigma)
=
\sum_{e\in E_u}|W_e|
\mathbf1_{\{e\text{ satisfait par }\sigma\}}.
```

Il n'existe donc aucun « réservoir global de liens faux » à ajouter au vote.
Une arête interne à un enfant de la fusion n'entre pas dans $`E_u`$ ; son
état de satisfaction est invariant lorsque cet enfant est retourné. Le bon
objet est la loi des marques **sur les frontières successives** du
dendrogramme.

La formalisation sépare trois niveaux.

1. **Algèbre déterministe, établie.** Les arêtes internes aux blocs retournés
   s'annulent exactement dans les quatre poids du heat bath. Une arête qui
   traverse la coupe courante reste en revanche pertinente, même si ses deux
   extrémités appartiendront à une même composante à un temps ultérieur.
2. **Marques de frontière, établies.** Conditionnellement à la partition
   complète au temps $t$, les marques des arêtes entre blocs restent
   indépendantes et possèdent une loi résiduelle explicite. La loi des arêtes
   internes est différente, mais elle n'est pas utilisée dans les
   $`\Lambda_v`$ de la coupe courante.
3. **Géométrie critique, ouverte.** Conditionner deux sommets lointains à
   appartenir à la même composante au seuil sélectionne une composante avec
   un biais de paires. Ce biais modifie les tailles et les formes des coupes,
   le LCA et toute sa chaîne ancestrale. Il ne modifie pas le noyau des
   marques de frontière après désintégration par la partition complète.

Tous les énoncés probabilistes ci-dessous sont sous la loi annealed de
Nishimori, conditionnellement au squelette non marqué lorsque ce
conditionnement est indiqué. Aucune conclusion globale de weak recovery
n'est déduite sans tenir compte de la masse de l'événement de paire. L'ancien
lemme de domination HF est faux en général multiport ; l'alternative active
est le transfert Feynman--Kac aux rangs réels.

## 1. Modèle, horloges et domaine critique

Dans le GSBM homogène binaire, écrivons

```math
O_{xy}=\Sigma_x\Sigma_y Z_{xy},
\qquad
\mathbb P(Z_{xy}=+1)=p,
\qquad
\mathbb P(Z_{xy}=-1)=1-p.
```

Un lien est **conforme** à la ground truth si $`Z_e=+1`$ et **faux** sinon.
Posons

```math
q:=1-p,
\qquad
u_p:=\log\frac pq.
```

Dans la configuration de référence $`\sigma=\Sigma`$, une arête conforme est
satisfaite et reçoit une horloge $`\xi_e\sim\mathrm{Exp}(u_p)`$ ; une arête
fausse est insatisfaite et reçoit $`\xi_e=+\infty`$. Posons

```math
A_e(t)
:=
\mathbf1_{\{Z_e=+1,\ \xi_e\le t\}},
\qquad
\Pi_t
:=
\mathrm{cc}\bigl(V,\{e:A_e(t)=1\}\bigr).
\qquad\text{(1.1)}
```

Les variables $`A_e(t)`$ sont i.i.d. de paramètre

```math
q_p(t)=p(1-e^{-u_pt}).
\qquad\text{(1.2)}
```

Sur la grille triangulaire, notons

```math
q_c=2\sin(\pi/18).
```

Pour tout $`p>1/2`$, l'équation $`q_p(\beta_c)=q_c`$ possède la solution
positive

```math
\boxed{
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
}
\qquad\text{(1.3)}
```

Mais le dendrogramme utilisé par la dynamique est censuré à $t=1$. Ainsi,

```math
\beta_c(p)\le1
\quad\Longleftrightarrow\quad
2p-1\ge q_c
\quad\Longleftrightarrow\quad
p\ge p_{\mathrm{SW}}:=\frac{1+q_c}{2}.
\qquad\text{(1.4)}
```

Pour $`1/2<p<p_{\mathrm{SW}}`$, le temps critique existe sur l'axe non
censuré, mais il est au-delà de la coupe $1$. Toute expérience critique dans
le dendrogramme censuré suppose donc $`p\ge p_{\mathrm{SW}}`$.

## 2. Ce qui entre réellement dans un $`\Lambda_v`$

Pour une partition $`\pi`$ de $V$, notons sa frontière

```math
\partial\pi
:=
\bigl\{\{x,y\}\in E:x\text{ et }y
\text{ appartiennent à deux blocs distincts de }\pi\bigr\}.
\qquad\text{(2.1)}
```

Au nœud $`u:C=C_1\mathbin{\dot\cup}C_2`$, la coupe des slides est
$`E_u=E(C_1,C_2)`$. Elle est une frontière **juste avant** la fusion. Elle ne
contient aucune arête interne à $`C_1`$ ou à $`C_2`$.

### Lemme 2.1 — annulation interne, statut : établi

Pour $`a,b\in\{0,1\}`$, retournons globalement $`C_1`$ lorsque $`a=1`$ et
$`C_2`$ lorsque $`b=1`$. Si les deux extrémités d'une arête $`e=\{x,y\}`$
appartiennent au même bloc retourné, alors

```math
\sigma_x^{ab}\sigma_y^{ab}
=
\sigma_x\sigma_y.
\qquad\text{(2.2)}
```

Son indicateur de satisfaction ne dépend donc pas de $`(a,b)`$. Il en
résulte trois annulations exactes.

1. Une arête interne à $`C_1`$ ou à $`C_2`$ n'entre pas dans $`E_u`$.
2. Les facteurs associés aux descendants stricts de $u$ sont identiques dans
   les quatre états et s'annulent dans la normalisation du heat bath.
3. Dans le facteur de base $`\mu_0(\sigma^{ab})`$, toute contribution interne
   aux blocs retournés est commune aux quatre poids.

La preuve est l'identité (2.2), arête par arête.

> **Distinction temporelle indispensable.** Une arête peut être interne à la
> composante parente $C$ après la fusion tout en traversant
> $`C_1\mid C_2`$ juste avant celle-ci. Elle appartient alors à $`E_u`$ et ne
> doit pas être supprimée. En particulier, une arête fausse traversante ne
> contribue pas à $`\Lambda_u(\Sigma)`$, mais elle contribue au taux d'un
> état de parité impaire. Le critère exact est l'appartenance à la coupe
> courante, pas l'appartenance ultérieure à une même composante.

Pour un ancêtre strict $`v\succ u`$, le même principe s'applique à la coupe
entre ses deux enfants. Le groupe $`E_v^{(0)}`$ qui ne touche aucun des deux
fils de $u$ est invariant sous $`(a,b)`$, mais il reste dans
$`\Lambda_v^{ab}`$. Il ne peut pas être retranché avant d'appliquer
$`F_v(x)=xe^{(1-\beta_v)x}`$, car cette fonction est non linéaire.

## 3. Loi exacte des marques sur une frontière

Fixons $`0\le t\le1`$ et une partition réalisable $`\pi`$. Pour une arête de
$`\partial\pi`$, l'événement $`\Pi_t=\pi`$ impose $`A_e(t)=0`$. Avant
normalisation, les trois catégories résiduelles ont les masses suivantes.

| catégorie sur $`\partial\pi`$ | masse | conforme ? | satisfaite par $`\Sigma`$ ? |
|---|---:|:---:|:---:|
| vraie tardive, $`t<\xi_e\le1`$ | $`pe^{-u_pt}-q`$ | oui | oui |
| vraie censurée, $`\xi_e>1`$ | $`q`$ | oui | oui |
| fausse, $`Z_e=-1`$ | $`q`$ | non | non |

La catégorie vraie précoce $`\xi_e\le t`$ est impossible sur la frontière,
car elle aurait relié les deux blocs.

### Lemme 3.1 — factorisation de frontière, statut : établi

Conditionnellement à $`\Pi_t=\pi`$, les catégories des arêtes de
$`\partial\pi`$ sont indépendantes. Leur loi commune est

```math
\boxed{
(\text{vraie tardive},\text{vraie censurée},\text{fausse})
\sim
\left(
h_p(t),
\frac{1-h_p(t)}2,
\frac{1-h_p(t)}2
\right),
}
\qquad\text{(3.1)}
```

où

```math
h_p(t)
:=
\frac{pe^{-u_pt}-q}{q+pe^{-u_pt}}
=
\tanh\left(\frac{u_p(1-t)}2\right),
\qquad\text{(3.2)}
```

et la probabilité d'être conforme vaut

```math
s_p(t)
:=
\frac{pe^{-u_pt}}{q+pe^{-u_pt}}
=
\frac{1+h_p(t)}2.
\qquad\text{(3.3)}
```

### Preuve

L'événement $`\Pi_t=\pi`$ est l'intersection de deux événements portant sur
des familles disjointes d'arêtes : chaque sous-graphe précoce induit par un
bloc de $`\pi`$ est connexe, et toutes les arêtes de $`\partial\pi`$ sont non
précoces. L'indépendance des marques par arête donne donc la factorisation
sur $`\partial\pi`$. Enfin $`e^{-u_p}=q/p`$ ; les trois masses du tableau,
normalisées par $`q+pe^{-u_pt}`$, donnent (3.1)--(3.3).

La loi des arêtes **internes** ne factorise généralement pas sous ce
conditionnement : elle est biaisée par la contrainte de connexité dans chaque
bloc. C'est précisément pourquoi il serait faux de leur appliquer (3.1).
Cette difficulté n'affecte pas le taux de la frontière courante, car ces
arêtes n'y figurent pas.

Le [fichier 15](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md) donne la désintégration
correcte de ces arêtes internes. Conditionnellement au graphe ouvert complet,
les seules arêtes internes encore fermées retrouvent la loi résiduelle ; la
connexité et le conditionnement de plus grande composante agissent uniquement
sur leur nombre. Il en déduit aussi la probabilité exacte des deux états de
parité paire au LCA critique.

### Spécialisation critique

Supposons $`p\ge p_{\mathrm{SW}}`$. L'identité
$`q_p(\beta_c)=q_c`$ donne

```math
pe^{-u_p\beta_c}=p-q_c.
\qquad\text{(3.4)}
```

Les masses non conditionnelles d'une arête sont

```math
\begin{array}{c|c}
\text{catégorie}&\text{masse}\cr
\hline
\text{vraie précoce }(\xi_e\le\beta_c)&q_c\cr
\text{vraie tardive }(\beta_c<\xi_e\le1)&2p-1-q_c\cr
\text{vraie censurée }(\xi_e>1)&1-p\cr
\text{fausse}&1-p.
\end{array}
\qquad\text{(3.5)}
```

Pour une arête de $`\partial\Pi_{\beta_c}`$, et seulement après
conditionnement par la partition complète,

```math
s_c(p)
=
\frac{p-q_c}{1-q_c},
\qquad
h_c(p)
=
\frac{2p-1-q_c}{1-q_c}.
\qquad\text{(3.6)}
```

La dérivée

```math
\frac{d}{dt}h_p(t)
=
-\frac{u_p}{2}\left[1-h_p(t)^2\right]
<0
\qquad\text{(3.7)}
```

montre que $`\beta_c`$ maximise la qualité **par arête de frontière et à
coupe fixée** parmi les temps $`t\ge\beta_c`$. Elle n'ordonne ni les tailles
des coupes, ni les partitions sélectionnées.

Le diagnostic conservateur qui compare, sur une frontière critique, les
seules vraies tardives aux fausses devient favorable à

```math
p_{\partial,\mathrm{late}}
:=
\frac{2+q_c}{3}
=
0.7824321184\ldots.
\qquad\text{(3.8)}
```

Il ne s'agit ni d'une proportion globale d'arêtes, ni d'un seuil de weak
recovery. Le paquet complet de frontière a une majorité conforme dès
$`p>p_{\mathrm{SW}}`$.

## 4. La composante critique vue depuis une paire lointaine

Soit $`G_L=(V_L,E_L)`$ une exhaustion triangulaire de diamètre d'ordre $L$,
$`n_L=|V_L|`$, et $`d_L`$ sa distance de graphe. Pour une composante
$`C\in\Pi_{\beta_c}`$ et $`\rho>0`$, posons

```math
N_{L,\rho}(C)
:=
\#\bigl\{(x,y)\in C^2:x\ne y,\
d_L(x,y)\ge\rho L\bigr\}.
\qquad\text{(4.1)}
```

La loi naturelle de la composante contenant une paire uniforme lointaine,
conditionnellement à sa connexion au seuil, est la loi de Palm

```math
\mathbb E_{L,\rho}^{\mathrm{pair}}[F]
=
\frac{
\mathbb E\left[
\displaystyle\sum_{C\in\Pi_{\beta_c}}
\sum_{\substack{x,y\in C\cr x\ne y,\ d_L(x,y)\ge\rho L}}
F(\Pi_{\beta_c},C,x,y)
\right]
}{
\mathbb E\left[
\displaystyle\sum_{C\in\Pi_{\beta_c}}N_{L,\rho}(C)
\right]
}.
\qquad\text{(4.2)}
```

### Proposition 4.1 — biais exact de la paire, statut : établi

Conditionnellement à une partition $`\pi`$, la composante sélectionnée dans
(4.2) a une probabilité proportionnelle à $`N_{L,\rho}(C)`$. La loi de la
partition elle-même est inclinée par
$`\sum_{C\in\pi}N_{L,\rho}(C)`$. Une composante uniforme, la plus grande
composante et une composante pondérée par sa taille ont donc trois lois
différentes de la loi pertinente.

La masse non conditionnelle de l'événement de paire est

```math
\rho_{L,\rho}^{c}
=
\frac1{n_L^2}
\mathbb E\sum_{C\in\Pi_{\beta_c}}N_{L,\rho}(C).
\qquad\text{(4.3)}
```

Cette identité est un simple double comptage des paires ordonnées.

### « Géante » ou « macroscopique » ?

Au paramètre critique de la percolation plane, il n'existe pas de composante
infinie de densité positive. Sous les hypothèses RSW usuelles pour une
exhaustion plane périodique,

```math
\frac1{n_L}
\max_{C\in\Pi_{\beta_c}}|C|
\longrightarrow0
\quad\text{en probabilité}.
\qquad\text{(4.4)}
```

Il est donc trompeur de parler d'une « composante géante au moment exact de
la percolation » si géante signifie $`|C|\ge\alpha n_L`$ avec
$`\alpha>0`$ fixé. L'événement naturel ici est une **composante critique
macroscopique**, par exemple de diamètre au moins $`\rho L`$. La présence de
la paire lointaine dans (4.2) impose déjà cette propriété.

Conditionner malgré tout par $`|C|\ge\alpha n_L`$ au seuil définit une
expérience de grande déviation, potentiellement beaucoup plus informative,
mais ce n'est plus la loi critique typique. Elle ne peut servir de borne sans
un théorème de comparaison séparé.

### Influence exacte du conditionnement de composante

Après désintégration par $`\Pi_{\beta_c}=\pi`$, l'événement de paire ne
dépend plus des marques cachées des arêtes de $`\partial\pi`$. Le lemme 3.1
reste donc valable, avec le même paramètre $`h_c(p)`$. Le conditionnement
change en revanche :

- la loi de $`\pi`$ et de la composante sélectionnée ;
- la taille et la forme de ses frontières ;
- le LCA des deux points et le facteur de biais entre ses deux fils ;
- les tailles des trois groupes et les temps de toute la chaîne ancestrale.

Si l'on conditionne seulement par l'événement grossier « même composante
macroscopique » sans révéler $`\pi`$, les marques de frontière sont une
mixture et ne sont plus indépendantes en général. Toute utilisation de la
loi binomiale doit donc annoncer la désintégration complète.

## 5. LCA critique et loi exacte de sa coupe

Pour deux sommets, posons

```math
\beta_{ij}
:=
\inf\{t:i\leftrightarrow j\text{ dans }\Pi_t\}.
\qquad\text{(5.1)}
```

L'implication déterministe fondamentale est

```math
i\leftrightarrow j\text{ dans }\Pi_{\beta_c}
\quad\Longrightarrow\quad
\beta_{ij}\le\beta_c.
\qquad\text{(5.2)}
```

Le cas favorable conforme à l'idée étudiée doit donc être formulé par la
fenêtre gauche

```math
\mathcal F_{L,\rho,\varepsilon}^{c}
:=
\left\{
d_L(I_L,J_L)\ge\rho L,\
\beta_c-\varepsilon\le\beta_{I_LJ_L}\le\beta_c
\right\}.
\qquad\text{(5.3)}
```

Il impose simultanément une composante critique macroscopique et une
séparation des deux branches au voisinage du seuil. La borne inférieure
$`\beta_c-\varepsilon`$ est asymptotiquement redondante : par la
proposition 5.1 ci-dessous, le conditionnement gauche
$`\{d_L\ge\rho L,\ \beta_{ij}\le\beta_c\}`$ suffit et définit la même
limite conditionnelle ; voir la
[note 42 §4.3](42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md) pour la version
sans fenêtre et le tableau des conventions.

### Proposition 5.1 — localisation critique grossière, statut : établi sous RSW

Pour tout $`\delta>0`$ fixé, sur le tore triangulaire ou dans une exhaustion
avec conditions de bord contrôlées,

```math
\mathbb P\left(
\beta_{I_LJ_L}\le\beta_c-\delta
\ \middle|\
d_L(I_L,J_L)\ge\rho L,\
\beta_{I_LJ_L}\le\beta_c
\right)
\longrightarrow0.
\qquad\text{(5.4)}
```

En effet, $`q_p(\beta_c-\delta)<q_c`$. La probabilité de relier deux points à
distance d'ordre $L$ décroît alors exponentiellement. Au seuil, les
inégalités de box-crossing donnent une minoration polynomiale suffisante du
dénominateur. Le quotient tend donc vers zéro.

Cette proposition localise le LCA à une distance $`o(1)`$ du seuil. Elle ne
donne pas la largeur optimale de la fenêtre proche-critique. L'échelle
$`L^{-3/4}`$ est rigoureusement reliée aux mesures pivotales pour la
percolation **par sites** triangulaire ; son transfert tel quel au modèle de
liens utilisé ici ne doit pas être supposé sans théorème d'universalité assez
fort.

### Contre-audit du mot « favorable »

Comme $`h_p(t)`$ est décroissant, une fusion strictement avant
$`\beta_c`$ possède une meilleure qualité par arête qu'une fusion critique.
Le seuil est le meilleur temps parmi les fusions **postcritiques**, mais le
moins bon parmi les fusions déjà présentes dans $`\Pi_{\beta_c}`$. L'idée
reste cohérente pour une paire lointaine grâce à (5.4) : son LCA ne peut pas
rester uniformément sous-critique. Le gain espéré au seuil vient alors de la
géométrie macroscopique des coupes, pas d'une monotonie scalaire globale.

### Théorème 5.2 — loi conditionnelle d'une coupe de Kruskal, statut : établi

Fixons le squelette non marqué, $`E_u=E(C_1,C_2)`$ de taille $`m\ge1`$ et le
temps $`\beta_u=t`$. L'arête gagnante est conforme. Pour les $`m-1`$ autres
arêtes de cette coupe, notons $R,S,U$ les nombres respectifs d'arêtes vraies
tardives, vraies censurées et fausses. Alors

```math
(R,S,U)
\sim
\mathrm{Mult}
\left(
m-1;
h_p(t),
\frac{1-h_p(t)}2,
\frac{1-h_p(t)}2
\right).
\qquad\text{(5.5)}
```

Le nombre d'arêtes de la coupe satisfaites par la ground truth vaut

```math
K_u
=
1+R+S
\stackrel d=
1+\mathrm{Bin}(m-1,s_p(t)).
\qquad\text{(5.6)}
```

Pour $`V_u:=2K_u-m`$,

```math
\mathbb E[V_u\mid m,t]
=
1+(m-1)h_p(t),
\qquad
\mathrm{Var}(V_u\mid m,t)
=
(m-1)[1-h_p(t)^2],
\qquad\text{(5.7)}
```

et

```math
\mathbb P(V_u>0\mid m,t)
=
\sum_{r=\lfloor m/2\rfloor}^{m-1}
\binom{m-1}{r}
s_p(t)^r[1-s_p(t)]^{m-1-r}.
\qquad\text{(5.8)}
```

Le paramètre de concentration est $`m h_p(t)^2`$. La loi (5.5) concerne la
coupe $`E_u`$, non toutes les arêtes dont les extrémités se trouvent dans la
composante critique finale. Le biais de paire modifie la loi de
$`(m,E_u,t)`$ ; conditionnellement à ces variables, il ne modifie pas le
noyau (5.5).

## 6. Tous les ancêtres : majorités groupées et quatre taux

Soit $`v\succ u`$. Avec les notations du fichier 08, décomposons

```math
E_v=E_v^{(0)}\mathbin{\dot\cup}E_v^{(1)}\mathbin{\dot\cup}E_v^{(2)},
```

où les groupes $1$ et $2$ touchent respectivement les deux fils $`C_1,C_2`$
de $u$, tandis que le groupe $0$ ne touche aucun d'eux. Ici encore, $`E_v`$
est exactement la coupe entre les deux enfants de $v$ : aucune arête interne
à l'un de ces enfants n'est comptée. Le groupe $0$ est une partie de cette
coupe ancestrale, pas un ensemble d'arêtes internes. Dans le cas homogène,
notons

```math
m_{v,r}:=|E_v^{(r)}|,
\qquad
K_{v,r}:=\#\{e\in E_v^{(r)}:Z_e=+1\},
\qquad
M_{v,r}:=2K_{v,r}-m_{v,r}.
```

Conditionnellement au squelette et à la catégorie gagnante
$`G_v\in\{0,1,2\}`$,

```math
\mathbb P(G_v=r\mid\mathscr D)
=
\frac{m_{v,r}}{m_{v,0}+m_{v,1}+m_{v,2}},
\qquad\text{(6.1)}
```

et, conditionnellement à $`G_v`$,

```math
\boxed{
K_{v,r}
=
\mathbf1_{\{G_v=r\}}
+
\mathrm{Bin}
\left(
m_{v,r}-\mathbf1_{\{G_v=r\}},
s_p(\beta_v)
\right).
}
\qquad\text{(6.2)}
```

Les quatre taux sous les flips $`(a,b)\in\{0,1\}^2`$ sont

```math
\boxed{
\frac{\Lambda_v^{ab}}{u_p}
=
K_{v,0}
+
\begin{cases}
K_{v,1},&a=0,\cr
m_{v,1}-K_{v,1},&a=1,
\end{cases}
+
\begin{cases}
K_{v,2},&b=0,\cr
m_{v,2}-K_{v,2},&b=1.
\end{cases}
}
\qquad\text{(6.3)}
```

Par conséquent,

```math
\Lambda_v^{00}-\Lambda_v^{10}=u_pM_{v,1},
\qquad
\Lambda_v^{00}-\Lambda_v^{01}=u_pM_{v,2},
\qquad\text{(6.4)}
```

```math
\Lambda_v^{00}-\Lambda_v^{11}
=
u_p(M_{v,1}+M_{v,2}).
\qquad\text{(6.5)}
```

Une majorité sur tout $`E_v`$ ne suffit donc pas. Pour que le taux de la
ground truth $`\Lambda_v^{00}`$ domine séparément les trois concurrents, il
suffit d'avoir des majorités strictes dans **les deux groupes affectés** :

```math
M_{v,1}>0,
\qquad
M_{v,2}>0.
\qquad\text{(6.6)}
```

Pour chaque groupe, l'échelle de concentration pertinente est

```math
m_{v,r}h_p(\beta_v)^2.
\qquad\text{(6.7)}
```

On a toujours $`\beta_v>\beta_u`$, donc

```math
h_p(\beta_v)<h_p(\beta_u).
\qquad\text{(6.8)}
```

Si $`\beta_v\ge\beta_c`$, alors $`h_p(\beta_v)\le h_c(p)`$. Mais un ancêtre
strict du LCA peut encore appartenir à la fenêtre gauche
$`(\beta_u,\beta_c]`$ ; sa qualité est alors comprise entre
$`h_c(p)`$ et $`h_p(\beta_u)`$. La composante critique et le biais de paire
déterminent combien d'ancêtres tombent de chaque côté du seuil. Attribuer
uniformément $`h_c(p)`$ à toute la chaîne n'est donc justifié dans aucun des
deux sens sans ce contrôle géométrique.

## 7. De la majorité des taux au heat bath de parité

Les slides 31--33 donnent, pour $`a,b\in\{0,1\}`$,

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
F_v(\Lambda_v^{ab}),
\qquad
F_v(x):=xe^{(1-\beta_v)x}.
\qquad\text{(7.1)}
```

Sous l'a priori uniforme, la probabilité que le heat bath conserve la parité
conforme entre $`i`$ et $`j`$ est

```math
P_u^{\mathrm{keep}}
=
\frac{q_u^{00}+q_u^{11}}
{q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}}.
\qquad\text{(7.2)}
```

Posons

```math
L_u
:=
\log\frac{q_u^{00}+q_u^{11}}
{q_u^{10}+q_u^{01}}.
\qquad\text{(7.3)}
```

Alors

```math
P_u^{\mathrm{keep}}>\frac12
\quad\Longleftrightarrow\quad
L_u>0,
\qquad\text{(7.4)}
```

et le biais signé du heat bath vaut $`\tanh(L_u/2)`$.

Lorsque les quatre poids sont strictement positifs, définissons

```math
D_{ab}:=\log\frac{q_u^{00}}{q_u^{ab}}.
```

Le critère exact devient

```math
\boxed{
1+e^{-D_{11}}
>
e^{-D_{10}}+e^{-D_{01}}.
}
\qquad\text{(7.5)}
```

Une condition suffisante simple est

```math
\min(D_{10},D_{01})>\log2.
\qquad\text{(7.6)}
```

Les cas où un taux est nul doivent être traités directement dans (7.2), sans
prendre de logarithme.

### Théorème 7.1 — certificat de majorité hiérarchique, statut : établi

Supposons l'a priori i.i.d. uniforme. Dans le cas pondéré, posons

```math
X_{v,r}:=2\lambda_{v,r}-T_{v,r}.
```

Si

```math
2\Lambda_u-T_u\ge0
\qquad\text{(7.7)}
```

et si, pour tout ancêtre strict $`v\succ u`$,

```math
X_{v,1}\ge0,
\qquad
X_{v,2}\ge0,
\qquad\text{(7.8)}
```

alors

```math
\boxed{
q_u^{00}+q_u^{11}
\ge
q_u^{10}+q_u^{01}.
}
\qquad\text{(7.9)}
```

Si la majorité locale (7.7) est stricte, alors (7.9) est stricte. Dans le cas
homogène, (7.8) est exactement $`M_{v,1},M_{v,2}\ge0`$.

### Preuve

Pour $`v\succ u`$, posons

```math
f_v(a,b):=F_v(\Lambda_v^{ab}),
\qquad
F_v(x)=xe^{(1-\beta_v)x}.
```

Comme $`0\le\beta_v\le1`$, la fonction $`F_v`$ est croissante et convexe sur
$`[0,+\infty)`$ :

```math
F_v'(x)
=
e^{(1-\beta_v)x}\left[1+(1-\beta_v)x\right]>0,
```

```math
F_v''(x)
=
(1-\beta_v)e^{(1-\beta_v)x}
\left[2+(1-\beta_v)x\right]
\ge0.
```

Écrivons les quatre coefficients de Walsh de $`f_v`$ :

```math
\widehat f_v(\varnothing)
=
\frac14(f_v^{00}+f_v^{10}+f_v^{01}+f_v^{11}),
```

```math
\widehat f_v(1)
=
\frac14(f_v^{00}+f_v^{01}-f_v^{10}-f_v^{11}),
```

```math
\widehat f_v(2)
=
\frac14(f_v^{00}+f_v^{10}-f_v^{01}-f_v^{11}),
```

```math
\widehat f_v(12)
=
\frac14(f_v^{00}+f_v^{11}-f_v^{10}-f_v^{01}).
\qquad\text{(7.10)}
```

Le coefficient constant est positif. Sous (7.8), les deux coefficients
linéaires sont non négatifs par croissance. Pour le dernier, les quatre
arguments de $`F_v`$ forment le
rectangle additif

```math
L,\quad L-X_{v,1},\quad L-X_{v,2},
\quad L-X_{v,1}-X_{v,2}.
```

La convexité donne

```math
F_v(L)+F_v(L-X_{v,1}-X_{v,2})
\ge
F_v(L-X_{v,1})+F_v(L-X_{v,2}),
```

donc $`\widehat f_v(12)\ge0`$.

Au nœud $u$, le facteur local prend une valeur sur les états pairs et une
autre sur les états impairs. Sous (7.7), ses seuls coefficients de Walsh non
nuls, $`\widehat f_u(\varnothing)`$ et $`\widehat f_u(12)`$, sont non
négatifs.

Enfin, la transformée de Walsh d'un produit pointwise est la convolution des
transformées. Le cône « quatre coefficients de Walsh non négatifs » est donc
stable par produit. Comme l'a priori uniforme est constant,

```math
\frac14
\left(
q_u^{00}+q_u^{11}-q_u^{10}-q_u^{01}
\right)
\ge0.
```

Si (7.7) est stricte, le coefficient local $`\widehat f_u(12)`$ est strictement
positif et le coefficient constant de chaque autre facteur est positif ; la
convolution donne la stricte positivité finale.

> **Portée.** Le théorème est un certificat suffisant, pas une équivalence.
> Le heat bath peut préférer la parité conforme même si un groupe ancestral
> a une majorité négative, grâce aux autres facteurs. L'a priori uniforme et
> la structure factorisée sont essentiels : pour quatre poids positifs
> arbitraires, la domination séparée de $`q^{00}`$ sur $`q^{10},q^{01}`$ ne
> suffit pas.

### Corollaire 7.2 — probabilité conditionnelle du certificat

Conditionnons par le squelette et par toutes les catégories gagnantes.
Pour $`n\ge0`$, $`g\in\{0,1\}`$ et $`s\in[0,1]`$, définissons les deux
queues binomiales exactes

```math
\mathcal A_{\ge}(n,g,s)
:=
\sum_{b=\lceil(n-g)/2\rceil}^{n}
\binom nb s^b(1-s)^{n-b},
\qquad\text{(7.11)}
```

et

```math
\mathcal A_{>}(n,g,s)
:=
\sum_{b=\lfloor(n-g)/2\rfloor+1}^{n}
\binom nb s^b(1-s)^{n-b}.
\qquad\text{(7.12)}
```

Pour le nœud local, posons $`n_u=m_u-1`$. Pour un groupe ancestral,

```math
g_{v,r}:=\mathbf1_{\{G_v=r\}},
\qquad
n_{v,r}:=m_{v,r}-g_{v,r}.
```

L'indépendance conditionnelle des comptes groupés donne exactement la
probabilité du certificat du théorème 7.1 :

```math
\mathcal C_u(\mathscr D,G)
=
\mathcal A_{>}\left(
n_u,1,s_p(\beta_u)
\right)
\prod_{v\succ u}\prod_{r=1}^2
\mathcal A_{\ge}\left(
n_{v,r},g_{v,r},s_p(\beta_v)
\right).
\qquad\text{(7.13)}
```

Par conséquent,

```math
\boxed{
\mathbb P\left(
P_u^{\mathrm{keep}}>\frac12
\;\middle|\;
\mathscr D,(G_v)_{v\succ u}
\right)
\ge
\mathcal C_u(\mathscr D,G).
}
\qquad\text{(7.14)}
```

Cette borne est exacte pour l'événement suffisant, mais pas nécessairement
pour l'événement de préférence lui-même.

Pour obtenir un certificat plus lisible, posons aussi

```math
\varepsilon(n,g,h)
:=
\begin{cases}
0,&n=0,\cr
\exp\left[-\dfrac{(nh+g)^2}{2n}\right],&n\ge1.
\end{cases}
\qquad\text{(7.15)}
```

Alors

```math
\mathbb P\left(
P_u^{\mathrm{keep}}>\frac12
\;\middle|\;
\mathscr D,(G_v)_{v\succ u}
\right)
\ge
\left[
1
-\varepsilon(n_u,1,h_p(\beta_u))
-\sum_{v\succ u}\sum_{r=1}^2
\varepsilon(n_{v,r},g_{v,r},h_p(\beta_v))
\right]_+.
\qquad\text{(7.16)}
```

Ici $`[x]_+:=\max(x,0)`$.

En effet, un groupe de taille $`m=n+g`$ a une marge

```math
M=g+\sum_{\ell=1}^n X_\ell,
\qquad
\mathbb E[X_\ell]=h,
\qquad
X_\ell\in\{-1,+1\}.
```

Pour $`n\ge1`$, l'inégalité de Hoeffding borne l'événement $`M\le0`$ par
(7.15) ; pour $`n=0`$, chacun des certificats pertinents est déterministe et
valide. Une union bound sur la majorité locale stricte et les deux majorités
ancestrales non négatives, suivie du théorème 7.1, donne (7.16).

Ce corollaire rend le verrou quantitatif explicite. Pour une chaîne de $H$
ancêtres, un régime suffisant est

```math
\min\left\{
(m_u-1)h_p(\beta_u)^2,\,
\min_{\substack{v\succ u\\r\in\{1,2\}}}
m_{v,r}h_p(\beta_v)^2
\right\}
\gg
\log H.
\qquad\text{(7.17)}
```

Il n'est pas nécessaire : le critère exact (7.5) peut rester favorable
lorsque le certificat de majorité échoue.

### Conséquence exacte au seul nœud $u$

Au nœud de fusion,

```math
\Lambda_u^{00}=\Lambda_u^{11}=u_pK_u,
\qquad
\Lambda_u^{10}=\Lambda_u^{01}=u_p(m-K_u).
```

Comme $`F_u`$ est strictement croissante, si l'a priori et les ancêtres sont
neutralisés, le heat bath préfère la parité conforme si et seulement si
$`K_u>m/2`$. Avec les ancêtres présents, cette équivalence locale disparaît :
seul (7.5) décide.

## 8. Ce que cela donne — et ne donne pas — pour la weak recovery

Les trois expériences suivantes ne doivent pas être confondues.

| expérience | contrainte sur le LCA | masse pour une paire uniforme | usage légitime |
|---|---|---:|---|
| connexion critique de Palm | $`\beta_{ij}\le\beta_c`$ | $`\rho_{L,\rho}^{c}`$ | loi naturelle de la composante critique vue par une paire |
| fenêtre critique favorable | $`\beta_c-\varepsilon\le\beta_{ij}\le\beta_c`$ | masse de $`\mathcal F_{L,\rho,\varepsilon}^{c}`$ | oracle géométrique à comparer aux autres paires |
| composante géante à $`t>\beta_c`$ | seulement $`\beta_{ij}\le t`$ | peut être d'ordre un | expérience supercritique, sans localisation automatique du LCA |

Par (4.4),

```math
\rho_{L,\rho}^{c}
\le
\mathbb E\left[
\frac1{n_L}\max_{C\in\Pi_{\beta_c}}|C|
\right]
\longrightarrow0.
\qquad\text{(8.1)}
```

Une paire uniforme appartient donc à une même composante critique
macroscopique avec une probabilité qui s'annule. La loi conditionnelle (4.2)
peut être très informative, mais sa contribution brute au score est

```math
\mathbb E\left[
\eta_{I_LJ_L}^{\mathrm{LCA}}
\mathbf1_{\mathcal F_{L,\rho,\varepsilon}^{c}}
\right]
=
\mathbb P(\mathcal F_{L,\rho,\varepsilon}^{c})
\mathbb E\left[
\eta_{I_LJ_L}^{\mathrm{LCA}}
\mid
\mathcal F_{L,\rho,\varepsilon}^{c}
\right].
\qquad\text{(8.2)}
```

Supprimer le premier facteur est un biais de sélection. La loi conditionnelle
reste un benchmark, mais elle ne domine pas en général les paires
postcritiques multiports. La petitesse de sa masse typique ne suffit donc ni
à la promouvoir en oracle, ni à la réfuter comme cas-test local.

Les seuils scalaires gardent seulement la portée suivante :

```math
\begin{array}{c|c|c}
\text{seuil}&\text{valeur}&\text{interprétation exacte}\cr
\hline
p_{\mathrm{SW}}&(1+q_c)/2=0.673648\ldots
&\beta_c\le1\text{ et majorité d'une frontière critique}\cr
p_{\partial,\mathrm{late}}&(2+q_c)/3=0.782432\ldots
&\text{vraies tardives de frontière majoritaires face aux fausses}\cr
p_{\mathrm{info}}&(1+\sqrt{q_c})/2=0.794659\ldots
&\text{baseline d'impossibilité information--percolation}\cr
p_{\mathrm N}^{(0)}&0.835805\ldots
&\text{calibration entropique exacte d'une face}
\end{array}
\qquad\text{(8.3)}
```

Le test de majorité tardive ne peut pas, seul, améliorer la baseline puisque
$`p_{\partial,\mathrm{late}}<p_{\mathrm{info}}`$. Une majorité locale ne
prouve pas davantage la weak recovery : le score LCA est un majorant de la
corrélation réelle, et le succès d'un oracle conditionnel ne construit aucun
estimateur non oracle.

Une meilleure borne par la dynamique hiérarchique reste possible seulement
si le calcul utilise simultanément :

1. la loi de Palm des tailles et des formes de coupes ;
2. les quatre taux de chaque ancêtre, y compris le baseline non linéaire du
   groupe $`E_v^{(0)}`$ ;
3. la dégradation exacte $`h_p(\beta_v)`$ le long de la chaîne ;
4. le critère pair contre impair (7.5), et non le seul signe d'une majorité ;
5. la masse des paires ou un véritable ordre de domination HF.

Aucun nouveau seuil numérique rigoureux de weak recovery ne découle encore
de ces résultats. Le noyau des marques est désormais fermé ; le verrou est
géométrique et hiérarchique.

## 9. Audit puis contre-audit

| affirmation proposée | audit exact | contre-audit |
|---|---|---|
| Les liens faux internes à un cluster ne doivent pas être comptés. | Vrai s'ils sont internes à un **enfant de la coupe courante** : lemme 2.1. | Faux si « interne » désigne la composante parente à un temps ultérieur : une arête traversant $`C_1\mid C_2`$ appartient à $`E_u`$. |
| Après $`\beta_c`$, il reste tous les faux et les vraies horloges tardives. | Vrai seulement arête par arête sur $`\partial\Pi_{\beta_c}`$, après conditionnement complet. | Les faux internes ont une autre loi et sont exclus ; les vraies censurées font aussi partie du paquet de frontière. |
| Le conditionnement de composante change la qualité de chaque lien. | Faux après désintégration par la partition : le noyau (3.1) est inchangé. | Vrai au niveau grossier : la mixture induit des dépendances et change toute la géométrie des coupes. |
| Deux points dans la composante critique fusionnent après le seuil. | Faux : (5.2) donne $`\beta_{ij}\le\beta_c`$. | Pour une paire lointaine, (5.4) localise néanmoins la fusion vers le seuil par la gauche. |
| Il existe une géante typique exactement au seuil. | Faux si « géante » signifie densité positive : (4.4). | Une composante de diamètre macroscopique existe avec probabilité non négligeable sous des événements de crossing, et la loi de paire la sélectionne fortement. |
| La fusion critique a les liens de meilleure qualité. | Vrai parmi les temps postcritiques, à coupe fixée. | Faux parmi les fusions déjà présentes à $`\beta_c`$ ; les fusions plus précoces ont un $`h_p`$ plus grand. |
| Tous les ancêtres ont la qualité critique. | Faux : seul $`\beta_v>\beta_u`$ est automatique. | Les ancêtres de la fenêtre gauche peuvent avoir $`h_p(\beta_v)\ge h_c`$ ; ceux après le seuil ont $`h_p(\beta_v)\le h_c`$. |
| Une majorité conforme dans chaque coupe suffit. | Une majorité globale par ancêtre ne suffit pas. | La majorité locale stricte et les deux majorités **groupées** non négatives suffisent sous a priori uniforme, par le théorème 7.1. |
| Conditionner par la paire critique donne une borne globale. | Faux sans le facteur de masse dans (8.2). | La domination HF uniforme est en outre réfutée en multiport ; une comparaison cible-spécifique exige une loi de bord contrôlée. |
| Le seuil de majorité est un seuil de weak recovery. | Faux. | Les deux seuils de majorité sont locaux et inférieurs à $`p_{\mathrm{info}}`$. |

## 10. Verrous mathématiques désormais isolés

Sous la loi de Palm (4.2), l'objet géométrique complet est

```math
\left(
\Pi_{\beta_c},C,I_L,J_L,u,
\left(
\beta_v,m_{v,0},m_{v,1},m_{v,2},G_v
\right)_{v\succeq u}
\right).
\qquad\text{(10.1)}
```

Conditionnellement à (10.1), les comptes $`K_{v,r}`$ sont donnés par (6.2),
les quatre taux par (6.3), les poids par (7.1) et le bit gardé par (7.2). La
loi des marques n'est donc plus le verrou principal.

Les difficultés restantes sont les suivantes.

1. **Exhaustion et événement macroscopique.** Fixer tore, boîte ou cylindre,
   métrique, conditions de bord, puis conserver la même définition de
   $`N_{L,\rho}`$ dans toutes les limites.
2. **Normalisation de Palm.** Estimer (4.3) et la loi de la composante sous le
   biais de paires. Remplacer ce biais par $`|C|`$, $`|C|^2`$ ou « plus grande
   composante » est une modification du problème.
3. **Fenêtre du LCA.** La localisation $`o(1)`$ de (5.4) est suffisante pour
   passer $`h_p(\beta_u)`$ à $`h_c(p)`$, mais pas pour décrire les nombres
   d'ancêtres et les interfaces dans la fenêtre. Il faut des estimations
   pivotales et multi-bras adaptées à la percolation de liens.
4. **Coupe du LCA.** Déterminer la loi de $`m_u=|E(C_1,C_2)|`$ sous le biais du
   nombre de paires lointaines séparées par $`C_1`$ et $`C_2`$. Ce biais est
   $`|C_1||C_2|`$ seulement sans contrainte de distance.
5. **Chaîne ancestrale.** Contrôler conjointement les trois tailles
   $`m_{v,r}`$, les temps $`\beta_v`$, le nombre d'ancêtres et les baselines
   $`K_{v,0}`$. Une estimation marginale de chaque $`\Lambda_v`$ ne contrôle
   pas le produit non linéaire.
6. **Coins nuls et queues.** Traiter exactement les états
   $`\Lambda_v^{ab}=0`$, puis montrer la sommabilité de la queue des facteurs
   dont $`m_{v,r}h_p(\beta_v)^2`$ reste petit.
7. **Transfert aux rangs réalisés.** Construire directement l'expérience
   multiport avec son état de bord. La domination HF uniforme étant fausse,
   une comparaison cible-spécifique ne peut servir qu'après contrôle explicite
   de cette loi de bord.
8. **Retour au score global.** Réinsérer la masse (8.2) et la composition du
   transfert réel avant toute conclusion de weak recovery.
9. **Seuil numérique.** Calculer le critère quatre états, puis seulement
   chercher une racine en $p$. Ni $`p_{\partial,\mathrm{late}}`$ ni
   $`p_{\mathrm N}^{(0)}`$ ne peut être déclaré seuil global par analogie.

### Programme de calcul certifiable

Historiquement, le premier calcul proposé était un cactus de triangles, puis
une bande triangulaire de largeur fixée. Le cactus reste un test unitaire ; la
bande n'est active qu'après construction d'une jauge de ports Markov-fermée.
Sans cette fermeture, l'état fidèle donne un déficit local nul. Sous **une
même loi de paire**, un éventuel calcul de bande doit comparer :

1. le critère exact (7.5) ;
2. le certificat groupé du théorème 7.1 ;
3. le diagnostic de frontière (3.8) ;
4. la même expérience sans biais de paire, afin de mesurer explicitement le
   changement de loi.

Le module
[`computations/critical_component_boundary.py`](../../computations/critical_component_boundary.py)
vérifie les identités scalaires, la loi binomiale de coupe et les quatre taux
groupés. Son test unitaire ajoute trois contre-audits finis : annulation des
arêtes internes, factorisation conditionnelle des seules arêtes de frontière
sur un triangle, et biais exact de sélection d'une composante par les paires.

## 11. Références primaires ciblées

- [Wierman, *Bond percolation on honeycomb and triangular lattices*](https://doi.org/10.2307/1426685) : valeur exacte de $`q_c`$ pour la percolation de liens triangulaire.
- [Grimmett--Manolescu, *Inhomogeneous bond percolation on square, triangular and hexagonal lattices*](https://doi.org/10.1214/11-AOP729) : transformation étoile--triangle et propriété de box-crossing.
- [Duminil-Copin--Tassion, *A new proof of the sharpness of the phase transition for Bernoulli percolation*](https://arxiv.org/abs/1502.03050) : décroissance exponentielle sous-critique.
- [Duminil-Copin--Sidoravicius--Tassion, *Absence of infinite cluster for critical Bernoulli percolation on slabs*](https://arxiv.org/abs/1401.7130) : absence de composante infinie au seuil dans le cadre planaire couvert.
- [Garban--Pete--Schramm, *The scaling limits of near-critical and dynamical percolation*](https://arxiv.org/abs/1305.5526) : fenêtre proche-critique et mesures pivotales pour la percolation par sites triangulaire ; cette référence ne justifie pas à elle seule le transfert au modèle de liens.
