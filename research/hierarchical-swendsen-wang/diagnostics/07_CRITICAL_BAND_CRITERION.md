# Bande critique et transmission de l'information

Ce fichier examine l'intuition suivante : une corrélation longue portée ne peut subsister que si la partie du dendrogramme située au-dessus du seuil de percolation transmet encore de l'information entre les deux sommets.

Le verdict est précis :

- la version « un chemin formé uniquement d'arêtes de niveau $(\beta,1]$ » est trop forte et n'est pas le bon événement ;
- la version après contraction des composantes de $`\Pi_\beta`$ est géométriquement exacte, mais redonne seulement la percolation Swendsen--Wang si elle n'est pas pondérée ;
- la version rigoureuse utile est une **naissance de connexion dans la bande**, pondérée par $`\eta_u`$, puis corrigée par la cohérence du signe du score LCA.

Pour le score signé de bande, la décomposition obtenue est

```math
Q_n^{>\beta}
=
\text{connexion quotient}
\times
\text{fiabilité locale}
\times
\text{cohérence signée}.
```

Lorsque $`S_n(\beta)\to0`$, ce score est équivalent au critère exact $`Q_n`$. Les deux premiers facteurs donnent une condition nécessaire explicite mais oracle. Le troisième est le verrou de la suffisance.

Dans tout ce fichier, l'a priori est binaire i.i.d. uniforme et les racines distinctes sont recolorées indépendamment et uniformément. Pour un potentiel $`\mu_0`$ général, une corrélation inter-racines peut subsister : il faut alors ajouter son terme aux formules de bande.

## 1. Filtration, temps de coalescence et bon quantificateur

Posons

```math
G_t=(V,E_t),
\qquad
E_t=\{e\in E:\xi_e\le t\},
\qquad
0\le t\le1,
```

et notons $`\Pi_t`$ la partition en composantes de $`G_t`$. Pour $i\ne j$,

```math
\beta_{ij}
=
\inf\{t:i\leftrightarrow j\text{ dans }G_t\}
=
\min_{\gamma:i\leadsto j}\max_{e\in\gamma}\xi_e.
```

Si $i,j$ ne sont pas reliés avant $1$, on pose $`\beta_{ij}=+\infty`$. Lorsque le LCA existe,

```math
\beta_{ij}=\beta_{u_{ij}}.
```

Avec des horloges continues, l'événement $`\beta_{ij}=\beta`$ a probabilité nulle pour un niveau déterministe $\beta$. Il faut donc employer :

- une bande $`a<\beta_{ij}\le b`$ ;
- une fenêtre critique $`|\beta_{ij}-\beta_c|\le\delta_n`$ ;
- ou la mesure complète des temps de fusion.

De même, « probabilité strictement positive » est trivial à volume fini dès qu'un chemin admissible existe. Le quantificateur pertinent est une masse macroscopique de paires. Pour $`I_n,J_n`$ uniformes et indépendants dans $V_n$, posons

```math
S_n(t)
=
\mathbb P(I_n\leftrightarrow J_n\text{ dans }G_t)
=
\frac1{n^2}
\mathbb E\sum_{C\in\Pi_t}|C|^2.
```

La condition $`S_n(t)\to0`$ signifie qu'aucune masse macroscopique de paires n'est déjà reliée au niveau $t$.

## 2. Cinq événements à ne pas confondre

Fixons un niveau déterministe $\beta$.

### 2.1 Connexion Swendsen--Wang

```math
\mathcal A_{ij}^{\mathrm{SW}}
=
\{\beta_{ij}\le1\}.
```

C'est l'événement où $i,j$ appartiennent à la même racine de $D$.

### 2.2 Naissance dans la bande

```math
\mathcal A_{ij}^{\mathrm{birth}}(\beta)
=
\{\beta<\beta_{ij}\le1\}.
```

La paire est séparée au niveau $\beta$, puis devient connectée avant la coupe $1$.

### 2.3 Connexion par sprinkling après contraction

Contractons chaque composante de $`\Pi_\beta`$ en un sommet. Ajoutons ensuite les arêtes telles que $`\beta<\xi_e\le1`$. Déterministiquement,

```math
\mathcal A_{ij}^{\mathrm{birth}}(\beta)
=
\left\{
[i]_\beta\ne[j]_\beta
\text{ et }
[i]_\beta\leftrightarrow[j]_\beta
\text{ dans le graphe quotient sprinklé}
\right\}.
```

Les chemins correspondants alternent : déplacements gratuits dans une composante précoce, puis arête tardive entre deux composantes. C'est la bonne lecture géométrique de « les arêtes de la bande connectent $i$ et $j$ ».

### 2.4 Chemin pur dans la bande

```math
\mathcal A_{ij}^{\mathrm{pure}}(\beta)
=
\left\{
i\not\leftrightarrow j\text{ dans }G_\beta,
\text{ et }
i\leftrightarrow j
\text{ par un chemin dont chaque arête vérifie }
\beta<\xi_e\le1
\right\}.
```

Cet événement interdit d'utiliser les composantes déjà construites sous $\beta$. Conditionnellement à la séparation au niveau $\beta$, il est plus exigeant que la naissance dans la bande et n'est pas un critère nécessaire universel de corrélation.

Un schéma de contre-exemple prend $`n=B_nL_n`$ sommets, avec $`B_n,L_n\to\infty`$, répartis en blocs de taille $`L_n`$. Des couplages internes forts alignent chaque bloc et le connectent avant $\beta$ avec probabilité $`1-o(1)`$. Un représentant par bloc porte ensuite un modèle de synchronisation récupérable dont les arêtes utiles ouvrent après $\beta$. Le quotient récupère les $`B_n`$ orientations et les propage aux blocs, tandis que le graphe pur tardif ne touche qu'une proportion $`O(L_n^{-1})`$ des sommets. Il faut formaliser la séparation des échelles de poids, mais cet exemple montre pourquoi la contraction des blocs appartient à l'énoncé même.

### 2.5 Redondance autour de l'arête de Kruskal

Au nœud $u$, on peut supprimer l'arête gagnante $e_u$, contracter les composantes de $`\Pi_{\beta_u^-}`$, puis demander si d'autres arêtes tardives reconnectent les deux fils. Cet événement mesure les cycles et les chemins alternatifs autour de la fusion. Il peut renforcer le LLR via tous les liens de $`E_u`$, mais sa simple positivité reste purement géométrique.

## 3. Théorème de réduction à la bande critique

Conservons le score LCA du fichier précédent :

```math
g_{ij}=\mathsf H_{ij}(\sigma_i\sigma_j),
\qquad
\eta_{ij}^{\mathrm{LCA}}=g_{ij}^2.
```

Définissons la masse informative des fusions nées dans une bande :

```math
\mathcal M_n((a,b])
=
\frac2{n^2}
\mathbb E\left[
\sum_{u:\,a<\beta_u\le b}
|C_{u,1}|\,|C_{u,2}|\,\eta_u
\right].
```

### Proposition — statut : établi conditionnellement à A1

Pour tout niveau déterministe $`0\le\beta\le1`$,

```math
\boxed{
Q_n
\le
S_n(\beta)
+
\mathcal M_n((\beta,1]).
}
```

En particulier, si

```math
S_n(\beta)\longrightarrow0,
```

alors la weak recovery impose

```math
\liminf_{n\to\infty}
\mathcal M_n((\beta,1])>0.
```

### Preuve

La borne LCA donne

```math
Q_n
\le
\frac1n
+
\frac2{n^2}
\mathbb E\sum_{u:\,\beta_u\le\beta}
|C_{u,1}||C_{u,2}|\eta_u
+
\mathcal M_n((\beta,1]).
```

Comme $`0\le\eta_u\le1`$, les deux premiers termes sont dominés par

```math
\frac1n
+
\frac2{n^2}
\mathbb E\sum_{u:\,\beta_u\le\beta}
|C_{u,1}||C_{u,2}|
=S_n(\beta).
```

### Conséquence événementielle quantitative

Supposons $`Q_n\ge q`$ et $`S_n(\beta)\le\varepsilon<q`$. Pour $`0<a<q-\varepsilon`$,

```math
\mathbb P\left(
\beta<\beta_{I_nJ_n}\le1,
\ \eta_{I_nJ_n}^{\mathrm{LCA}}\ge a
\right)
\ge
\frac{q-\varepsilon-a}{1-a}.
```

La weak recovery exige donc une masse positive de paires qui :

1. naissent au-dessus de la coupe critique ;
2. ont une fiabilité LCA non négligeable.

C'est une version rigoureuse de l'intuition proposée.

## 4. Critère nécessaire et suffisant : le troisième facteur

La masse $`\mathcal M_n`$ reste une borne oracle, car $D$ est révélé dans le heat bath LCA. Pour isoler exactement la perte restante, posons

```math
Z_{ij}^{>\beta}
=
g_{ij}\,
\mathbf1_{\{\beta<\beta_{ij}\le1\}},
```

avec la convention $`Z_{ii}^{>\beta}=0`$, et

```math
Q_n^{>\beta}
=
\frac1{n^2}
\sum_{i,j}
\mathbb E_O\left[
\left(
\mathbb E_{\nu_O}[Z_{ij}^{>\beta}]
\right)^2
\right].
```

Par Jensen,

```math
Q_n^{>\beta}
\le
H_n^{>\beta}
:=
\frac1{n^2}
\sum_{i,j}
\mathbb E[(Z_{ij}^{>\beta})^2]
=
\mathcal M_n((\beta,1]).
```

La partie précoce est dominée par $`S_n(\beta)`$. Une application de Cauchy--Schwarz donne

```math
\left|Q_n-Q_n^{>\beta}\right|
\le
S_n(\beta)+2\sqrt{S_n(\beta)}.
```

En effet, posons $`Y_{ij}=g_{ij}-Z_{ij}^{>\beta}`$, avec $`Y_{ii}=1`$. Alors

```math
\frac1{n^2}\sum_{i,j}
\mathbb E[(Y_{ij})^2]
\le S_n(\beta),
```

et la différence des carrés des moyennes est contrôlée par Cauchy--Schwarz. Cette convention incorpore exactement la diagonale $`1/n`$.

Donc, sous $`S_n(\beta)\to0`$,

```math
\boxed{
\liminf Q_n>0
\quad\Longleftrightarrow\quad
\liminf Q_n^{>\beta}>0.
}
```

Ce critère de bande est nécessaire et suffisant. Il emploie toutefois le **score signé**, pas seulement la connectivité.

### Factorisation en trois verrous

Posons

```math
R_n^{>\beta}
=
\mathbb P(\beta<\beta_{I_nJ_n}\le1),
```

puis, lorsque les dénominateurs sont non nuls,

```math
\overline\eta_n^{>\beta}
=
\frac{H_n^{>\beta}}{R_n^{>\beta}},
\qquad
\kappa_n^{>\beta}
=
\frac{Q_n^{>\beta}}{H_n^{>\beta}}.
```

On a exactement

```math
\boxed{
Q_n^{>\beta}
=
R_n^{>\beta}
\,\overline\eta_n^{>\beta}
\,\kappa_n^{>\beta}.
}
```

- $`R_n^{>\beta}`$ est la naissance géométrique dans le quotient sprinklé ;
- $`\overline\eta_n^{>\beta}`$ est la fiabilité locale moyenne des fusions ;
- $`\kappa_n^{>\beta}`$ est la cohérence signée après marginalisation de $D$.

Le défaut de cohérence vaut

```math
H_n^{>\beta}-Q_n^{>\beta}
=
\frac1{n^2}
\sum_{i,j}
\mathbb E_O
\mathrm{Var}_{\nu_O}(Z_{ij}^{>\beta}).
```

La conjecture fondée sur la seule connexion de bande ne conserve que le premier facteur.

## 5. Réinterprétation dendrogramme de l'information-percolation

Dans le GSBM triangulaire homogène, posons

```math
u_p=\log\frac p{1-p},
\qquad
q_p(t)=p(1-e^{-u_pt}),
\qquad
q_p(1)=2p-1.
```

Notons

```math
q_c=2\sin(\pi/18)
=0.3472963553\ldots
```

le seuil de percolation par arêtes de la grille triangulaire. Le niveau géométrique critique est

```math
\beta_c(p)
=q_p^{-1}(q_c)
=-
\frac1{u_p}
\log\left(1-\frac{q_c}{p}\right).
```

Il appartient à $`[0,1]`$ exactement lorsque

```math
p\ge p_{\mathrm{SW}}
:=
\frac{1+q_c}{2}
=0.6736481777\ldots
```

La contraction $`\chi^2`$ du canal binaire d'une arête vaut

```math
\gamma_p=(2p-1)^2.
```

Introduisons le **temps informationnel équivalent**

```math
t_\chi(p)
=q_p^{-1}(\gamma_p)
=-
\frac1{u_p}
\log\left(
1-\frac{(2p-1)^2}{p}
\right).
```

Le graphe de Bernoulli de paramètre $`\gamma_p`$ utilisé par information-percolation a la même loi marginale que $`G_{t_\chi}`$. C'est un couplage de la borne informationnelle avec l'échelle des horloges, pas l'affirmation que $D$ devient une nouvelle observation.

Ici, par définition, $`q_p(t_\chi)=\gamma_p`$ ; ce n'est pas la convention différente $`q_p(t_\chi)^2=q_c`$.

Lorsque $`t_\chi\ge\beta_c`$, après contraction de $`\Pi_{\beta_c}`$, toute connexion nouvellement créée dans $`G_{t_\chi}`$ utilise exactement les composantes critiques et les arêtes de la sous-bande

```math
\beta_c<\xi_e\le t_\chi.
```

Le théorème d'Abbe--Boix donne la condition nécessaire

```math
t_\chi(p)>\beta_c(p)
\quad\Longleftrightarrow\quad
(2p-1)^2>q_c.
```

Ainsi,

```math
\boxed{
p>
p_{\mathrm{info}}
:=
\frac{1+\sqrt{q_c}}2
=0.7946592758\ldots
}
```

est nécessaire à la weak recovery. Au seuil,

```math
t_\chi(p_{\mathrm{info}})
=
\beta_c(p_{\mathrm{info}})
=0.4245677743\ldots
```

Cette identité confirme l'intuition d'une bande au-dessus du seuil, mais avec trois corrections : contraction préalable, borne supérieure $`t_\chi`$ plutôt que $1$, et caractère seulement nécessaire.

## 6. Trois graphes de bande, trois seuils triangulaires

Au niveau $`\beta_c`$, la densité marginale des seules arêtes tardives est

```math
r_p
=
q_p(1)-q_c
=
2p-1-q_c.
```

| Objet | Condition de percolation | Seuil en $p$ | Portée |
|---|---:|---:|---|
| Quotient de $`\Pi_{\beta_c}`$ avec toute la bande jusqu'à $1$ | $`q_p(1)>q_c`$ | $`0.673648\ldots`$ | exactement SW, trop faible |
| Quotient informationnel jusqu'à $`t_\chi`$ | $`(2p-1)^2>q_c`$ | $`0.794659\ldots`$ | borne information-percolation de référence |
| Graphe formé uniquement des arêtes $`\beta_c<\xi_e\le1`$ | $`r_p>q_c`$ | $`0.847296\ldots`$ | trop exigeant, non nécessaire en général |

Le dernier seuil vaut

```math
p_{\mathrm{pure}}
=
\frac12+q_c
=0.8472963553\ldots
```

Il est supérieur au point multicritique de Nishimori numérique $`0.835806\ldots`$. Ce dernier n'étant pas un seuil rigoureusement établi dans le présent dossier, la comparaison est un diagnostic fort, pas une réfutation autonome sur la grille triangulaire.

## 7. Ce que les arêtes de bande apportent à une coupe

Considérons un bucket de $m$ liens, fusionnant au temps $\beta$, conditionnellement au squelette de Kruskal non marqué. Une arête gagnante latente est uniforme ; conditionnellement au fait qu'un autre lien soit encore fermé à $\beta$, posons

```math
h_p(\beta)
=
\mathbb P(\beta<\xi_e\le1\mid\xi_e>\beta)
=
\frac{q_p(1)-q_p(\beta)}{1-q_p(\beta)}
=
\tanh\left(\frac{u_p(1-\beta)}2\right).
```

Les $m-1$ autres liens se répartissent en :

- $R$ : satisfaits et dans la bande $(\beta,1]$ ;
- $S$ : satisfaits mais avec $`\xi_e>1`$ ;
- $U$ : insatisfaits.

Conditionnellement à ce squelette, y compris lorsque le bucket est sélectionné par Kruskal, on a exactement

```math
(R,S,U)
\sim
\mathrm{Mult}\left(
m-1;
h_p(\beta),
\frac{1-h_p(\beta)}2,
\frac{1-h_p(\beta)}2
\right).
```

Le nombre total de liens satisfaits est

```math
k=1+R+S.
```

Conditionnellement à $R=r$,

```math
k
\ \stackrel d=\
1+r+\mathrm{Bin}(m-1-r,1/2).
```

Les arêtes de bande constituent donc exactement la partie biaisée du vote. Les variables $S$ et $U$ forment le bruit symétrique.

### Connexion contre information

L'existence d'au moins une arête redondante de bande est gouvernée par

```math
\mathbb P(R\ge1)
=
1-(1-h_p(\beta))^{m-1},
```

donc par l'échelle $`m h_p(\beta)`$. En revanche, le rapport signal sur bruit du vote est gouverné par

```math
m h_p(\beta)^2.
```

Plus précisément, pour $`m\ge2`$ et $`X=2k-m`$,

```math
\frac{(\mathbb E X)^2}{\mathrm{Var}(X)}
=
\frac{[1+(m-1)h_p(\beta)]^2}
{(m-1)(1-h_p(\beta)^2)}.
```

Lorsque $`m h_p(\beta)\to\infty`$, l'échelle dominante est donc $`m h_p(\beta)^2`$, à facteur $`1/(1-h_p(\beta)^2)`$ près. Si $`h_p(\beta)`$ est plus petit, le terme résiduel dû à l'arête gagnante reste visible ; à $`h_p(\beta)=0`$, son rapport signal sur bruit vaut $`1/(m-1)`$.

Une simple connexion apparaît donc bien avant une transmission statistiquement stable. C'est le pendant local du passage d'une probabilité d'ouverture à une contraction $`\chi^2`$.

### Fenêtre terminale d'une grande coupe

Dans le modèle local $`B_u=0`$,

```math
L_{m,k,\beta}^{\mathrm{loc}}
=
\log\frac{k}{m-k}
+
u_p(1-\beta)(2k-m).
```

Pour $\beta<1$ fixé, la fiabilité tend vers $1$ lorsque $`m\to\infty`$. La fenêtre où une grande coupe reste non triviale est

```math
1-\beta\asymp m^{-1/2}.
```

Plus précisément, si $`u_p(1-\beta)=a/\sqrt m`$, alors

```math
\frac{2k-m}{\sqrt m}
\Longrightarrow
Z+\frac a2,
\qquad
L_{m,k,\beta}^{\mathrm{loc}}
\Longrightarrow
aZ+\frac{a^2}{2},
```

avec $`Z\sim\mathcal N(0,1)`$. Ainsi,

```math
\mathbb E\eta_u
\longrightarrow
\mathbb E\left[
\tanh^2\left(
\frac{aZ+a^2/2}{2}
\right)
\right].
```

À $`\beta=1`$, on retrouve $`\mathbb E\eta_u=1/m`$.

Le [fichier 09](09_CRITICAL_MERGER_ORACLE.md) applique cette fenêtre au temps
$`\beta_c(p)`$. Il obtient les identités exactes

```math
h_p(\beta_c)
=
\frac{2(p-p_{\mathrm{SW}})}{1-q_c},
\qquad
u_p(1-\beta_c)
=
2\,\mathrm{artanh}(h_p(\beta_c)),
```

puis la limite

```math
\Gamma_m^c\left(
p_{\mathrm{SW}}+\frac{(1-q_c)\alpha}{2\sqrt m}
\right)
\longrightarrow
\mathbb E\left[\tanh^2(\alpha Z+\alpha^2)\right].
```

Cette calibration reste locale et oracle ; elle doit être multipliée par la
masse de la fenêtre parmi les paires.

## 8. Attention : $\eta_u$ est une quantité oracle

Conditionnellement au dendrogramme révélé, une coupe à une seule arête a $`\eta_u=1`$ dès qu'elle fusionne. Pourtant, la contraction informationnelle réelle d'un canal arête vaut

```math
(2p-1)^2,
```

et non $`2p-1`$. C'est exactement le gap déjà visible sur un chemin :

```math
A_{ij}^{(1)}=(2p-1)^\ell,
\qquad
c_{ij}^2=(2p-1)^{2\ell}.
```

Il ne faut donc pas utiliser directement $`\eta_u`$ comme probabilité d'une percolation suffisante. Il faut marginaliser ou rafraîchir $D$, ou contrôler le facteur $`\kappa_n`$.

Pour $m$ observations BSC indépendantes d'une même parité, la contraction de bloc non oracle est

```math
\gamma_m^{\mathrm{BSC}}
=
\sum_{k=0}^m
\binom mk p^k(1-p)^{m-k}
\tanh^2\left(
\frac{u_p}{2}(2k-m)
\right).
```

Par l'identité de Nishimori, la même quantité s'écrit aussi sans le carré, comme espérance du biais postérieur sous la vérité fixée à $`+1`$. Elle vérifie

```math
\gamma_1^{\mathrm{BSC}}=(2p-1)^2
```

et fournit la calibration obligatoire de toute capacité de quotient. Les coupes choisies par Kruskal ne sont toutefois pas des blocs déterministes indépendants : leur contraction conditionnelle reste à dériver.

## 9. Sprinkling pivotal pondéré

Dans le couplage annealed homogène,

```math
\mathbb P(\beta_{ij}\le t)
=
\tau_{ij}(q_p(t)).
```

Sur un graphe fini, la formule de Russo donne

```math
\frac d{dt}
\tau_{ij}(q_p(t))
=
q_p'(t)
\sum_{e\in E}
\mathbb P_{q_p(t)}
\left(
e\text{ pivotal pour }\{i\leftrightarrow j\}
\right),
```

avec

```math
q_p'(t)=p u_p e^{-u_pt}.
```

Par conséquent,

```math
\mathbb P(\beta<\beta_{ij}\le1)
=
\int_\beta^1
q_p'(t)
\sum_e
\mathbb P_{q_p(t)}(e\text{ pivotal})
\,dt.
```

En désintégrant

```math
h_{ij}(t)
=
\mathbb E[\eta_{u_{ij}}\mid\beta_{ij}=t],
```

la borne LCA devient formellement

```math
\mathbb E[c_{ij}(O)^2]
\le
\int_0^1
h_{ij}(t)
\,d\tau_{ij}(q_p(t)).
```

La version au-dessus de la coupe critique est l'intégrale sur $`(\beta_c,1]`$. Cet objet est un **flux pivotal informatif** : fréquence des fusions pivotales, taille des blocs fusionnés et fiabilité de leur canal.

La conditionnelle $`h_{ij}(t)`$ incorpore le biais du LCA sélectionné. Elle ne peut pas être remplacée sans preuve par la loi d'une coupe déterministe.

## 10. Pourquoi les LCA ne sont pas tous exactement critiques

Pour des sommets lointains dans un modèle homogène,

```math
\mathbb P(\beta_{ij}\le t)
=
\tau_{ij}(q_p(t)).
```

Dans un régime supercritique avec unicité et mélange spatial, la limite attendue pour des sommets qui s'éloignent est

```math
\tau_{ij}(q_p(t))
\longrightarrow
\theta(q_p(t))^2,
\qquad t>\beta_c.
```

Conditionnellement à une connexion avant $1$, cela produit une distribution généralement non dégénérée au-dessus de $`\beta_c`$, et non un atome au seuil. Les attaches locales peuvent aussi introduire des goulots au-dessus de la fenêtre critique.

Il existe en outre une obstruction de masse entièrement finie. Pour
$`\mathcal A_{n,\delta}^c=\{\beta_c<\beta_{I_nJ_n}\le\beta_c+\delta\}`$,

```math
\mathbb E\left[
\eta_{I_nJ_n}^{\mathrm{LCA}}
\mathbf1_{\mathcal A_{n,\delta}^c}
\right]
\le
S_n(\beta_c+\delta)-S_n(\beta_c).
```

Ainsi, une fiabilité conditionnelle parfaite dans une fenêtre où le membre de
droite tend vers zéro ne peut pas soutenir un overlap macroscopique. La preuve
et les conséquences sont détaillées dans le
[contre-audit de l'oracle critique](09_CRITICAL_MERGER_ORACLE.md).

Les résultats de percolation proche-critique et de MSF indiquent plutôt que la **géométrie macroscopique** est codée par les pivots proche-critiques. La formulation correcte est donc :

1. contracter les détails microscopiques et les attaches locales ;
2. étudier une fenêtre autour de $`\beta_c`$ pour les goulots macroscopiques ;
3. conserver toute la mesure $`\mathcal M_n(dt)`$ pour les paires ponctuelles.

Le [fichier 12](../archive/roadmaps/12_FAVORABLE_HIERARCHICAL_REDUCTION.md) formalise aussi une
ancienne réduction conditionnelle HF : si le log-rapport complet de toute
paire postcritique était dominé par celui de l'oracle critique, l'intégrale
ci-dessus serait majorée par la seule fiabilité critique. Le contre-exemple
multiport du [fichier 29](29_AUDIT_FROID_PIVOT_RANGS_REELS.md) montre que cet
ordre **uniforme est faux**. Seule une comparaison cible-spécifique sous une
loi de bord déjà contrôlée peut encore être recherchée ; elle ne remplace pas
le transfert direct aux rangs réalisés.

### Échelle de fenêtre à tester en volume fini

Sur un tore de diamètre linéaire $L$, la théorie proche-critique de la percolation de sites triangulaire utilise l'échelle

```math
q-q_c\asymp L^{-3/4}.
```

Transposée comme guide de calcul au couplage par arêtes présent, elle suggère

```math
t_L(\lambda)
=
\beta_c
+
\frac{\lambda}{q_p'(\beta_c)L^{3/4}},
\qquad
q_p'(\beta_c)=p u_p e^{-u_p\beta_c}.
```

Le nombre de pivots macroscopiques est alors d'ordre $`L^{3/4}`$, de sorte qu'un sprinkling de largeur $`L^{-3/4}`$ produit un nombre d'événements macroscopiques d'ordre un. Cette échelle est rigoureuse dans le cadre proche-critique de sites traité par Garban--Pete--Schramm ; son transfert aux horloges d'arêtes et à la pondération signée $`\eta_u`$ est un objectif, pas un résultat acquis ici.

## 11. Pistes de recherche maintenant bien séparées

### Piste A — fermer le théorème de bande

Formaliser la proposition de la section 3 avec l'exhaustion exacte de la grille et prouver $`S_n(\beta_c)\to0`$ sous les conditions de bord choisies.

### Piste B — capacité du quotient critique

Contracter $`\Pi_{\beta_c}`$, puis associer à chaque bundle entre deux blocs une contraction conditionnelle non oracle. La cible est une capacité $L^2$, un flot ou un rayon spectral qui :

- redonne $`(2p-1)^2`$ pour une arête ;
- utilise $`\gamma_m^{\mathrm{BSC}}`$ pour un bundle indépendant de taille $m$ ;
- intègre le message $`B_u`$ et les cycles pour une fusion réelle ;
- domine rigoureusement les corrélations après marginalisation de $D$.

### Piste C — fenêtre terminale et interfaces

Étudier conjointement la multiplicité $m$, le temps $t$ et le message $B$ sous la mesure pivotale de Kruskal. Le paramètre local naturel est $`m h_p(t)^2`$.

### Piste D — itérations et cohérence

Employer la chaîne pair-spécifique avec rafraîchissement de $D$ pour faire décroître la variance conditionnelle. Le but est de contrôler uniformément

```math
H_n^{(m_n)}-Q_n.
```

### Piste E — score signé calculable

Estimer

```math
\widehat c_{ij}^{>\beta}
=
\mathbf1_{\{\beta<\beta_{ij}\le1\}}
\sigma_i\sigma_jm_{u_{ij}},
```

puis agréger cette matrice implicitement. Une borne inférieure sur sa composante spectrale cohérente fournirait la partie suffisante.

### Piste F — enhancement essentiel

Si une contraction de bloc peut être encadrée par un enhancement local monotone, les théorèmes d'Aizenman--Grimmett peuvent servir à prouver un déplacement strict du seuil. La frustration, les signes et la sélection adaptative de Kruskal empêchent pour l'instant une application directe.

## 12. Statut des énoncés

| Énoncé | Statut |
|---|---|
| $`Q_n\le S_n(\beta)+\mathcal M_n((\beta,1])`$ | Établi conditionnellement à A1 |
| Réduction nécessaire et suffisante au score signé de bande lorsque $`S_n(\beta)\to0`$ | Établi conditionnellement à A1 |
| Connexion nouvellement créée dans le quotient par la bande équivaut à $`\beta<\beta_{ij}\le1`$ | Identité déterministe |
| Bande pure nécessaire à la weak recovery | Pas un critère universel ; la contraction préalable est indispensable |
| Réécriture $`t_\chi>\beta_c`$ de la borne information-percolation | Établi dans le modèle homogène annealed |
| Loi multinomiale conditionnelle au squelette non marqué | Établie ; le biais restant porte sur la géométrie du squelette |
| Fiabilité locale d'une fusion exactement critique | Établie ; bord oracle $`p_{\mathrm{SW}}`$ et fenêtre $`m^{-1/2}`$ |
| Masse informative d'une fenêtre critique | Majorée exactement par $`S_n(\beta_c+\delta)-S_n(\beta_c)`$ |
| Domination de toute l'expérience postcritique par l'oracle critique HF | Fausse en général pour une fusion multiport ; seule une comparaison cible-spécifique sous la vraie loi de bord reste ouverte |
| Capacité de quotient pondérée par une contraction non oracle | À construire |
| Seuil strictement supérieur à $`0.794659\ldots`$ | À prouver |

## 13. Références directement utiles

- [Abbe--Boix, information-percolation](https://arxiv.org/abs/1806.03227) : contraction $`\chi^2`$ et condition nécessaire.
- [Polyanskiy--Wu, méthode information-percolation](https://arxiv.org/abs/1806.04195) : SDPI et comparaison avec une percolation d'effacement.
- [Wierman, percolation par arêtes triangulaire](https://doi.org/10.2307/1426685) : valeur exacte $`q_c=2\sin(\pi/18)`$.
- [Garban--Pete--Schramm, percolation proche-critique](https://arxiv.org/abs/1305.5526) : fusions par pivots dans la fenêtre critique.
- [Garban--Pete--Schramm, mesures pivotales](https://arxiv.org/abs/1008.1378) : mesure des pivots macroscopiques.
- [Nolin, percolation proche-critique](https://arxiv.org/abs/0711.4948) : longueur de corrélation et fenêtre critique bidimensionnelle.
- [Garban--Pete--Schramm, MST et invasion](https://arxiv.org/abs/1309.0269) : géométrie macroscopique de la MSF proche-critique.
- [Damron--Sapozhnikov, outlets d'invasion](https://arxiv.org/abs/0903.4496) : distinction entre attaches locales et approche du seuil critique.
- [Evans--Kenyon--Peres--Schulman](https://doi.org/10.1214/aoap/1019487349) : capacité $L^2$ et reconstruction sur arbres.
- [Aizenman--Grimmett](https://doi.org/10.1007/BF01029985) : déplacement strict sous enhancement essentiel.

La prochaine étape mathématique n'est donc pas de tester une connectivité
brute supplémentaire. Pour la voie hiérarchique, elle consiste à calculer le
transfert multiport $`K,U`$ avec tous les $`\Lambda_v`$ aux rangs réalisés,
puis à prouver une décroissance annealed de Feynman--Kac. La bande critique
reste un benchmark local, pas une enveloppe globale.
