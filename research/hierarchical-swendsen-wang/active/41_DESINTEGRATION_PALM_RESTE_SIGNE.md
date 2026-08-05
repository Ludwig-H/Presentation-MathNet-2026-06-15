# Désintégration exacte du reste signé : moyenner en $D$ avant le carré

**Statut : identités exactes en volume fini ; Palm géométrique et
normalisation du scalaire fermées ; embedding récursif des ports et
contraction encore ouverts ; aucune nouvelle borne de weak recovery.**

Cette note simplifie la
[cible double-géante](38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md). Le reste signé
avait été écrit avec deux dendrogrammes et le raffinement commun de leurs
deux coupes critiques. Cette représentation est exacte, mais elle masque une
désintégration plus simple.

Le verdict est :

```math
\boxed{
\text{moyenner le transfert signé en }D
\quad\text{puis}\quad
\text{prendre son carré}.
}
\qquad\text{(0.1)}
```

À une erreur contrôlée par la masse quadratique des blocs critiques, il
n'est pas nécessaire de construire une Palm jointe abstraite. La bonne Palm
est la loi de **deux corridors cross-block indépendants conditionnellement à
la même observation et aux mêmes endpoints**.

Cette réduction ferme la porte TRI1 au niveau scalaire et géométrique. Elle
ne ferme pas TRI2 : pour composer les transferts, il faut encore plonger les
espaces de ports dépendant de $D$ dans des coordonnées communes et contrôler
la norme $L^2$ avant toute fermeture non linéaire.

## 1. Trois fonctions d'un dendrogramme

Fixons une observation $O$ et une paire ordonnée $`(i,j)`$. Sous
$`D\sim\rho_O`$, posons

```math
m_D(i,j)
=
\pi_{O,D}(\sigma_i\sigma_j).
\qquad\text{(1.1)}
```

Soit $`R_\star(D)`$ la plus grande racine finale, avec une règle de départage
déterministe, et $`A_D(i)`$ le bloc de $i$ à la coupe critique. Définissons

```math
g_D(i,j)
=
\mathbf1_{\{i,j\in R_\star(D)\}},
\qquad
s_D(i,j)
=
\mathbf1_{\{A_D(i)=A_D(j)\}}.
\qquad\text{(1.2)}
```

Les trois moyennes single-$D$ pertinentes sont

```math
\begin{aligned}
a_{O,ij}
&=
\mathbb E_{D\mid O}[g_Dm_D],
\\
b_{O,ij}
&=
\mathbb E_{D\mid O}[g_Ds_Dm_D],
\\
d_{O,ij}
&=
\mathbb E_{D\mid O}[g_D(1-s_D)m_D]
=
a_{O,ij}-b_{O,ij}.
\end{aligned}
\qquad\text{(1.3)}
```

Ici $b$ porte les paires dans un même bloc critique et $d$ les paires dans
deux blocs critiques distincts d'un seul dendrogramme. Les signes de
$m_D$ sont conservés dans les trois définitions.

## 2. Différence de carrés exacte

Tirons $`D^{(1)},D^{(2)}`$ indépendamment sous $`\rho_O`$. Les endpoints
sont dans la double géante et dans deux cellules distinctes du raffinement
commun exactement sur l'événement

```math
g_1g_2(1-s_1s_2)=1.
\qquad\text{(2.1)}
```

L'indépendance conditionnelle des deux dendrogrammes donne immédiatement

```math
\begin{aligned}
&
\mathbb E_{D^{(1)},D^{(2)}\mid O}
\left[
g_1g_2(1-s_1s_2)m_1m_2
\right]
\\
&\hspace{2cm}
=
\left(
\mathbb E_{D\mid O}[g_Dm_D]
\right)^2
-
\left(
\mathbb E_{D\mid O}[g_Ds_Dm_D]
\right)^2
\\
&\hspace{2cm}
=
a_{O,ij}^2-b_{O,ij}^2
=
d_{O,ij}^2+2b_{O,ij}d_{O,ij}.
\end{aligned}
\qquad\text{(2.2)}
```

Par conséquent, le reste signé de la note 38 vaut exactement

```math
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
\left(
a_{O,ij}^2-b_{O,ij}^2
\right).
\qquad\text{(2.3)}
```

Cette identité explique simultanément deux faits.

1. Un tirage fini du reste à deux $D$ peut être négatif.
2. Le bon objet n'est pas la moyenne des carrés
   $`\mathbb E_D[m_D^2]`$, mais le carré de la moyenne signée
   $`(\mathbb E_D[m_D])^2`$.

Le premier objet est l'enveloppe single-$D$ de Jensen. Le second est la
quantité exacte de weak recovery.

## 3. Élimination quantitative du terme same-block

Posons la masse quadratique critique dans la grande racine

```math
S_L^c(p)
=
\frac1{n_L^2}
\mathbb E_{O,D}
\left[
\sum_{\substack{
A\in\Pi_{\beta_c}(D)\\
A\subseteq R_\star(D)
}}
|A|^2
\right].
\qquad\text{(3.1)}
```

Comme $`|m_D|\le1`$,

```math
|b_{O,ij}|
\le
\mathbb E_{D\mid O}[g_Ds_D].
\qquad\text{(3.2)}
```

La convexité, puis la somme sur les endpoints, donnent

```math
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
b_{O,ij}^2
\le
S_L^c(p).
\qquad\text{(3.3)}
```

Définissons maintenant le carré cross-block

```math
\mathcal D_L^\times(p)
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
d_{O,ij}^2
\ge0.
\qquad\text{(3.4)}
```

Puisque la même moyenne normalisée de $d_{O,ij}^2$ est au plus un,
Cauchy--Schwarz appliqué à (2.2) donne

```math
\left|
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
-
\mathcal D_L^\times(p)
\right|
\le
2\sqrt{S_L^c(p)}.
\qquad\text{(3.5)}
```

À la coupe de percolation critique, $`S_L^c(p)\to0`$. Avec la réduction hors
géante de la note 38,

```math
\boxed{
Q_L(p)\longrightarrow0
\quad\Longleftrightarrow\quad
\mathcal D_L^\times(p)\longrightarrow0.
}
\qquad\text{(3.6)}
```

Ainsi, les cancellations restent indispensables, mais elles peuvent être
prises **dans la moyenne single-$D$ $d_{O,ij}$**, avant de former un carré
positif. La valeur absolue ne doit toujours pas être appliquée avant cette
moyenne.

## 4. Palm cross--cross non circulaire

Pour $O,i,j$ fixés, posons la probabilité géométrique

```math
\delta_{O,ij}
=
\mathbb E_{D\mid O}
\left[
g_D(i,j)(1-s_D(i,j))
\right].
\qquad\text{(4.1)}
```

La masse annealed de deux corridors cross-block est

```math
\alpha_L^\times
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
\delta_{O,ij}^2.
\qquad\text{(4.2)}
```

Si $`I,J`$ sont uniformes et ordonnés dans $`V_L`$, cette quantité est
exactement la probabilité de l'événement

```math
\mathcal A_\times
=
\left\{
g_1(I,J)(1-s_1(I,J))
g_2(I,J)(1-s_2(I,J))
=1
\right\}.
\qquad\text{(4.2a)}
```

Si $`\alpha_L^\times=0`$, alors $`\mathcal D_L^\times=0`$ et aucune Palm
n'est requise. Sinon, la Palm $`\mathbb P_L^\times`$ se construit sans
utiliser $m_D$, $Q_L$ ni la conclusion recherchée :

1. biaiser $`(O,i,j)`$ par $`\delta_{O,ij}^2`$ ;
2. conditionnellement à $`(O,i,j)`$, tirer $`D^{(1)},D^{(2)}`$
   indépendamment selon

```math
\rho_{O,ij}^{\times}(dD)
=
\frac{
g_D(i,j)(1-s_D(i,j))
}{
\delta_{O,ij}
}
\rho_O(dD).
\qquad\text{(4.3)}
```

Les triplets $`(O,i,j)`$ tels que $`\delta_{O,ij}=0`$ ont une masse Palm
nulle ; la conditionnelle (4.3) peut y être fixée arbitrairement. Plus
explicitement, pour toute fonction bornée $F$,

```math
\mathbb E_{\mathbb P_L^\times}[F]
=
\frac{
\mathbb E\left[
g_1(I,J)(1-s_1(I,J))
g_2(I,J)(1-s_2(I,J))
F
\right]
}{
\alpha_L^\times
}.
\qquad\text{(4.3a)}
```

On a alors l'identité exacte

```math
\mathcal D_L^\times(p)
=
\alpha_L^\times
\mathbb E_{\mathbb P_L^\times}
\left[
m_{D^{(1)}}(I,J)m_{D^{(2)}}(I,J)
\right].
\qquad\text{(4.4)}
```

Les deux corridors ne sont donc pas indépendants sous la loi annealed
totale : ils partagent $O,i,j$. Ils sont toutefois réellement indépendants
conditionnellement à ces variables. C'est la factorisation utile pour TRI2.

## 5. La Palm du reste original comme audit

La Palm cross--cross suffit asymptotiquement grâce à (3.5). La Palm exacte
de l'événement original (2.1) reste utile pour auditer les normalisations.
Posons

```math
\gamma_{O,ij}
=
\mathbb E_{D\mid O}[g_D],
\qquad
\kappa_{O,ij}
=
\mathbb E_{D\mid O}[g_Ds_D].
\qquad\text{(5.1)}
```

Sa masse conditionnelle est

```math
Z_{O,ij}^{\mathrm{off}}
=
\gamma_{O,ij}^2-\kappa_{O,ij}^2.
\qquad\text{(5.2)}
```

Sa masse annealed, pour des endpoints uniformes ordonnés, est

```math
\alpha_L^{\mathrm{off}}
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
Z_{O,ij}^{\mathrm{off}}.
\qquad\text{(5.2a)}
```

Lorsque $`Z_{O,ij}^{\mathrm{off}}>0`$, notons
$`\rho^\times=\rho(\cdot\mid g(1-s)=1)`$ et
$`\rho^\circ=\rho(\cdot\mid gs=1)`$, avec les composantes de poids nul
simplement omises. La loi jointe conditionnelle est le mélange

```math
\frac{\gamma-\kappa}{\gamma+\kappa}
\rho^\times\otimes\rho^\times
+
\frac{\kappa}{\gamma+\kappa}
\left(
\rho^\times\otimes\rho^\circ
+
\rho^\circ\otimes\rho^\times
\right).
\qquad\text{(5.3)}
```

La loi annealed biaise ensuite $`(O,i,j)`$ par
$`Z_{O,ij}^{\mathrm{off}}`$. Moyenner uniformément les Palm conditionnelles
serait une erreur de normalisation. Lorsque
$`Z_{O,ij}^{\mathrm{off}}=0`$, la conditionnelle jointe peut être choisie
arbitrairement puisqu'elle reçoit une masse Palm nulle ; si
$`\alpha_L^{\mathrm{off}}=0`$, le reste original est déjà nul.

## 6. Ordre exact des moyennes pour l'opérateur

Supposons chaque transfert à $D$ fixé plongé dans des coordonnées physiques
communes et notons $`K_{O,D}^{\times}`$ son opérateur cross-block signé.
L'opérateur à moyenner est

```math
\overline K_O^\times
=
\mathbb E_{D\mid O}
\left[
K_{O,D}^{\times}
\right].
\qquad\text{(6.1)}
```

Le carré exact à observation fixée est

```math
\left(
\overline K_O^\times
\right)^{\!*}
\overline K_O^\times,
\qquad\text{puis on moyenne en }O.
\qquad\text{(6.2)}
```

Deux ordres incorrects doivent rester exclus :

```math
\mathbb E_{D\mid O}
\left[
\left(K_{O,D}^{\times}\right)^{\!*}
K_{O,D}^{\times}
\right]
\qquad\text{et}\qquad
\left(
\mathbb E_{O,D}K_{O,D}^{\times}
\right)^{\!*}
\left(
\mathbb E_{O,D}K_{O,D}^{\times}
\right).
\qquad\text{(6.3)}
```

Le premier est l'enveloppe de Jensen à dendrogramme partagé. Le second
moyenne l'observation avant le carré et perd précisément la norme
$`L^2(O)`$ de weak recovery.

L'écriture (6.1) n'est pas encore une construction d'opérateur. Les espaces
de ports, les cellules et les corridors varient avec $D$. Il faut définir
un embedding commun, ou employer la Palm cross--cross de la section 4 comme
représentation échantillonnable du carré.

## 7. État récursif minimal

Une récursion locale ne peut pas utiliser les deux corridors nus comme état
de Markov. Deux réalisations ayant les mêmes arbres visibles et des
potentiels extérieurs différents peuvent avoir des Jacobiennes différentes.
Un état suffisant doit conserver au minimum :

- les mêmes endpoints physiques $`i,j`$ ;
- les rangs réels, buckets complets et incidences des deux corridors ;
- l'identification des sommets et arêtes physiques partagés ;
- tous les ports du séparateur courant dans chaque réplique ;
- la loi de bord extérieure exacte sur leurs configurations, modulo
  normalisation et flip global ;
- la loi conditionnelle commune de la partie encore inexplorée, ou une
  statistique démontrée suffisante de l'observation $O$.

L'espace tangent d'une réplique est donc en général un secteur centré
$`\mathcal H_{r,x}\subset L_0^2(\ell_{r,x})`$ de la loi de ports
$`\ell_{r,x}`$. Une Jacobienne descendante est un opérateur

```math
\mathbf J_{r,c}:
\mathcal H_{r,X_c}
\longrightarrow
\mathcal H_{r,x},
\qquad
\mathbf J_{1,c}\otimes\mathbf J_{2,c}
\text{ agit dans le secteur overlap.}
\qquad\text{(7.1)}
```

La réduction à un scalaire doit être démontrée ; elle ne peut pas être
postulée lorsque les buckets sont multiports.

## 8. Opérateur signé et normalisation positive

Supposons qu'un état $x$ représente une famille de paires physiques et que
ses descendants $c$ en forment une partition. Si $`N_c(x)`$ est le nombre de
paires du descendant, la sélection Palm emploie uniquement les poids
positifs

```math
\omega_c(x)
=
\frac{N_c(x)}{\sum_{c'}N_{c'}(x)}.
\qquad\text{(8.1)}
```

L'opérateur répliqué correctement normalisé est

```math
(\mathcal L^{(2)}h)(x)
=
\mathbb E
\left[
\sum_c
\omega_c(x)
\left(
\mathbf J_{1,c}\otimes\mathbf J_{2,c}
\right)
h(X_c)
\ \middle|\
X=x
\right].
\qquad\text{(8.2)}
```

Dans une réduction scalaire, le tenseur est remplacé par le produit signé
$`J_{1,c}J_{2,c}`$. Ce produit ne doit être ni transformé en probabilité,
ni normalisé par sa somme, ni remplacé par sa valeur absolue. Toute la
normalisation probabiliste appartient à $`\omega_c`$ et au noyau positif de
l'état ; toute la cancellation appartient aux Jacobiennes.

Un test immédiat est

```math
\mathbf J_{1,c}=\mathbf J_{2,c}=I
\text{ pour tout }c
\quad\Longrightarrow\quad
\mathcal L^{(2)}\mathbf1=\mathbf1.
\qquad\text{(8.3)}
```

Une somme nue sur les descendants introduirait un faux facteur de
branchement.

L'état doit aussi retenir toute variable commune qui corrèle la Jacobienne
présente avec la réponse future. Le contre-exemple minimal est une variable
cachée $`H\in\{-1,+1\}`$ telle que

```math
A_0=H,
\qquad
h(X_1)=H.
\qquad\text{(8.4)}
```

Après oubli de $H$, on trouverait $`\mathbb E[A_0]=0`$, alors que
$`\mathbb E[A_0h(X_1)]=1`$. L'observation commune et les messages extérieurs
peuvent jouer exactement ce rôle. Un rayon spectral favorable après
compression non suffisante de l'état serait donc artificiel.

## 9. Lemme de contraction candidat

Soit $`H_L`$ le nombre de couches sans source entre les endpoints et la
région LCA sous la Palm cross--cross. Un lemme suffisant et falsifiable
prend la forme suivante sur les volumes où $`\alpha_L^\times>0`$. Si cette
masse tend elle-même vers zéro, $`\mathcal D_L^\times\le
\alpha_L^\times`$ ferme déjà la cible sans opérateur.

1. $`H_L\to\infty`$ en probabilité Palm.
2. Les états augmentés forment un noyau de Markov positif suffisant.
3. L'overlap exact satisfait, avant la source terminale du LCA,

```math
q_k
=
\left(
\mathbf J_{1,k}\otimes\mathbf J_{2,k}
\right)
q_{k+1}
+
r_k.
\qquad\text{(9.1)}
```

   La donnée terminale vérifie uniformément
   $`\|q_{H_L}\|_{\mathcal B}\le C_0`$.

4. Pour un espace de sections $`\mathcal B`$, les produits inhomogènes
   vérifient, avec $`r<1`$,

```math
\left\|
\mathcal L_{L,0}^{(2)}
\cdots
\mathcal L_{L,k-1}^{(2)}
\right\|_{\mathcal B\to\mathcal B}
\le
Cr^k.
\qquad\text{(9.2)}
```

5. Les erreurs de compression ou de linéarisation satisfont

```math
\sum_{k<H_L}
\left\|
\mathcal L_{L,0}^{(2)}
\cdots
\mathcal L_{L,k-1}^{(2)}
\mathbb E[r_k\mid X_k]
\right\|_{\mathcal B}
\longrightarrow0.
\qquad\text{(9.3)}
```

Ces cinq propriétés impliquent

```math
\mathbb E_{\mathbb P_L^\times}[m_1m_2]
\longrightarrow0,
\qquad
\mathcal D_L^\times(p)\longrightarrow0,
\qquad
Q_L(p)\longrightarrow0.
\qquad\text{(9.4)}
```

Sur les lois complètes de ports, le transfert est linéaire et $`r_k=0`$,
mais l'espace d'état n'est pas de dimension bornée. Une compression par
Jacobiennes doit contrôler un reste non linéaire, par exemple

```math
\|\Phi_x(q)-\mathbf J_xq\|
\le
C_x\|q\|^{1+\varepsilon}.
\qquad\text{(9.5)}
```

Un simple énoncé $`\rho(\mathcal L^{(2)})<1`$ ne suffit que si un opérateur
limite homogène est effectivement construit et si les opérateurs finis y
convergent assez uniformément. Sinon, la bonne porte est la borne de
produits (9.2). Il faut aussi vérifier $`\Phi_x(0)=0`$ dans les couches
itérées : une source affine répétée peut maintenir une réponse non nulle
malgré un rayon spectral inférieur à un.

## 10. Conséquences pour le programme

La porte TRI1 se sépare désormais en deux morceaux.

| porte | statut exact |
|---|---|
| TRI1-s | désintégration scalaire, biais de paire et Palm cross--cross | **fermée** par (2.2)--(4.4) |
| TRI1-o | espace commun de ports et opérateur récursif | **ouvert** |
| TRI2 | rayon spectral sous le biais $`L^2(O)`$ | **ouvert** |
| TRI3 | enveloppe non linéaire et fermeture de weak recovery | **ouvert** |

Le prochain diagnostic doit donc estimer $d_{O,ij}$, pas
$`\mathbb E_D[|g(1-s)m_D|]`$ et pas
$`\mathbb E_D[g(1-s)m_D^2]`$. Deux méthodes restent fidèles :

1. moyenner plusieurs dendrogrammes signés pour chaque même
   $`(O,i,j)`$, puis utiliser un estimateur sans biais de leur carré ;
2. tirer directement deux corridors selon la Palm cross--cross (4.3).

Dans les deux cas, les erreurs standards doivent être clusterisées par
observation. Le
[diagnostic $L=4$](../diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)
reste un audit exact du raffinement commun ; son reste parfois négatif est
compatible avec (2.2). L'extension à des volumes supérieurs doit utiliser
une élimination exacte avec largeur exposée ou, à défaut, signaler
explicitement tout biais de mélange.

## 11. Verdict

La double géante n'exige plus une loi jointe mystérieuse au niveau du
scalaire de weak recovery. Après élimination de la diagonale critique, elle
se réduit à

```math
\boxed{
\mathcal D_L^\times(p)
=
\frac1{n_L^2}
\mathbb E_O
\sum_{i,j}
\left(
\mathbb E_{D\mid O}
\left[
g_D(1-s_D)m_D
\right]
\right)^2.
}
\qquad\text{(11.1)}
```

La difficulté restante est maintenant localisée : construire et contracter
le transfert signé **moyenné en $D$ à observation fixée**, dans un espace
qui conserve tous les ports et facteurs postcritiques. C'est une cible plus
précise que l'ancien opérateur joint schématique, mais ce n'est pas encore
une preuve à $`p=0.81`$.
