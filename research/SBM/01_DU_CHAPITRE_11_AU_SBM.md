# Du chapitre 11 au critère à deux répliques

## 1. Même point de départ bayésien

Le chapitre 11 procède dans l'ordre suivant :

1. définir la loi jointe de la vérité cachée et de l'observation ;
2. calculer la postérieure ;
3. définir les algorithmes, le random guess et le recouvrement ;
4. remplacer la vérité par une réplique postérieure grâce à Nishimori ;
5. construire une transition qui laisse la postérieure invariante ;
6. borner le recouvrement au moyen de la géométrie produite par cette
   transition.

On garde cette architecture bayésienne, puis on l'étend par une quantité à
deux répliques et par les dendrogrammes. Les sections de référence sont
§§11.1.1–11.1.5 et §11.3.1 du
[manuscrit](../../Manuscrit_de_thèse.pdf#page=134).

La correspondance est volontairement presque terme à terme :

| chapitre 11 | rôle dans ce dossier |
|---|---|
| Déf. 33, (11.1)–(11.4) | loi jointe du SBM et postérieure Ising exacte |
| Déf. 34–38, (11.7)–(11.15) | algorithmes, random guess, overlap et weak recovery |
| §11.1.4 | séparation weak / almost exact / exact |
| Th. 8 et Cor. 2–3, (11.16)–(11.18) | Nishimori et transition postérieure invariante |
| §11.2.5, (11.25) | extrémité Swendsen--Wang des horloges |
| Th. 10 et Cor. 4, (11.27)–(11.31) | obstruction par percolation à relever au niveau quadratique |
| Prop. 17 et Cor. 5, (11.32)–(11.36) | postérieure signée et comparaison avec la règle de gel classique |

Soient $X_i\in\{-1,+1\}$ des labels i.i.d. uniformes. Conditionnellement à
$X=(X_i)_{i\le n}$, les arêtes sont indépendantes et

```math
\mathbb P(A_{ij}=1\mid X_iX_j=+1)=\frac an,
\qquad
\mathbb P(A_{ij}=1\mid X_iX_j=-1)=\frac bn,
\qquad a>b>0.
\qquad\text{(1.1)}
```

L'observation est $O=A$ et la vérité à recouvrer est $X$. Avec

```math
h_1=\frac12\log\frac ab,
\qquad
h_0=\frac12\log\frac{1-a/n}{1-b/n},
\qquad\text{(1.2)}
```

la postérieure exacte sous l'a priori produit s'écrit

```math
\mu_A(x)
\propto
\exp\left[
h_0\sum_{i<j}x_ix_j
+
(h_1-h_0)
\sum_{\{i,j\}:\,A_{ij}=1}x_ix_j
\right].
\qquad\text{(1.3)}
```

Le terme dense est un potentiel de magnétisation, car

```math
\sum_{i<j}x_ix_j
=
\frac12
\left[
\left(\sum_i x_i\right)^2-n
\right].
\qquad\text{(1.4)}
```

Dans la bisection plantée exactement équilibrée, (1.4) est constant mais la
contrainte $\sum_i x_i=0$ couple les orientations. Dans le modèle à labels
i.i.d., il n'y a pas cette contrainte mais le potentiel dense reste présent.
Cette distinction disparaît dans la limite locale de broadcast ; elle ne
doit pas être oubliée dans une preuve sur le graphe fini.

## 2. Overlap et trois régimes

Comme dans les équations (11.7)–(11.8) du manuscrit, pour
$\sigma,\tau\in\{-1,+1\}^n$,

```math
\mathrm{ov}_n(\sigma,\tau)
=
\frac12
\left(
1+
\left|
\frac1n\sum_{i=1}^n\sigma_i\tau_i
\right|
\right).
\qquad\text{(2.1)}
```

Le score signé sous-jacent est

```math
R_n(\sigma,\tau)
=
\frac1n\sum_{i=1}^n\sigma_i\tau_i.
\qquad\text{(2.2)}
```

Pour l'a priori uniforme binaire :

| régime | exigence |
|---|---|
| weak recovery | $|R_n(X,\widehat X)|$ reste séparé de $0$ |
| almost exact recovery | $|R_n(X,\widehat X)|\to1$ en probabilité |
| exact recovery | $|R_n(X,\widehat X)|=1$ avec probabilité tendant vers $1$ |

Le chapitre 11 privilégie la weak recovery : elle détecte une information
globale sans exiger que la postérieure se concentre sur une seule orbite
$\{X,-X\}$.

Deux précautions du cadre général restent utiles. Hors de l'a priori
i.i.d. uniforme, il faut conserver le niveau $\mathrm{RG}_n(s)$ du meilleur
random guess : l'équilibre des proportions et l'invariance par permutation
ne suffisent pas seuls à imposer le seuil $1/2$. Par ailleurs, dans la preuve
du Théorème 10, les petits clusters de taille inférieure à $\delta n$ ont
une variance totale majorée par $\delta$ ; l'ordre correct est donc de
prendre d'abord la limite en $n$, puis $\delta\downarrow0$.

## 3. Nishimori : remplacer la vérité

Conditionnellement à $A$, tirons une réplique
$\sigma^{(1)}\sim\mu_A$, indépendante de l'aléa d'un algorithme
$\widehat X(A)$. Le Théorème 8 du manuscrit donne

```math
(A,X,\widehat X)
\overset{\mathrm{loi}}{=}
(A,\sigma^{(1)},\widehat X)
\quad
\text{pour tout score mesurable de ces variables}.
\qquad\text{(3.1)}
```

Ainsi, pour étudier la possibilité de récupérer $X$, on peut étudier une
réplique postérieure. Cette substitution exige la **vraie postérieure**.
Le paramètre $\beta$ utilisé plus loin est donc un temps d'horloge
auxiliaire ; il ne tempère pas $\mu_A$.

Le Corollaire 2 du manuscrit permet aussi de produire la réplique par un
noyau $K_A$ qui laisse $\mu_A$ invariante. C'est le rôle de la dynamique
hiérarchique. L'invariance garantit la bonne loi marginale d'une sortie
partant d'un état déjà stationnaire ; elle ne la rend pas indépendante de
cet état initial. Les deux tirages i.i.d. utilisés ci-dessous exigent un
échantillonnage exact, ou deux chaînes indépendantes effectivement amenées à
l'équilibre. La preuve de ce mélange reste ouverte pour le sweep
hiérarchique.

## 4. La quantité exacte est quadratique

Tirons maintenant, conditionnellement à $A$,

```math
\sigma^{(1)},\sigma^{(2)}
\overset{\mathrm{i.i.d.}}{\sim}
\mu_A.
\qquad\text{(4.1)}
```

Posons

```math
Q_n
=
\mathbb E
\left\langle
R_n\left(\sigma^{(1)},\sigma^{(2)}\right)^2
\right\rangle.
\qquad\text{(4.2)}
```

Si

```math
C_A
=
\left\langle
\sigma\sigma^\top
\right\rangle_{\mu_A},
\qquad\text{(4.3)}
```

alors

```math
Q_n
=
\frac1{n^2}
\mathbb E\,\mathrm{tr}(C_A^2)
=
\frac1{n^2}
\sum_{i,j}
\mathbb E
\left[
\left\langle\sigma_i\sigma_j\right\rangle_{\mu_A}^2
\right].
\qquad\text{(4.4)}
```

L'équation (4.4) explique pourquoi un unique dendrogramme **partagé et
échantillonné** ne représente pas le carré sans biais : le carré postérieur
est un produit de deux espérances conditionnelles indépendantes. On peut
bien sûr calculer ce carré analytiquement sans simuler deux dendrogrammes ;
les deux copies sont nécessaires à sa représentation répliquée.

## 5. Critère informationnel à deux répliques

### Proposition 5.1

Pour la version de la weak recovery avec un avantage positif avec une
probabilité positive,

```math
\liminf_{n\to\infty}Q_n>0
\qquad\Longleftrightarrow\qquad
\text{weak recovery informationnellement possible}.
\qquad\text{(5.1)}
```

Le passage à une probabilité de succès tendant vers $1$ demande ensuite une
concentration ou une amplification. Dans le SBM binaire symétrique, elle est
fournie par les théorèmes de seuil rappelés dans la
[note 03](03_PREUVE_DU_SEUIL_WEAK_RECOVERY.md).

### Preuve, sens direct

Supposons qu'un algorithme $\tau=\tau(A)$ et des constantes
$\varepsilon,\eta>0$ vérifient

```math
\mathbb P\left(
\left|R_n(X,\tau)\right|\ge\varepsilon
\right)
\ge\eta+o(1).
\qquad\text{(5.2)}
```

Par Nishimori,

```math
\varepsilon^2\eta+o(1)
\le
\mathbb E
\left[
R_n(\sigma,\tau)^2
\right].
\qquad\text{(5.3)}
```

Conditionnellement à $(A,\tau)$,

```math
\mathbb E
\left[
R_n(\sigma,\tau)^2
\mid A,\tau
\right]
=
\frac{\tau^\top C_A\tau}{n^2}
\le
\frac{\lambda_{\max}(C_A)}n
\le
\frac{\sqrt{\mathrm{tr}(C_A^2)}}n.
\qquad\text{(5.4)}
```

Après espérance et Jensen, le dernier membre est au plus $\sqrt{Q_n}$.
Donc (5.2) impose $\liminf Q_n>0$.

### Preuve, sens réciproque

L'algorithme théorique tire une réplique postérieure
$\tau=\sigma^{(1)}$. Nishimori et (4.2) donnent

```math
\mathbb E
\left[
R_n(X,\tau)^2
\right]
=
Q_n.
\qquad\text{(5.5)}
```

Sous l'hypothèse $\liminf Q_n>0$, il existe $q>0$ tel que
$Q_n\ge q$ pour tout $n$ assez grand. La variable
$Z=R_n(X,\tau)^2\in[0,1]$ vérifie

```math
\mathbb P(Z\ge q/2)
\ge
\frac{q}{2-q}
\ge
\frac q2.
\qquad\text{(5.6)}
```

On obtient donc un overlap non trivial avec une probabilité uniformément
positive. Comme l'a priori est ici i.i.d. uniforme, la Prop. 13 du
manuscrit donne en outre
$\mathrm{RG}_n(1/2+\varepsilon)\to0$ pour tout $\varepsilon>0$.
Cette sortie bat donc bien le meilleur random guess au sens de la Déf. 38.

## 6. Extension exacte du raisonnement de percolation

Le Théorème 10 du chapitre 11 considère une partition gelée et exploite
l'indépendance des recoloriages de ses composantes. Il remplace ensuite
chaque composante par sa taille.

Le lift hiérarchique garde davantage d'information :

```math
\text{partition gelée}
\quad\leadsto\quad
\text{dendrogramme entier},
```

```math
\text{masse des composantes}
\quad\leadsto\quad
\mathbb E
\left[
\left\langle\sigma_i\sigma_j\right\rangle_A^2
\right],
```

```math
\text{une recoloration}
\quad\leadsto\quad
\text{deux Gibbs exacts indépendants}.
```

On reste donc au plus près du chapitre 11 : la vérité est encore remplacée
par une sortie postérieure invariante. La nouveauté est que la weak recovery
est testée dans le secteur quadratique exact plutôt que par la seule masse
macroscopique d'une partition.

## 7. Ce qu'il faut maintenant augmenter

Pour une réplique, introduisons un dendrogramme $D$ conditionnel au spin et
la mesure jointe

```math
\nu_A(d\sigma,dD)
=
\mu_A(d\sigma)R_A(dD\mid\sigma).
\qquad\text{(7.1)}
```

L'augmentation correcte de (4.1) est

```math
(\sigma^{(1)},D^{(1)}),
(\sigma^{(2)},D^{(2)})
\overset{\mathrm{i.i.d.}}{\sim}
\nu_A
\quad\text{conditionnellement à }A.
\qquad\text{(7.2)}
```

Partager $D$ remplacerait le carré d'une moyenne par la moyenne d'un carré.
La [note suivante](02_DEUX_DENDROGRAMMES_A_BETA_C.md) calcule exactement
l'écart entre ces deux expériences.
