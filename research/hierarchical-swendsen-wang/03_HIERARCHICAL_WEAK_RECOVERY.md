# Vers un critère hiérarchique nécessaire et suffisant

Ce fichier fixe la cible informationnelle exacte, puis décrit ce que la dynamique hiérarchique peut contrôler à temps fini.

## 1. Critère exact à deux répliques

On se place dans le cas binaire symétrique, avec a priori i.i.d. uniforme. Pour


$$
R_n(\sigma,\tau)=\frac1n\sum_{i=1}^n\sigma_i\tau_i,
$$


soient $\sigma^{(1)},\sigma^{(2)}$ deux tirages indépendants de $`\mu_O`$ conditionnellement à $O=(X,W)$. Posons


$$
Q_n
=
\mathbb E\left\langle
R_n(\sigma^{(1)},\sigma^{(2)})^2
\right\rangle.
$$


Si


$$
C_O=\left\langle\sigma\sigma^\top\right\rangle_{\mu_O},
$$


alors


$$
Q_n
=
\frac1{n^2}\,
\mathbb E\operatorname{tr}(C_O^2).
$$


### Proposition — statut : immédiat à formaliser

La propriété


$$
\exists\,\varepsilon,\eta>0,\quad
\liminf_{n\to\infty}
\mathbb P\bigl(|R_n(\Sigma_n,\tau_n)|\ge\varepsilon\bigr)
\ge\eta
$$


pour un algorithme $`\tau_n`$ est équivalente à


$$
\liminf_{n\to\infty}Q_n>0.
$$


Cette formulation correspond à la weak recovery avec avantage de probabilité positif. La version « succès avec probabilité tendant vers $1$ » demande en plus un argument de concentration ou d'amplification.

### Preuve courte

Si $`\tau_n`$ récupère avec paramètres $(\varepsilon,\eta)$, l'identité de Nishimori donne


$$
\varepsilon^2\eta
\le
\mathbb E R_n(\sigma,\tau_n)^2.
$$


Conditionnellement à $O$,


$$
\mathbb E\left[R_n(\sigma,\tau_n)^2\mid O,\tau_n\right]
=
\frac{\tau_n^\top C_O\tau_n}{n^2}
\le
\frac{\lambda_{\max}(C_O)}n.
$$


Comme


$$
\frac{\lambda_{\max}(C_O)}n
\le
\frac{\sqrt{\operatorname{tr}(C_O^2)}}n,
$$


Jensen donne


$$
\varepsilon^2\eta\le\sqrt{Q_n}.
$$


Réciproquement, prendre comme algorithme une réplique postérieure $`\tau_n=\sigma^{(1)}`$. Par Nishimori,


$$
\mathbb E R_n(\Sigma_n,\tau_n)^2=Q_n.
$$


Puisque $`0\le R_n^2\le1`$, une borne $`Q_n\ge q>0`$ implique par exemple


$$
\mathbb P\left(
|R_n(\Sigma_n,\tau_n)|\ge\sqrt{q/2}
\right)
\ge \frac q2.
$$


Ce critère est la cible nécessaire et suffisante. La percolation n'est qu'une manière de majorer $`Q_n`$.

## 2. Obstruction fournie par un parcours hiérarchique

Partons de $`\sigma\sim\mu_O`$, tirons $D\mid\sigma$, puis appliquons à $D$ fixé un programme $S$ de heat baths internes. Soit $\sigma'$ la configuration obtenue. Comme le programme conserve $`\nu_O(\cdot\mid D)`$, $\sigma'$ a encore pour marginale $`\mu_O`$.

Posons


$$
\zeta_i=\sigma_i\sigma_i',
\qquad
H_S(O,\sigma,D)
=
\mathbb E\left[\zeta\zeta^\top\mid O,\sigma,D\right],
$$


l'espérance portant sur le hasard des heat baths. La matrice $`H_S`$ est positive semi-définie et sa diagonale vaut $1$.

Définissons


$$
h_n(S)
=
\mathbb E\left[\frac{\lambda_{\max}(H_S)}n\right],
\qquad
h_n^\star=\inf_{S\in\mathcal S_n}h_n(S),
$$


où $`\mathcal S_n`$ est une classe explicitement fixée de parcours admissibles.

### Théorème d'obstruction — statut : immédiat à formaliser

Si une weak recovery existe avec paramètres $(\varepsilon,\eta)$, alors, pour tout parcours admissible $S$,


$$
h_n(S)\ge\varepsilon^2\eta+o(1).
$$


En particulier,


$$
h_n^\star\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
$$


### Preuve courte

Conditionnellement à $`O,\sigma,D,\tau_n`$, posons $`a_i=\sigma_i\tau_{n,i}`$. Le hasard du parcours est indépendant de celui de l'algorithme et


$$
R_n(\sigma',\tau_n)=\frac1n a^\top\zeta.
$$


Donc


$$
\mathbb E\left[R_n(\sigma',\tau_n)^2
\mid O,\sigma,D,\tau_n\right]
=
\frac{a^\top H_Sa}{n^2}
\le
\frac{\lambda_{\max}(H_S)}n,
$$


car $`\|a\|_2^2=n`$. La marginale de $\sigma'$ est postérieure ; l'identité de Nishimori donne la même probabilité de succès contre $\sigma'$ que contre la vérité.

## 3. Swendsen–Wang est exactement le cas racine

Pour une recoloration **globale, uniforme et indépendante** des composantes $C$ de $`\Pi_1`$,


$$
H_S(i,j)=
\mathbf 1_{\{i,j\text{ appartiennent au même }C\}}.
$$


Après permutation des sommets, $`H_S`$ est bloc-diagonale avec un bloc de $1$ de taille $|C|$ pour chaque composante. Par conséquent,


$$
\lambda_{\max}(H_S)=\max_C|C|.
$$


Ainsi,


$$
\frac{\max_C|C|}{n}\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
$$


Le théorème $\theta^{\max}$ du chapitre 11 est plus quantitatif pour la fraction récupérable ; la matrice $`H_S`$ est en revanche la bonne extension à des flips corrélés et multi-échelles.

Cette identité concerne la recoloration globale de chaque composante, et non le heat bath à quatre états d'un nœud interne supérieur.

## 4. Fermeture à l'équilibre

Considérons maintenant la chaîne alternée complète :

1. rafraîchir $D\mid\sigma$ ;
2. mettre à jour $\sigma\mid D$ ;
3. répéter jusqu'au mélange.

Si cette chaîne est ergodique et mélangée, $\sigma'$ devient une réplique postérieure indépendante de $\sigma$, conditionnellement à $O$. Après avoir éliminé $D$,


$$
H_\infty(O,\sigma)
=
\operatorname{diag}(\sigma)\,
C_O\,
\operatorname{diag}(\sigma),
$$


donc


$$
h_{n,\infty}
=
\mathbb E\frac{\lambda_{\max}(C_O)}n.
$$


Comme $`C_O\succeq0`$ et $`\operatorname{tr}(C_O)=n`$,


$$
Q_n
\le h_{n,\infty}
\le\sqrt{Q_n}.
$$


Par le critère à deux répliques,


$$
\liminf h_{n,\infty}>0
\quad\Longleftrightarrow\quad
\text{weak recovery}.
$$


Cette identité sépare deux questions :

- le **seuil informationnel**, décrit par $`Q_n`$ ou $`h_{n,\infty}`$ ;
- le **coût dynamique**, décrit par la vitesse à laquelle $`h_n(S_m)`$ approche $`h_{n,\infty}`$ lorsque le nombre de mises à jour croît.

À $D$ fixé, même après mélange, on n'obtient que $`\nu_O(\cdot\mid D)`$. Le rafraîchissement de $D$ est indispensable dans cette section.

## 5. Coefficient informationnel d'un nœud

Pour un nœud $u$, les log-odds de l'orientation relative sont


$$
L_u=\log\frac{q_u^{00}+q_u^{11}}{q_u^{10}+q_u^{01}}.
$$


Deux coefficients naturels sont


$$
\rho_u=\left|\tanh\frac{L_u}{2}\right|,
\qquad
\eta_u=\rho_u^2.
$$


$`\eta_u`$ est la contraction $\chi^2$ du canal binaire symétrique ayant les mêmes odds. Dans un modèle idéal où les nœuds seraient des canaux indépendants sur un arbre,


$$
K_D(i,j)
=
\prod_{u\in\operatorname{path}_D(i,j)}\rho_u
$$


serait la corrélation entre deux feuilles, et la positivité d'une capacité $L^2$ ou d'un rayon spectral associé donnerait un critère de reconstruction de type Kesten–Stigum.

Dans la dynamique réelle, cette factorisation est fausse en général :

- $`q_u^{ab}`$ contient les facteurs des ancêtres ;
- plusieurs fusions dépendent des mêmes arêtes du graphe initial ;
- Kruskal sélectionne les coupes de façon dépendante des horloges ;
- les cycles créent des canaux multi-terminaux.

## 6. Conjecture de capacité hiérarchique

**Conjecture.** Il existe une matrice positive $`K_D^{\mathrm{info}}`$, calculée récursivement à partir de contractions conditionnelles de blocs, telle que


$$
H_S\preceq K_D^{\mathrm{info}}
$$


pour un parcours hiérarchique convenable et telle que


$$
\lim_{n\to\infty}
\mathbb E\frac{\lambda_{\max}(K_D^{\mathrm{info}})}n=0
$$


soit une condition d'impossibilité. Sous des hypothèses de localité et d'homogénéité, la positivité de la capacité correspondante devrait aussi permettre de construire un estimateur et devenir suffisante.

La première moitié doit être cherchée avec des strong data-processing inequalities conditionnelles, idéalement multi-terminales. La seconde demande un algorithme explicite, par exemple une belief propagation sur des blocs du dendrogramme, et non l'usage oracle de $D$.

## 7. Hiérarchie des objectifs

| Niveau | Énoncé | Statut |
|---|---|---|
| 0 | $`Q_n>0`$ caractérise la weak recovery binaire à probabilité positive | Immédiat à formaliser |
| 1 | $`h_n(S)\to0`$ pour un parcours invariant implique l'impossibilité | Immédiat à formaliser |
| 2 | Borne calculable de $`H_S`$ par contractions de nœuds/blocs | À prouver |
| 3 | La borne domine l'information-percolation arête par arête | À prouver |
| 4 | Positivité de la capacité $\Rightarrow$ estimateur récupérant | À prouver sous hypothèses |
| 5 | Capacité nulle/positive donne le seuil exact sur la grille triangulaire | Conjecture |

## 8. Extensions

- $K>2$ : remplacer les corrélations binaires par les représentations non triviales de $`\mathfrak S_K`$ ou par un opérateur de confusion centré.
- A priori non uniforme : centrer par le meilleur random guess et intégrer $`\mu_0`$ dans chaque noyau.
- Interactions triangulaires : employer une SDPI multi-terminale au niveau d'un triangle ou d'une fusion, plutôt qu'un produit de contractions d'arêtes.
- Mélange : comparer random scan, parcours bas-haut, haut-bas et choix adaptatif des nœuds via $`h_n(S)`$.
