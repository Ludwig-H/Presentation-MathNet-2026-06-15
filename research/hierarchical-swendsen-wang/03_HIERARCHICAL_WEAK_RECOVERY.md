# Vers un critère hiérarchique nécessaire et suffisant

Ce fichier fixe la cible informationnelle exacte, puis décrit ce que la dynamique hiérarchique peut contrôler à temps fini.

## 1. Critère exact à deux répliques

On se place dans le cas binaire symétrique, avec a priori i.i.d. uniforme. Pour
```math
R_n(\sigma,\tau)=\frac1n\sum_{i=1}^n\sigma_i\tau_i,
```
soient $\sigma^{(1)},\sigma^{(2)}$ deux tirages indépendants de $`\mu_O`$ conditionnellement à $O=(X,W)$. Posons
```math
Q_n
=
\mathbb E\left\langle
R_n(\sigma^{(1)},\sigma^{(2)})^2
\right\rangle.
```
Si
```math
C_O=\left\langle\sigma\sigma^\top\right\rangle_{\mu_O},
```
alors
```math
Q_n
=
\frac1{n^2}\,
\mathbb E\mathrm{tr}(C_O^2).
```
### Proposition — statut : immédiat à formaliser

La propriété
```math
\exists\,\varepsilon,\eta>0,\quad
\liminf_{n\to\infty}
\mathbb P\bigl(|R_n(\Sigma_n,\tau_n)|\ge\varepsilon\bigr)
\ge\eta
```
pour un algorithme $`\tau_n`$ est équivalente à
```math
\liminf_{n\to\infty}Q_n>0.
```
Cette formulation correspond à la weak recovery avec avantage de probabilité positif. La version « succès avec probabilité tendant vers $1$ » demande en plus un argument de concentration ou d'amplification.

### Preuve courte

Si $`\tau_n`$ récupère avec paramètres $(\varepsilon,\eta)$, l'identité de Nishimori donne
```math
\varepsilon^2\eta
\le
\mathbb E R_n(\sigma,\tau_n)^2.
```
Conditionnellement à $O$,
```math
\mathbb E\left[R_n(\sigma,\tau_n)^2\mid O,\tau_n\right]
=
\frac{\tau_n^\top C_O\tau_n}{n^2}
\le
\frac{\lambda_{\max}(C_O)}n.
```
Comme
```math
\frac{\lambda_{\max}(C_O)}n
\le
\frac{\sqrt{\mathrm{tr}(C_O^2)}}n,
```
Jensen donne
```math
\varepsilon^2\eta\le\sqrt{Q_n}.
```
Réciproquement, prendre comme algorithme une réplique postérieure $`\tau_n=\sigma^{(1)}`$. Par Nishimori,
```math
\mathbb E R_n(\Sigma_n,\tau_n)^2=Q_n.
```
Puisque $`0\le R_n^2\le1`$, une borne $`Q_n\ge q>0`$ implique par exemple
```math
\mathbb P\left(
|R_n(\Sigma_n,\tau_n)|\ge\sqrt{q/2}
\right)
\ge \frac q2.
```
Ce critère est la cible nécessaire et suffisante. La percolation n'est qu'une manière de majorer $`Q_n`$.

## 2. Obstruction fournie par un parcours hiérarchique

Partons de $`\sigma\sim\mu_O`$, tirons $D\mid\sigma$, puis appliquons à $D$ fixé un programme $S$ de heat baths internes. Soit $\sigma'$ la configuration obtenue. Comme le programme conserve $`\nu_O(\cdot\mid D)`$, $\sigma'$ a encore pour marginale $`\mu_O`$.

Posons
```math
\zeta_i=\sigma_i\sigma_i',
\qquad
H_S(O,\sigma,D)
=
\mathbb E\left[\zeta\zeta^\top\mid O,\sigma,D\right],
```
l'espérance portant sur le hasard des heat baths. La matrice $`H_S`$ est positive semi-définie et sa diagonale vaut $1$.

Définissons
```math
h_n(S)
=
\mathbb E\left[\frac{\lambda_{\max}(H_S)}n\right],
\qquad
h_n^\star=\inf_{S\in\mathcal S_n}h_n(S),
```
où $`\mathcal S_n`$ est une classe explicitement fixée de parcours admissibles.

### Théorème d'obstruction — statut : immédiat à formaliser

Si une weak recovery existe avec paramètres $(\varepsilon,\eta)$, alors, pour tout parcours admissible $S$,
```math
h_n(S)\ge\varepsilon^2\eta+o(1).
```
En particulier,
```math
h_n^\star\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
### Preuve courte

Conditionnellement à $`O,\sigma,D,\tau_n`$, posons $`a_i=\sigma_i\tau_{n,i}`$. Le hasard du parcours est indépendant de celui de l'algorithme et
```math
R_n(\sigma',\tau_n)=\frac1n a^\top\zeta.
```
Donc
```math
\mathbb E\left[R_n(\sigma',\tau_n)^2
\mid O,\sigma,D,\tau_n\right]
=
\frac{a^\top H_Sa}{n^2}
\le
\frac{\lambda_{\max}(H_S)}n,
```
car $`\|a\|_2^2=n`$. La marginale de $\sigma'$ est postérieure ; l'identité de Nishimori donne la même probabilité de succès contre $\sigma'$ que contre la vérité.

## 3. Swendsen–Wang est exactement le cas racine

Pour une recoloration **globale, uniforme et indépendante** des composantes $C$ de $`\Pi_1`$,
```math
H_S(i,j)=
\mathbf 1_{\{i,j\text{ appartiennent au même }C\}}.
```
Après permutation des sommets, $`H_S`$ est bloc-diagonale avec un bloc de $1$ de taille $|C|$ pour chaque composante. Par conséquent,
```math
\lambda_{\max}(H_S)=\max_C|C|.
```
Ainsi,
```math
\frac{\max_C|C|}{n}\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
Le théorème $\theta^{\max}$ du chapitre 11 est plus quantitatif pour la fraction récupérable ; la matrice $`H_S`$ est en revanche la bonne extension à des flips corrélés et multi-échelles.

Cette identité concerne la recoloration globale de chaque composante, et non le heat bath à quatre états d'un nœud interne supérieur.

## 4. Borne pairwise au nœud de coalescence

Pour une paire $i\ne j$, soit $`u_{ij}=\mathrm{LCA}_D(i,j)`$ lorsqu'elle appartient au même arbre. Le heat bath des deux fils de $`u_{ij}`$ est la projection conditionnelle qui efface leurs orientations absolues tout en gardant leurs relations internes. Avec
```math
L_u=\log\frac{q_u^{00}+q_u^{11}}{q_u^{10}+q_u^{01}},
\qquad
\eta_u=\tanh^2\frac{L_u}{2},
```
le lemme de projection donne, sous a priori binaire uniforme,
```math
\eta_{ij}^{\mathrm{LCA}}
:=
\begin{cases}
\eta_{u_{ij}},&i,j\text{ dans le même arbre},\\
0,&i,j\text{ dans deux racines distinctes},
\end{cases}
```
et
```math
\left\langle\sigma_i\sigma_j\right\rangle_O^2
\le
\mathbb E_{\nu_O}[\eta_{ij}^{\mathrm{LCA}}].
```
Sur la diagonale, on pose $`\eta_{ii}^{\mathrm{LCA}}=1`$.

En regroupant les paires par LCA,
```math
\boxed{
Q_n
\le
H_n^{\mathrm{LCA}}
:=
\frac1{n^2}\mathbb E\left[
n+2\sum_{u\in D}|C_{u,1}||C_{u,2}|\eta_u
\right].
}
```
L'identité d'arbre
```math
n+2\sum_u|C_{u,1}||C_{u,2}|
=
\sum_{R\text{ racine}}|R|^2
```
montre que
```math
H_n^{\mathrm{LCA}}
\le
\frac1{n^2}\mathbb E\sum_R|R|^2.
```
Le second moment des composantes Swendsen--Wang, qui implique l'obstruction qualitative du chapitre 11 lorsque les composantes sont sous-macroscopiques, est donc la version non pondérée. Le facteur $`\eta_u`$ retient la fiabilité exacte de la fusion où la paire se rencontre. Cette comparaison ne remplace pas la borne quantitative $`\theta^{\max}`$ sur la fraction récupérable.

### Itération pair-spécifique

Rafraîchissons $D\mid\sigma$, appliquons le heat bath au LCA de la paire, oublions $D$, puis répétons. Le noyau marginal $`K_{ij}^{\mathrm{LCA}}`$ est auto-adjoint et positif. Ainsi
```math
A_{ij}^{(m)}
=
\langle f_{ij},(K_{ij}^{\mathrm{LCA}})^m f_{ij}\rangle_O,
\qquad f_{ij}=\sigma_i\sigma_j,
```
décroît avec $m$, reste supérieur à $`c_{ij}(O)^2`$, et converge vers cette valeur si le noyau est ergodique pour $`f_{ij}`$. On obtient donc une suite de bornes
```math
H_n^{(1)}=H_n^{\mathrm{LCA}}\ge H_n^{(2)}\ge\cdots\ge Q_n,
```
exacte à la limite sous contrôle uniforme du mélange. Cette famille est pair-spécifique : elle ne doit pas être confondue avec une unique matrice de persistance $`H_S`$ issue d'un seul parcours commun à toutes les paires.

Les identités, les quatre événements face à la vérité et la formule exacte contenant $`\beta_u=\xi_{e_u}`$ sont développés dans [06_LCA_SPIN_CORRELATION.md](06_LCA_SPIN_CORRELATION.md).

## 5. Fermeture à l'équilibre

Considérons maintenant la chaîne alternée complète :

1. rafraîchir $D\mid\sigma$ ;
2. mettre à jour $\sigma\mid D$ ;
3. répéter jusqu'au mélange.

Si cette chaîne est ergodique et mélangée, $\sigma'$ devient une réplique postérieure indépendante de $\sigma$, conditionnellement à $O$. Après avoir éliminé $D$,
```math
H_\infty(O,\sigma)
=
\mathrm{diag}(\sigma)\,
C_O\,
\mathrm{diag}(\sigma),
```
donc
```math
h_{n,\infty}
=
\mathbb E\frac{\lambda_{\max}(C_O)}n.
```
Comme $`C_O\succeq0`$ et $`\mathrm{tr}(C_O)=n`$,
```math
Q_n
\le h_{n,\infty}
\le\sqrt{Q_n}.
```
Par le critère à deux répliques,
```math
\liminf h_{n,\infty}>0
\quad\Longleftrightarrow\quad
\text{weak recovery}.
```
Cette identité sépare deux questions :

- le **seuil informationnel**, décrit par $`Q_n`$ ou $`h_{n,\infty}`$ ;
- le **coût dynamique**, décrit par la vitesse à laquelle $`h_n(S_m)`$ approche $`h_{n,\infty}`$ lorsque le nombre de mises à jour croît.

À $D$ fixé, même après mélange, on n'obtient que $`\nu_O(\cdot\mid D)`$. Le rafraîchissement de $D$ est indispensable dans cette section.

## 6. Coefficient informationnel d'un nœud

Pour un nœud $u$, les log-odds de l'orientation relative sont
```math
L_u=\log\frac{q_u^{00}+q_u^{11}}{q_u^{10}+q_u^{01}}.
```
Deux coefficients naturels sont
```math
\rho_u=\left|\tanh\frac{L_u}{2}\right|,
\qquad
\eta_u=\rho_u^2.
```
$`\eta_u`$ est la contraction $\chi^2$ du canal binaire symétrique ayant les mêmes odds. Dans un modèle idéal où les nœuds seraient des canaux indépendants sur un arbre,
```math
K_D(i,j)
=
\prod_{u\in\mathrm{path}_D(i,j)}\rho_u
```
serait la corrélation entre deux feuilles, et la positivité d'une capacité $L^2$ ou d'un rayon spectral associé donnerait un critère de reconstruction de type Kesten–Stigum.

Dans la dynamique réelle, cette factorisation est fausse en général :

- $`q_u^{ab}`$ contient les facteurs des ancêtres ;
- plusieurs fusions dépendent des mêmes arêtes du graphe initial ;
- Kruskal sélectionne les coupes de façon dépendante des horloges ;
- les cycles créent des canaux multi-terminaux.

## 7. Conjecture de capacité hiérarchique

**Conjecture.** Il existe une matrice positive $`K_D^{\mathrm{info}}`$, calculée récursivement à partir de contractions conditionnelles de blocs, telle que
```math
H_S\preceq K_D^{\mathrm{info}}
```
pour un parcours hiérarchique convenable et telle que
```math
\lim_{n\to\infty}
\mathbb E\frac{\lambda_{\max}(K_D^{\mathrm{info}})}n=0
```
soit une condition d'impossibilité. Sous des hypothèses de localité et d'homogénéité, la positivité de la capacité correspondante devrait aussi permettre de construire un estimateur et devenir suffisante.

La première moitié doit être cherchée avec des strong data-processing inequalities conditionnelles, idéalement multi-terminales. La seconde demande un algorithme explicite, par exemple une belief propagation sur des blocs du dendrogramme, et non l'usage oracle de $D$.

## 8. Hiérarchie des objectifs

| Niveau | Énoncé | Statut |
|---|---|---|
| 0 | $`\liminf_n Q_n>0`$ caractérise la weak recovery binaire à probabilité positive | Immédiat à formaliser |
| 1 | $`h_n(S)\to0`$ pour un parcours invariant implique l'impossibilité | Immédiat à formaliser |
| 2 | $`Q_n\le H_n^{\mathrm{LCA}}\le`$ second moment FK | Établi conditionnellement à A1 |
| 3 | $`A_{ij}^{(m)}\downarrow c_{ij}^2`$ à volume fini | Établi sous ergodicité de l'observable |
| 4 | Contrôle uniforme de $`H_n^{(m_n)}-Q_n`$ | À prouver |
| 5 | Comparaison stricte à l'information-percolation | À prouver |
| 6 | Non-disparition d'une capacité calculable $\Rightarrow$ estimateur récupérant | À prouver sous hypothèses |
| 7 | Capacité nulle/positive donne le seuil exact sur la grille triangulaire | Conjecture |

Le niveau intermédiaire désormais concret est : $`H_n^{\mathrm{LCA}}\to0`$ implique l'impossibilité, tandis que $`H_n^{(m)}\downarrow Q_n`$ à volume fini sous ergodicité. La difficulté est de rendre la convergence uniforme en $n$ et calculable.

## 9. Extensions

- $K>2$ : remplacer les corrélations binaires par les représentations non triviales de $`\mathfrak S_K`$ ou par un opérateur de confusion centré.
- A priori non uniforme : centrer par le meilleur random guess et intégrer $`\mu_0`$ dans chaque noyau.
- Interactions triangulaires : employer une SDPI multi-terminale au niveau d'un triangle ou d'une fusion, plutôt qu'un produit de contractions d'arêtes.
- Mélange : comparer random scan, parcours bas-haut, haut-bas et choix adaptatif des nœuds via $`h_n(S)`$.
