# Point de départ : chapitre 11 et borne de percolation

## 1. Définition du recouvrement

Le manuscrit utilise
```math
\mathrm{ov}_n(\sigma,\tau)
=
\max_{\pi\in\mathfrak S_K}
\frac1n\sum_{i=1}^n
\mathbf 1_{\{\tau_i=\pi(\sigma_i)\}}.
```
Pour $K=2$, avec des spins dans $\{-1,+1\}$,
```math
\mathrm{ov}_n(\sigma,\tau)
=
\frac12\left(
1+
\left|\frac1n\sum_i\sigma_i\tau_i\right|
\right).
```
Le niveau général de random guess doit rester $`\mathrm{RG}_n(s)`$. La simplification au seuil $1/K$ exige notamment un a priori i.i.d. uniforme. « Équilibré et invariant par permutation des noms de labels » ne suffit pas : une loi uniforme sur $\{a,-a\}$, avec $a$ équilibrée, est un contre-exemple pour $K=2$.

## 2. Théorème $\theta^{\max}$ corrigé

Soit $`\kappa_n`$ l'objet gelé produit par une transition de clusters et $`\mathcal C(\kappa_n)`$ sa famille de composantes. Posons
```math
\theta_\delta(\kappa_n)
=
\frac1n\sum_{C\in\mathcal C(\kappa_n)}
|C|\,\mathbf 1_{\{|C|\ge\delta n\}},
```
puis
```math
\theta^{\max}
=
\inf\left\{
\theta\in[0,1]:
\lim_{\delta\downarrow0}\limsup_{n\to\infty}
\mathbb P\left[\theta_\delta(\kappa_n)>\theta\right]=0
\right\}.
```
### Hypothèses

- $K$ est fixé.
- L'a priori est i.i.d. uniforme, ou une hypothèse distincte donne explicitement le bon random guess.
- La transition laisse exactement la postérieure invariante.
- Son hasard est indépendant de l'algorithme $`\tau_n`$, conditionnellement à l'observation et à la configuration de départ.
- Conditionnellement à $`\kappa_n`$, les composantes sont recolorées indépendamment et uniformément.

Sous ces hypothèses, pour tout $\eta>0$ et tout algorithme $`\tau_n`$,
```math
\mathbb P\left[
\mathrm{ov}_n(\Sigma_n,\tau_n)
\ge
\frac1K+\frac{K-1}{K}
\bigl(\theta^{\max}+\eta\bigr)
\right]\longrightarrow0.
```
Ainsi, une récupération au niveau
```math
\frac1K+\frac{K-1}{K}\varepsilon
```
avec probabilité tendant vers $1$ impose
```math
\varepsilon\le\theta^{\max}.
```
### Correction de la preuve pour les petits clusters

À $\delta$ fixé,
```math
\sum_{C:\,|C|<\delta n}
\left(\frac{|C|}{n}\right)^2
\le\delta.
```
La fluctuation due aux petits clusters est donc $`O_{\mathbb P}(\sqrt\delta)`$, et non $`o_{\mathbb P}(1)`$ directement à $\delta$ fixé. L'ordre correct est :

1. obtenir la borne conditionnelle ;
2. prendre $`\limsup_{n\to\infty}`$ ;
3. faire ensuite $\delta\downarrow0$.

Une union sur les $K!$ permutations termine la borne uniforme.

## 3. Ce que ce théorème mesure — et ce qu'il perd

Le théorème exploite une propriété très forte du pas aux racines : les orientations des composantes sont indépendantes et uniformes. Il remplace donc toute la postérieure à l'intérieur d'une composante par la seule taille de cette composante.

Cela perd :

- les poids $`|W_e|`$ au-delà de la décision binaire « gelé/non gelé » ;
- les temps de fusion $`\beta_u`$ ;
- les nombreux liens entre les deux fils d'une fusion ;
- les cycles et triangles qui renforcent ou contredisent une orientation relative ;
- la transmission graduelle de l'information le long du dendrogramme.

Avec des heat baths internes, les descendants ne sont plus recolorés indépendamment. On ne peut donc pas remplacer mécaniquement $\theta^{\max}$ par la masse des composantes à une autre coupe.

## 4. Baseline plus forte déjà connue : information-percolation

Avant d'utiliser la nouvelle dynamique, il faut comparer au meilleur résultat général déjà disponible. Pour un canal binaire sur une arête,
```math
\mathbb P(H_e=1\mid\Sigma_i=\Sigma_j)=f_{\mathrm{in}}(r_e),
\qquad
\mathbb P(H_e=1\mid\Sigma_i\ne\Sigma_j)=f_{\mathrm{out}}(r_e),
```
la contraction $\chi^2$ du canal, sous entrée uniforme, vaut
```math
\boxed{
\eta_e
=
\frac{(f_{\mathrm{in}}(r_e)-f_{\mathrm{out}}(r_e))^2}
{(f_{\mathrm{in}}(r_e)+f_{\mathrm{out}}(r_e))
(2-f_{\mathrm{in}}(r_e)-f_{\mathrm{out}}(r_e))}.
}
```
Les bornes d'information-percolation de Polyanskiy–Wu et Abbe–Boix dominent les corrélations par une percolation indépendante de paramètres $`\eta_e`$. Si cette percolation ne porte pas de masse macroscopique — avec le passage en volume fini correctement justifié — la weak recovery est impossible.

Dans le modèle homogène
```math
f_{\mathrm{in}}=p,
\qquad
f_{\mathrm{out}}=1-p,
```
on obtient
```math
\eta=(2p-1)^2.
```
Sur la grille triangulaire, comme
```math
p_c^{\mathrm{bond}}(\mathbb T)=2\sin(\pi/18),
```
la borne sous-critique donne
```math
p<
\frac{1+\sqrt{2\sin(\pi/18)}}2
=0.794659\ldots
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
Cette borne est déjà plus forte que la borne par arêtes du chapitre 11 et que sa dynamique triangulaire. Le premier résultat hiérarchique utile doit donc :

1. retrouver l'information-percolation comme cas simple ou comme borne dominée ;
2. utiliser les fusions multi-arêtes pour la dépasser ;
3. rester une vraie borne sur l'overlap, pas seulement sur la connectivité.

## 5. Premier raffinement hiérarchique obtenu

Le score LCA du [fichier dédié](06_LCA_SPIN_CORRELATION.md) remplace le second moment brut des composantes par
```math
H_n^{\mathrm{LCA}}
=
\frac1{n^2}\mathbb E\left[
n+2\sum_u|C_{u,1}||C_{u,2}|\eta_u
\right],
\qquad
\eta_u=\tanh^2(L_u/2).
```
On a
```math
Q_n
\le H_n^{\mathrm{LCA}}
\le
\frac1{n^2}\mathbb E\sum_{R\text{ racine}}|R|^2.
```
Le dernier membre est la version second-moment de la percolation Swendsen--Wang. Le nouveau critère est donc déjà un raffinement démontré au niveau de l'algèbre finie, conditionnellement à la formalisation complète de la loi jointe A1. Il reste à le comparer quantitativement à l'information-percolation, qui demeure la baseline rigoureuse à battre sur la grille triangulaire.

## 6. A priori non uniforme

Pour $`\mu_0`$ général, deux modifications sont indispensables :

1. conserver la définition $`\mathrm{RG}_n(s)`$ au lieu de remplacer sans preuve le niveau de référence par $1/K$ ;
2. imposer au noyau de recoloriage $`R_\kappa`$ la balance
```math
    e^{-U_0(\sigma)}R_\kappa(\sigma,\sigma')
    =
    e^{-U_0(\sigma')}R_\kappa(\sigma',\sigma).
```
La loi jointe $`\nu_O(\sigma,D)`$ du fichier précédent est conçue pour intégrer proprement cette seconde contrainte.

## 7. Statut

- Théorème $\theta^{\max}$ : **établi après correction des hypothèses et de l'ordre des limites**.
- Passage précis de la percolation infinie à $\theta^{\max}=0$ selon l'exhaustion finie : **à formaliser**.
- Borne information-percolation : **établie dans les références primaires** ; adaptation exacte aux conventions du manuscrit à rédiger.
- Borne LCA $`Q_n\le H_n^{\mathrm{LCA}}`$ : **algèbre finie établie, conditionnelle à la finalisation de A1**.
- Comparaison de $`H_n^{\mathrm{LCA}}`$ à l'information-percolation et seuil strictement meilleur : **à prouver**.
