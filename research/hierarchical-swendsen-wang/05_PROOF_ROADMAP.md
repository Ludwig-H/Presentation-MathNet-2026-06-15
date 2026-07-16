# Feuille de route des preuves

## Bloc A — Stabiliser les identités finies

### A1. Mesure jointe

**But.** Rédiger une preuve autonome de la loi de $(\sigma,D)$, y compris les arêtes satisfaites qui ne fusionnent jamais avant $1$, les égalités d'horloges et les composantes déjà reliées.

**Livrable.** Un théorème fini : marginale en $\sigma$, conditionnelle $D\mid\sigma$, conditionnelle $\sigma\mid D$.

### A2. Balance des heat baths

**But.** Montrer état par état que les quatre poids $`q_u^{ab}`$ donnent la conditionnelle exacte. Vérifier que des mises à jour successives, avec répétition des nœuds, préservent $\nu(\cdot\mid D)$.

**Livrable.** Un lemme de composition couvrant random scan, parcours déterministe et choix adaptatif dépendant de $D$, avec les restrictions nécessaires.

### A3. Extrémités

**But.** Isoler précisément les hypothèses sous lesquelles :

- racines = Swendsen–Wang ;
- feuilles = heat bath de Glauber ;
- feuilles + acceptation appropriée = Metropolis–Hastings.

**Critère de clôture.** Aucun « égal à » sans noyau et mesure cible explicités.

## Bloc B — Formaliser le critère informationnel

### B1. Deux répliques

Rédiger le théorème
```math
\text{weak recovery à probabilité positive}
\quad\Longleftrightarrow\quad
\liminf Q_n>0
```
dans les conventions exactes du chapitre 11.

Points à traiter :

- algorithmes randomisés ;
- distinction $\liminf/\limsup$ ;
- passage de la corrélation à $`\mathrm{ov}_n`$ ;
- version succès avec probabilité $1-o(1)$ ;
- extension à $K>2$.

### B2. Projection au LCA

Pour chaque paire $i\ne j$, formaliser le noyau qui conserve $D$ et rééchantillonne les deux orientations au nœud $`u_{ij}`$ lorsqu'il existe, ou les deux racines distinctes sinon. Le livrable central est
```math
\langle\sigma_i\sigma_j\rangle_O^2
\le
\mathbb E_{\nu_O}\left[
\mathbf1_{\{i,j\text{ dans le même arbre}\}}\eta_{u_{ij}}
\right],
\qquad
\eta_u=\tanh^2(L_u/2),
```
avec la convention de recoloration des racines distinctes sous a priori uniforme.

La preuve doit contenir explicitement :

1. les quatre événements $00,01,10,11$ ;
2. sur l'événement de connexion, l'identité $`m_u=f_{ij}\,\mathbb E[f_{ij}\mid\mathcal G_u]`$ ;
3. pour la quantité étendue, $`\mathbb E[m_{ij}^{\mathrm{LCA}}\mid\mathcal G_{ij}]=\eta_{ij}^{\mathrm{LCA}}`$, avec valeur nulle pour les racines distinctes ;
4. la différence entre survie depuis la vérité et accord avec une réplique indépendante.

### B3. Somme globale LCA

Après finalisation de A1, mettre sous forme de théorème autonome l'identité déjà dérivée
```math
Q_n
\le
H_n^{\mathrm{LCA}}
=
\frac1{n^2}\mathbb E\left[
n+2\sum_u|C_{u,1}||C_{u,2}|\eta_u
\right]
\le
\frac1{n^2}\mathbb E\sum_R|R|^2.
```
Le second membre doit être identifié exactement à la version second-moment de la borne Swendsen--Wang, et non seulement comparé qualitativement.

### B4. Chaîne LCA pair-spécifique

Formaliser le noyau marginal déjà identifié
```math
K_{ij}^{\mathrm{LCA}}=A^*\mathsf H_{ij}A
```
et sa factorisation comme opérateur auto-adjoint, positif et contractant. Pour
```math
A_{ij}^{(m)}=\langle f_{ij},(K_{ij}^{\mathrm{LCA}})^m f_{ij}\rangle,
```
établir la décroissance spectrale et identifier les conditions exactes sous lesquelles
```math
A_{ij}^{(m)}\downarrow c_{ij}^2.
```
Le verrou asymptotique est un contrôle uniforme permettant de choisir $`m=m_n`$. Le verrou algorithmique distinct est que le noyau dépend de la paire.

### B5. Matrice de persistance

Rédiger le théorème
```math
h_n(S)\to0
\quad\Longrightarrow\quad
\text{impossibilité}.
```
Vérifier explicitement le cas racine
```math
H_S(i,j)=\mathbf1_{\{i\leftrightarrow j\text{ dans }\Pi_1\}}.
```
### B6. Comparaison avec $\theta^{\max}$

Déterminer une version quantitative de $`H_S`$ qui retrouve non seulement la taille de la plus grande composante, mais la borne sur la fraction récupérable. Candidats :
```math
\frac1{n^2}\mathrm{tr}(H_S^2),
\qquad
\frac1n\lambda_{\max}(H_S),
\qquad
\text{profil des valeurs propres de }H_S.
```
## Bloc C — Passer des nœuds à une capacité globale

### C1. Canal exact d'une fusion

Pour chaque $u$, calculer le noyau à quatre états et sa contraction sur l'orientation relative :

- contraction $\chi^2$ ;
- corrélation maximale ;
- coefficient de Dobrushin ;
- second coefficient singulier.

Comparer ces choix sur les petits graphes.

### C2. Ancêtres

L'obstacle central est
```math
q_u^{ab}
\propto
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
```
Trois voies sont à tester :

1. **SDPI conditionnelle directe**, sans factorisation ;
2. **domination** par un canal symétrique indépendant plus informatif ;
3. **élargissement de l'état** du nœud pour rendre le processus markovien sur l'arbre.

La décomposition de départ est additive :
```math
L_u
=
B_u
+\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
```
Sur cactus, $`B_u`$ doit obéir à une récursion de canal binaire symétrique ; sur bandes, à une matrice de transfert ; sur la grille, il faudra contrôler ou dominer sa loi.

### C3. Cycles et multi-terminal

Une fusion voit tous les liens entre $`C_1`$ et $`C_2`$. Employer une SDPI multi-terminale pour éviter la multiplication naïve des contractions des arêtes. Les triangles doivent être le premier cas non trivial.

### C4. Matrice dominante

Construire $`K_D^{\mathrm{info}}`$ telle que
```math
H_S\preceq K_D^{\mathrm{info}}.
```
**Critère de succès minimal.** Pour une recoloration aux racines, retrouver exactement les blocs de composantes. Pour un canal arête par arête, retrouver au moins la borne d'information-percolation.

## Bloc D — Condition suffisante

Une quantité oracle dépendant de $D$ ne suffit pas. Il faut un estimateur fonction de $O=(X,W)$.

Options :

1. belief propagation sur une approximation locale du dendrogramme ;
2. échantillonnage postérieur puis agrégation de plusieurs répliques ;
3. méthode spectrale appliquée à une estimation de $`C_O`$ ;
4. estimateur multiscale sur blocs géométriques, inspiré de la synchronisation sur grilles.

La preuve doit séparer :

- existence informationnelle de l'estimateur ;
- calcul effectif de l'estimateur ;
- temps de mélange de la dynamique proposée.

## Bloc E — Grille triangulaire

### E1. Baselines

Reproduire rigoureusement :
```math
p_c^{\mathrm{edge}}=0.673648\ldots,
\quad
p_c^\triangle=0.719224\ldots,
\quad
p_c^{\mathrm{info}}=0.794659\ldots
```
avec les hypothèses et conditions de bord.

### E2. Coupes déterministes

Calculer les contractions pour $m=1,2,3,\ldots$ liens transverses, conditionnellement à $k$ liens satisfaits et à $\beta$.

Inclure le contrôle exact au temps de fusion
```math
k\mid(m,\beta=t)
\stackrel d=
1+\mathrm{Bin}\left(m-1,
\mathrm{logistic}(u_p(1-t))\right)
```
pour une coupe annealed conditionnée, puis quantifier précisément ce qui reste biaisé par le choix aléatoire de $`E_u`$.

### E3. Biais de Kruskal

Décrire la loi d'une fusion conditionnellement à la filtration passée. Aucun remplacement de $k$ par une variable binomiale sans preuve.

### E4. Cas intermédiaires

Prouver d'abord le critère sur :

- arbres ;
- cactus triangulaires ;
- bandes de largeur fixe.

Utiliser ces cas pour choisir la bonne définition de capacité avant d'attaquer le plan.

Calibrations obligatoires :

- chemin : $`A_1=(2p-1)^\ell`$ contre $c^2=(2p-1)^{2\ell}$ ;
- triangle : calcul exact de $`A_m`$ et de $`\eta_\triangle`$ ;
- cactus de triangles : factorisation de la valeur exacte et convergence de la chaîne LCA.

### E5. Objectif

Obtenir une borne rigoureuse
```math
p_\star>0.794659\ldots.
```
Le point $0.8358058\ldots$ reste un repère conjectural, non une cible à annoncer comme acquise.

## Bloc F — Mélange et interpolation algorithmique

Pour un nombre $m$ de mises à jour, comparer :

- racines seulement ;
- feuilles seulement ;
- random scan uniforme ;
- random scan pondéré par $`|C_u|`$, $`\rho_u`$ ou $`\beta_u`$ ;
- parcours bas-haut et haut-bas.

Mesurer et, si possible, borner
```math
h_n(S_m)-h_{n,\infty}.
```
Cette partie répond à une question algorithmique distincte du seuil de weak recovery : quel parcours fait perdre l'information le plus vite tout en conservant Gibbs ?

## Garde-fous

1. Une composante géante FK n'est pas une condition suffisante d'ordre magnétique en présence de frustration.
2. Une capacité calculée avec $D$ révélé ne donne directement qu'une borne oracle.
3. Une valeur issue de dualité/répliques non rigoureuse reste étiquetée conjecture.
4. Toute affirmation asymptotique doit préciser l'exhaustion, les conditions de bord et l'ordre des limites.
5. Les égalités exactes sont d'abord vérifiées sur un graphe fini par énumération.
6. Une mise à jour LCA pair-spécifique ne définit pas un seul parcours commun à toutes les paires ; ne pas lui attribuer automatiquement une matrice PSD échantillon par échantillon.
7. Une probabilité de survie depuis $\Sigma$ sous le couplage de Nishimori est une autocorrélation Markovienne, pas la précision d'une nouvelle réplique indépendante.

## Ordre conseillé des premiers travaux

1. Finaliser A1 puis écrire B1--B3 comme premier paquet de résultats réutilisables.
2. Vérifier A1--A3 et le théorème LCA sur tous les graphes à au plus quatre sommets.
3. Prouver B4 et calculer $`A_m`$ sur un chemin et un triangle.
4. Calculer C1 et le message $`B_u`$ sur un cactus de triangles.
5. Comparer quantitativement $`H_n^{\mathrm{LCA}}`$ à l'information-percolation.
6. Traiter les bandes triangulaires avant de viser une nouvelle constante sur la grille entière.
