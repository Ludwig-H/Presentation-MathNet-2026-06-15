# Feuille de route des preuves

## Bloc A — Stabiliser les identités finies

### A1. Mesure jointe

**But.** Rédiger une preuve autonome de la loi de \((\sigma,D)\), y compris les arêtes satisfaites qui ne fusionnent jamais avant \(1\), les égalités d'horloges et les composantes déjà reliées.

**Livrable.** Un théorème fini : marginale en \(\sigma\), conditionnelle \(D\mid\sigma\), conditionnelle \(\sigma\mid D\).

### A2. Balance des heat baths

**But.** Montrer état par état que les quatre poids \(q_u^{ab}\) donnent la conditionnelle exacte. Vérifier que des mises à jour successives, avec répétition des nœuds, préservent \(\nu(\cdot\mid D)\).

**Livrable.** Un lemme de composition couvrant random scan, parcours déterministe et choix adaptatif dépendant de \(D\), avec les restrictions nécessaires.

### A3. Extrémités

**But.** Isoler précisément les hypothèses sous lesquelles :

- racines = Swendsen–Wang ;
- feuilles = heat bath de Glauber ;
- feuilles + acceptation appropriée = Metropolis–Hastings.

**Critère de clôture.** Aucun « égal à » sans noyau et mesure cible explicités.

## Bloc B — Formaliser le critère informationnel

### B1. Deux répliques

Rédiger le théorème

\[
\text{weak recovery à probabilité positive}
\quad\Longleftrightarrow\quad
\liminf Q_n>0
\]

dans les conventions exactes du chapitre 11.

Points à traiter :

- algorithmes randomisés ;
- distinction \(\liminf/\limsup\) ;
- passage de la corrélation à \(\operatorname{ov}_n\) ;
- version succès avec probabilité \(1-o(1)\) ;
- extension à \(K>2\).

### B2. Matrice de persistance

Rédiger le théorème

\[
h_n(S)\to0
\quad\Longrightarrow\quad
\text{impossibilité}.
\]

Vérifier explicitement le cas racine

\[
H_S(i,j)=\mathbf1_{\{i\leftrightarrow j\text{ dans }\Pi_1\}}.
\]

### B3. Comparaison avec \(\theta^{\max}\)

Déterminer une version quantitative de \(H_S\) qui retrouve non seulement la taille de la plus grande composante, mais la borne sur la fraction récupérable. Candidats :

\[
\frac1{n^2}\operatorname{tr}(H_S^2),
\qquad
\frac1n\lambda_{\max}(H_S),
\qquad
\text{profil des valeurs propres de }H_S.
\]

## Bloc C — Passer des nœuds à une capacité globale

### C1. Canal exact d'une fusion

Pour chaque \(u\), calculer le noyau à quatre états et sa contraction sur l'orientation relative :

- contraction \(\chi^2\) ;
- corrélation maximale ;
- coefficient de Dobrushin ;
- second coefficient singulier.

Comparer ces choix sur les petits graphes.

### C2. Ancêtres

L'obstacle central est

\[
q_u^{ab}
\propto
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
\]

Trois voies sont à tester :

1. **SDPI conditionnelle directe**, sans factorisation ;
2. **domination** par un canal symétrique indépendant plus informatif ;
3. **élargissement de l'état** du nœud pour rendre le processus markovien sur l'arbre.

### C3. Cycles et multi-terminal

Une fusion voit tous les liens entre \(C_1\) et \(C_2\). Employer une SDPI multi-terminale pour éviter la multiplication naïve des contractions des arêtes. Les triangles doivent être le premier cas non trivial.

### C4. Matrice dominante

Construire \(K_D^{\mathrm{info}}\) telle que

\[
H_S\preceq K_D^{\mathrm{info}}.
\]

**Critère de succès minimal.** Pour une recoloration aux racines, retrouver exactement les blocs de composantes. Pour un canal arête par arête, retrouver au moins la borne d'information-percolation.

## Bloc D — Condition suffisante

Une quantité oracle dépendant de \(D\) ne suffit pas. Il faut un estimateur fonction de \(O=(X,W)\).

Options :

1. belief propagation sur une approximation locale du dendrogramme ;
2. échantillonnage postérieur puis agrégation de plusieurs répliques ;
3. méthode spectrale appliquée à une estimation de \(C_O\) ;
4. estimateur multiscale sur blocs géométriques, inspiré de la synchronisation sur grilles.

La preuve doit séparer :

- existence informationnelle de l'estimateur ;
- calcul effectif de l'estimateur ;
- temps de mélange de la dynamique proposée.

## Bloc E — Grille triangulaire

### E1. Baselines

Reproduire rigoureusement :

\[
p_c^{\mathrm{edge}}=0.673648\ldots,
\quad
p_c^\triangle=0.719224\ldots,
\quad
p_c^{\mathrm{info}}=0.794659\ldots
\]

avec les hypothèses et conditions de bord.

### E2. Coupes déterministes

Calculer les contractions pour \(m=1,2,3,\ldots\) liens transverses, conditionnellement à \(k\) liens satisfaits et à \(\beta\).

### E3. Biais de Kruskal

Décrire la loi d'une fusion conditionnellement à la filtration passée. Aucun remplacement de \(k\) par une variable binomiale sans preuve.

### E4. Cas intermédiaires

Prouver d'abord le critère sur :

- arbres ;
- cactus triangulaires ;
- bandes de largeur fixe.

Utiliser ces cas pour choisir la bonne définition de capacité avant d'attaquer le plan.

### E5. Objectif

Obtenir une borne rigoureuse

\[
p_\star>0.794659\ldots.
\]

Le point \(0.8358058\ldots\) reste un repère conjectural, non une cible à annoncer comme acquise.

## Bloc F — Mélange et interpolation algorithmique

Pour un nombre \(m\) de mises à jour, comparer :

- racines seulement ;
- feuilles seulement ;
- random scan uniforme ;
- random scan pondéré par \(|C_u|\), \(\rho_u\) ou \(\beta_u\) ;
- parcours bas-haut et haut-bas.

Mesurer et, si possible, borner

\[
h_n(S_m)-h_{n,\infty}.
\]

Cette partie répond à une question algorithmique distincte du seuil de weak recovery : quel parcours fait perdre l'information le plus vite tout en conservant Gibbs ?

## Garde-fous

1. Une composante géante FK n'est pas une condition suffisante d'ordre magnétique en présence de frustration.
2. Une capacité calculée avec \(D\) révélé ne donne directement qu'une borne oracle.
3. Une valeur issue de dualité/répliques non rigoureuse reste étiquetée conjecture.
4. Toute affirmation asymptotique doit préciser l'exhaustion, les conditions de bord et l'ordre des limites.
5. Les égalités exactes sont d'abord vérifiées sur un graphe fini par énumération.

## Ordre conseillé des premiers travaux

1. Écrire proprement B1 et B2 : ce sont les deux premiers résultats réutilisables.
2. Vérifier A1–A3 sur tous les graphes à au plus quatre sommets.
3. Calculer C1 sur une coupe déterministe.
4. Montrer que C4 retrouve l'information-percolation.
5. Traiter un arbre puis un cactus de triangles.
6. Seulement ensuite viser une nouvelle constante sur la grille triangulaire.
