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

### B5. Réduction à la bande critique

Pour une coupe déterministe $\beta$, formaliser

```math
S_n(\beta)
=
\frac1{n^2}\mathbb E\sum_{C\in\Pi_\beta}|C|^2,
```

```math
\mathcal M_n((\beta,1])
=
\frac2{n^2}
\mathbb E\sum_{u:\,\beta<\beta_u\le1}
|C_{u,1}||C_{u,2}|\eta_u,
```

et démontrer comme théorème autonome

```math
Q_n\le S_n(\beta)+\mathcal M_n((\beta,1]).
```

Le résultat doit inclure :

1. la distinction entre bande pure, naissance de connexion et sprinkling sur le quotient ;
2. la réduction nécessaire et suffisante au score signé de bande lorsque $`S_n(\beta)\to0`$ ;
3. la factorisation connexion $\times$ fiabilité $\times$ cohérence ;
4. les bons quantificateurs pour une paire macroscopique.

### B5 bis. Oracle de fusion critique — calcul local fermé

Le [fichier 09](09_CRITICAL_MERGER_ORACLE.md) établit, pour un bucket critique
homogène et $`B_u=0`$,

```math
\Gamma_m^c(p_{\mathrm{SW}})=\frac1m,
\qquad
\Gamma_m^c(p)\longrightarrow1
\quad(p>p_{\mathrm{SW}}),
```

avec une borne exponentielle et la limite exacte dans la fenêtre
$`p-p_{\mathrm{SW}}\asymp m^{-1/2}`$. Il établit aussi le contre-audit

```math
\mathcal C_{n,\delta}^c
\le
S_n(\beta_c+\delta)-S_n(\beta_c).
```

Ce bloc local est clos. Il ne doit pas être promu en seuil global : le travail
restant consiste à contrôler simultanément la masse des fusions, le message
ancestral, la contraction après marginalisation de $D$ et la cohérence signée.

### B5 ter. Probabilité paire critique — du bucket à la paire lointaine

Le [fichier 15](15_CRITICAL_GIANT_PAIR_FLIP.md) établit pour le bucket local

```math
\overline P_m^c(p)
=
\frac{1+\Gamma_m^c(p)}2,
```

ainsi que

```math
\overline P_m^c(p_{\mathrm{SW}})
=
\frac12+\frac1{2m},
```

et, pour $`p>p_{\mathrm{SW}}`$ fixé,

```math
1-\overline P_m^c(p)
\sim
\frac{C_{m\bmod2}(p)}{\sqrt m}e^{-mI_c(p)},
\qquad
I_c(p)=-\frac12\log(1-h_c(p)^2).
```

Les deux constantes $`C_0(p),C_1(p)`$ sont des séries absolument
convergentes données explicitement dans le fichier 15 ; distinguer la parité
est nécessaire pour avoir un véritable équivalent.

Deux lemmes distincts sont nécessaires pour transporter cette calibration à
une paire lointaine de la grille.

1. **CUT.** La taille $`M_L=|E_{u_{I_LJ_L}}|`$ de la coupe critique tend vers
   l'infini sous le biais de Palm retenu, ou sa loi limite est calculée
   explicitement si elle reste tendue.
2. **ANC.** Le message ancestral vérifie $`B_L/M_L\to0`$, ou plus
   généralement possède une limite jointe avec $`M_L`$ permettant de décider
   le signe de $`a_ch_c+B_L/M_L`$.

Sous CUT et ANC, la probabilité paire complète tend vers $`1`$. Ces hypothèses
sont **à prouver**. La distance des sommets et leur appartenance à la plus
grande composante ne suffisent pas : une interface critique peut être réduite
à une arête pivotale. Un exposant pour la moyenne exige en plus une grande
déviation jointe de $`(K_L,B_L,M_L)`$.

### B5 quater. Probabilités de flip et chemin descendant

Le [fichier 16](16_FLIP_PROBABILITIES_DESCENDANT_PATH.md) ferme les
probabilités élémentaires : racine à deux états, feuille à deux états et
nœud interne à quatre états. Il donne aussi, à tout niveau $`t`$,

```math
K\stackrel d=1+\mathrm{Bin}(m-1,s_p(t)),
```

```math
\ell_{m,k}(t;p)
=
\log\frac{k}{m-k}
+u_p(1-t)(2k-m).
```

Pour un balayage de clusters, la relation d'une paire est exactement le
produit des signes de flip sur les deux bras vers son LCA. Le programme
descendant se sépare alors en trois étapes.

1. **PATH-JOINT.** Estimer la corrélation jointe des décisions de heat bath,
   sans remplacer son espérance par un produit de marginales.
2. **SIDE-MSG.** Marginaliser les branches latérales par des messages de
   frontière, exactement sur cactus puis sur bandes.
3. **PATH-COMP.** Comparer le vrai canal de chemin à l'oracle factorisé
   PATH-FAC. Aucun ordre stochastique n'est actuellement établi.

Une récursion tordue sur un état de frontière calcule déjà exactement cette
corrélation sur cactus et séparateurs bornés. Dans l'oracle factorisé, la
condition $`m_{\min}\ge(c_c(p)^{-1}+\varepsilon)\log H`$ suffit à
préserver la relation le long d'un chemin critique de longueur $`H`$ ; elle
n'est pas encore transférée à la grille complète.

Le chemin de la MSF avec gagnants marqués est exclu de cette comparaison : il
révèle exactement la relation de la réplique génératrice et définit une
variable auxiliaire différente.

### B5 quinquies. Seuil de décorrélation du chemin

Le [fichier 17](17_PATH_DECORRELATION_THRESHOLD.md) remplace la question
qualitative par l'atténuation exacte

```math
A_L(p)=-\sum_{w\in\mathcal P_L}\log\Gamma_{m_w}(t_w;p).
```

Les sous-objectifs géométriques sont désormais ordonnés.

1. Montrer ou réfuter que, pour un $`M`$ fixé,
   $`N_{L,M}=\#\{w:2\le m_w\le M\}`$ diverge sous la loi de paire critique.
   Une réponse positive force déjà PATH-FAC vers $`1/2`$ pour tout
   $`p<1`$ fixé.
2. Si les petites coupes sont rares, mesurer d'abord la fonction de partition
   $`\Phi_L(I)=\sum_w m_w^{-1/2}e^{-Im_w}`$. Son abscisse de transition est le
   seuil géométrique général. Dans le sous-cas
   $`m_L\sim\alpha\log H_L`$, l'oracle régulier possède le seuil explicite

   ```math
   p_{\mathrm{path}}(\alpha)
   =
   \frac{1+q_\triangle+(1-q_\triangle)\sqrt{1-e^{-2/\alpha}}}{2}.
   ```

   Pour le vrai chemin descendant, remplacer cette somme par
   $`\sum_w m_w^{-1/2}e^{-m_wI(t_w;p)}`$ et mesurer le nombre de niveaux dans
   la fenêtre $`\beta_c-t_w=O(m_w^{-1})`$ ; hors de cette fenêtre, leur poids
   relatif est exponentiellement plus petit.

3. Sur cactus puis bandes, calculer les normes $`\kappa_r`$ des opérateurs
   tordus et prouver une contraction de bloc sommable. C'est la version jointe
   rigoureuse du critère $`A_L\to\infty`$.

La constante de Nishimori correspond formellement à
$`\alpha=7.053596192884\ldots`$. Cette égalité est seulement une cible à
contre-auditer par la géométrie ; elle ne doit pas être utilisée comme entrée
du raisonnement.

### B6. Matrice de persistance

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

### B7. Comparaison avec $\theta^{\max}$

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

La voie prioritaire est maintenant le lemme HF du
[fichier 12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md). Pour une paire lointaine
du même arbre, il demande de coupler le log-rapport postcritique complet avec
celui de l'oracle où la séparation a lieu en $`\beta_c`$, de sorte que

```math
|L^{\mathrm{post}}|
\le
|L^{\mathrm c}|+o(1)
```

avec probabilité $`1-o(1)`$. Sous ce lemme, l'annulation de la fiabilité de
l'oracle critique interdit la weak recovery sans supposer que les temps LCA
réels se concentrent au seuil.

La décomposition de départ est additive :
```math
L_u
=
B_u
+\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
```
La décomposition en trois groupes du
[fichier 08](08_ANCESTRAL_LAMBDA_CHAIN.md) donne exactement les quatre
$`\Lambda_v(\sigma^{ab})`$ pour une réalisation fixée. Le
[fichier 10](10_ANCESTRAL_LAMBDA_ESTIMATION.md) donne en plus, conditionnellement
au squelette non marqué, la course pondérée exacte, les moments des quatre
taux, leur concentration et la fonctionnelle $`\mathcal R_u`$ qui contrôle la
queue du message ancestral. Le verrou n'est donc plus le canal des marques,
mais la loi du squelette
$`(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succ u}`$ sous le biais du LCA d'une
paire lointaine critique.

Sur cactus, il faut propager cette loi exacte ; sur bandes, construire une
matrice de transfert certifiée ; sur la grille, établir la convergence des
premiers ancêtres, la sommabilité de $`\mathcal R_u`$ et le contrôle des quatre
coins proches de zéro. La borne

```math
|\eta_u-\eta_u^{(K)}|
\le
\min\left(1,\frac{2\mathcal R_u^{(>K)}}{3\sqrt3}\right)
```

transporte ensuite directement ces résultats vers la quantité de weak
recovery.

Lorsque $`\beta_u\simeq\beta_c`$, le fichier 09 fournit en plus le sandwich
stochastique $`1/2\le s_p(\beta_v)<s_p(\beta_c)`$ pour tous les comptes
ancestraux. Ce sandwich doit être transporté séparément dans chacun des quatre
états : complémenter un groupe renverse ses bornes, et aucune monotonie de
$`\eta_u`$ ne peut être supposée.

### C3. Cycles et multi-terminal

Une fusion voit tous les liens entre $`C_1`$ et $`C_2`$. Employer une SDPI multi-terminale pour éviter la multiplication naïve des contractions des arêtes. Le [fichier 11](11_TRIANGLE_BLOCK_SDPI.md) ferme le premier audit :

- la contraction uniforme d'un triangle est $`\eta_\triangle`$ ;
- sa SDPI globale est $`\gamma_2=2q^2/(1+q^2)`$, donc le facteur scalaire est
  moins bon que la baseline par arêtes ;
- le canal d'effacement multi-état satisfait l'inégalité $`\chi^2`$ voulue
  pour tout a priori $\mu$ tel que $`\max_x\mu_x\le1/2`$ ;
- le cas d'un atome dominant est le lemme $`P_\star`$ encore à prouver.

Ce calcul reste un audit auxiliaire. La prochaine étape de la voie
hiérarchique n'est pas de recalculer $`\eta_\triangle`$ ni de fermer d'abord
$`P_\star`$, mais de contrôler la chaîne des $`\Lambda_v`$ sous le biais de la
paire critique et de prouver HF.

### C4. Matrice dominante

Construire $`K_D^{\mathrm{info}}`$ telle que
```math
H_S\preceq K_D^{\mathrm{info}}.
```

**Critère de succès minimal.** Pour une recoloration aux racines, retrouver exactement les blocs de composantes. Pour un canal arête par arête, retrouver au moins la borne d'information-percolation.

### C5. Capacité du quotient critique

Contracter $`\Pi_{\beta_c}`$, puis remplacer chaque bundle entre deux blocs par sa contraction **après marginalisation du dendrogramme**. Les calibrations obligatoires sont

```math
\gamma_1=(2p-1)^2
```

pour une arête et, pour $m$ observations BSC indépendantes,

```math
\gamma_m^{\mathrm{BSC}}
=
\sum_{k=0}^m
\binom mk p^k(1-p)^{m-k}
\tanh^2\left(\frac{u_p}{2}(2k-m)\right).
```

Employer directement $`\eta_u`$ conditionnellement à $D$ est interdit pour une conclusion suffisante : sur une coupe à une arête, cette quantité vaut $1$ conditionnellement à la fusion et surestime la contraction réelle.

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

conditionnellement au squelette non marqué. Le choix de Kruskal biaise encore la géométrie $`(E_u,m,\beta)`$, mais ne biaise plus les marques résiduelles une fois ce squelette fixé.

### E3. Bande critique et flux pivotal

Établir dans les conventions triangulaires

```math
\beta_c=q_p^{-1}(q_c),
\qquad
t_\chi=q_p^{-1}((2p-1)^2),
```

puis vérifier que $`t_\chi>\beta_c`$ redonne exactement $`p>0.794659\ldots`$. Désintégrer la mesure des fusions avec la formule de Russo et pondérer les pivots par la fiabilité conditionnelle.

Employer la loi multinomiale exacte conditionnellement au squelette non marqué et vérifier le crossover $`m h_p(\beta)^2`$. Sur une fusion choisie par Kruskal, il reste à contrôler la loi du squelette groupé $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)`$, pas à recorriger les marques conditionnelles.

Le fichier 14 fixe trois tests à ne pas confondre :

1. majorité des vrais tardifs contre les faux ;
2. majorité conforme dans chacun des deux groupes affectés d'un ancêtre ;
3. préférence quatre états exacte
   $`q_u^{00}+q_u^{11}>q_u^{10}+q_u^{01}`$.

Le premier test a son seuil à $`0.782432\ldots`$ et ne peut donc pas améliorer
la baseline. Le théorème du cône de Walsh montre que le test 2, joint à la
majorité locale, implique le test 3 sous a priori uniforme. La réciproque
n'est pas nécessaire ; les calculs de cactus et de bandes doivent donc
enregistrer à la fois le certificat 2 et le critère exact 3.

### E4. Géométrie biaisée de Kruskal

Décrire la loi du squelette d'une fusion conditionnellement à la filtration passée. Le noyau des marques est maintenant exact après conditionnement ; le verrou est la distribution des tailles de groupes et des temps le long de la chaîne ancestrale.

### E5. Cas intermédiaires

Prouver d'abord le critère sur :

- arbres ;
- cactus triangulaires ;
- bandes de largeur fixe.

Utiliser ces cas pour choisir la bonne définition de capacité avant d'attaquer le plan.

Calibrations obligatoires :

- chemin : $`A_1=(2p-1)^\ell`$ contre $c^2=(2p-1)^{2\ell}$ ;
- triangle : calcul exact de $`A_m`$ et de $`\eta_\triangle`$ ;
- cactus de triangles : factorisation de la valeur exacte et convergence de la chaîne LCA.

### E6. Objectif

Obtenir une borne rigoureuse

```math
p_\star>0.794659\ldots.
```
Le calcul auxiliaire du fichier 11 produit le candidat conditionnel

```math
p_\star^{\mathrm{cond}}=0.8099092892\ldots,
```

racine de son enveloppe affine multi-état. Il reste conditionnel au lemme
$`P_\star`$ et ne doit pas être cité comme résultat établi ni comme sortie de
la dynamique hiérarchique. La cible prioritaire est la valeur obtenue, ou
l'absence d'amélioration constatée, après calcul de
$`\Gamma_{L,\varepsilon}^{\mathrm{fav}}`$ avec tous les ancêtres. Le point
$0.8358058\ldots$ reste un repère conjectural, non une cible à annoncer comme
acquise.

Le fichier 13 fournit désormais une seconde calibration exacte : l'équation
de Nishimori--Ohzeki est une entropie conditionnelle de face égale à un bit et
se représente par quatre horloges exponentielles. Cela ne modifie pas le
statut du seuil. Le fichier 14 remplace la troncature autoduale non canonique
par la décomposition exacte des liens postcritiques. Le calcul prioritaire sur
le cactus de deux triangles doit produire, pour chaque ancêtre, les trois
comptes conformes, les quatre $`\Lambda_v^{ab}`$, puis le log-rapport de
parité complet. Une énumération indépendante doit contre-auditer la loi
conditionnelle de Kruskal et le critère pair/impair. Sans contrôle séparé du
squelette sélectionné et sans HF, ni l'autodualité de face ni le score oracle
ne donnent une borne de weak recovery.

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
5. Fermer B5 puis retrouver information-percolation par $`t_\chi>\beta_c`$.
6. Calculer C5 sur des bundles déterministes, puis sur un cactus.
7. Comparer quantitativement $`H_n^{\mathrm{LCA}}`$ à l'information-percolation.
8. Traiter les bandes triangulaires avant de viser une nouvelle constante sur la grille entière.
9. Prouver HF sur cactus puis sur bandes en couplant les vecteurs quatre états
   de la chaîne ancestrale critique et postcritique.
10. Conserver le lemme $`P_\star`$ comme piste auxiliaire indépendante ; ne
    pas substituer sa constante conditionnelle au résultat hiérarchique.
