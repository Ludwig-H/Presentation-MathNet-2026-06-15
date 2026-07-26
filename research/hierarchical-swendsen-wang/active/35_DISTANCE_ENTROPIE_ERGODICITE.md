# Moteur conditionnel : distance, entropie et ergodicité

**Statut : actif mais subordonné au contrôle du reste signé sur la double
géante ; le premier test spectral à un dendrogramme est défavorable à
$`L=4`$ ; aucun nouveau seuil revendiqué.**

> [!IMPORTANT]
> La [cible prioritaire](38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) emploie deux
> Gibbs exacts sur deux dendrogrammes entiers indépendants. Avant d'adapter le
> présent moteur, il fallait d'abord mesurer l'enveloppe spectrale à un
> dendrogramme. Le
> [diagnostic à p = 0,81](../diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)
> la trouve macroscopique à $`L=4`$ ; la cible expérimentale devient donc le
> reste signé inter-cellules de la double géante. Un opérateur à deux
> dendrogrammes ne devient une cible légitime qu'après avoir spécifié sa loi
> jointe sans réinjecter l'overlap à démontrer.

**Rôle :** transformer l'intuition « beaucoup d'occasions critiques de
mélange le long des deux bras » en une suite de lemmes falsifiables.

**Prérequis :** [cadre exact](../foundations/01_MATHEMATICAL_FRAMEWORK.md),
[critère pairwise](../foundations/03_HIERARCHICAL_WEAK_RECOVERY.md),
[projections collapsed](../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md) et
[audit des no-go](../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md).

**Ne démontre pas :** une abondance de cellules, une contraction uniforme,
ni une nouvelle borne de weak recovery.

## 1. Décision scientifique

L'ergodicité de la grille triangulaire est utile, mais elle ne remplace pas
à elle seule le contrôle sous la mesure inclinée par l'énergie. La combinaison
candidate est

```math
\boxed{
\text{grande déviation géométrique sous une loi de paire ordinaire}
\quad+
\text{contrôle KL du changement de mesure énergétique}
\quad+
\text{contraction collapsed par blocs}.
}
\qquad\text{(1.1)}
```

Cette route est plus directe que la construction complète d'une mesure
stationnaire vue depuis l'ancêtre. Elle respecte aussi le diagnostic du petit
tore : la dissipation existe, mais elle est portée par une famille rare de
cellules near-critical.

## 2. Quel champ est stationnaire ?

Sous la loi annealed de Nishimori, soit $`\sigma`$ une réplique postérieure
et, pour une arête $`e=\{x,y\}`$,

```math
R_e=O_e\sigma_x\sigma_y.
\qquad\text{(2.1)}
```

Le champ $`(R_e)_e`$ a la loi du bruit planté **sous la loi annealed**. Dans
le modèle homogène, il est i.i.d. En ajoutant des marques d'horloge i.i.d. et
en construisant le dendrogramme de manière équivariante, on obtient un champ
marqué stationnaire ergodique ; le dendrogramme non marqué en est un facteur
équivariant.

Ce fait autorise des théorèmes ergodiques pour des observables locales
translatées. Il ne donne pas automatiquement une loi stationnaire de
« l'environnement vu depuis le prochain ancêtre » : le passage d'une fusion
à son parent est plusieurs-vers-un. Il ne donne pas non plus d'ergodicité
quenched à observation fixée, de contrôle uniforme pour des déplacements
$`z_L\to\infty`$, ni de stationnarité du champ global d'énergie.

## 3. Trois distances à ne pas confondre

### 3.1 Ultramétrique de coalescence

La quantité la plus canonique est

```math
\tau_D(x,y)
=
\inf\{\beta:x\text{ et }y\text{ appartiennent à la même composante de }\Pi_\beta\}.
\qquad\text{(3.1)}
```

Elle est symétrique, mesurable depuis le dendrogramme non marqué et localise
exactement le niveau du LCA. Elle ne compte toutefois aucune accumulation :
pour des points lointains, elle se concentre près du niveau critique au lieu
de diverger.

### 3.2 Compteur symétrique enraciné dans la paire

Pour une paire non ordonnée $`\{i,j\}`$, on explore les deux bras avec une
règle déterministe qui ne distingue pas $`i`$ de $`j`$. La règle sélectionne
au plus une arête d'ancrage par cellule candidate et ne lit que le dendrogramme
non marqué. Si $`\mathcal A_D(i,j)\subseteq\mathcal P_D(i,j)`$ est l'ensemble
des ancrages retenus sur le chemin unique entre les deux feuilles, on pose

```math
N_D^{\mathrm{geo}}(i,j)
=
\sum_{e\in\mathcal P_D(i,j)}
\mathbf 1_{\{e\in\mathcal A_D(i,j)\}}.
\qquad\text{(3.2)}
```

Le compteur mesure les occasions candidates sur les deux bras ; leur
contraction énergétique est un théorème séparé. Il est symétrique si
$`\mathcal A_D(i,j)=\mathcal A_D(j,i)`$, mais ce n'est en général ni une
distance ni une pseudo-distance, car les ancrages peuvent dépendre de la
paire. Si une famille globale d'arêtes admissibles, indépendante de la paire,
suffit ultérieurement, la même somme définit alors une pseudo-distance
d'arbre. Prouver une croissance de $`N_D^{\mathrm{geo}}`$ avec une queue assez
forte est le véritable lemme géométrique.

### 3.3 Distance de premier passage critique

Une métrique $`0/1`$ construite sur les indicatrices d'ouverture **calibrées
au rang critique** est une bonne source d'intuition multiscalaire, mais une
mauvaise candidate à une loi des grands nombres linéaire. Les arêtes
simplement vérifiées ont ici une densité voisine de $`0.81`$ et sont
supercritiques : elles ne définissent pas cette FPP critique. Dans la FPP
critique de sites sur la grille triangulaire, le temps de passage est d'ordre
logarithmique, pas linéaire ; voir
[Yao](https://arxiv.org/abs/1310.1247).

Ce résultat ne se transfère pas directement à la percolation d'arêtes ni au
dendrogramme du GSBM. Il indique seulement que $`K\asymp\log d(i,j)`$ est une
échelle plausible et que Kingman, normalisé par la distance euclidienne,
verrait une constante nulle.

## 4. Définition minimale d'une bonne cellule

Une cellule candidate relie deux projections collapsed consécutives sur un
bras. Elle est **géométriquement admissible** si elle satisfait les trois
premières propriétés ci-dessous. Elle est **analytiquement bonne** si le
quatrième point est ensuite démontré.

1. Son rang supérieur appartient à une fenêtre near-critical définie par des
   probabilités de traversée, sans importer d'exposant non disponible.
2. Des événements de séparation de bras isolent deux routes concurrentes
   entre les interfaces intérieure et extérieure d'un annulus.
3. L'événement d'admissibilité et la règle d'ancrage sont mesurables depuis le
   seul dendrogramme non marqué ; ni observation supplémentaire, ni spin, ni
   identité d'arête gagnante ne sont ajoutés silencieusement.
4. Le bloc possède une perte énergétique positive pour le potentiel entrant
   effectivement atteint.

Le point 4 n'entre donc pas dans $`\ell_{\mathrm{geo}}`$ : sinon le coût
deviendrait tautologique et ne serait plus mesurable depuis la seule
géométrie. On définit d'abord une cellule géométriquement admissible, puis on
prouve une contraction sous la loi des potentiels atteints.

## 5. Lemme géométrique cible

Tirons $`I_L,J_L`$ indépendamment et uniformément, éventuellement sous le
conditionnement $`d(I_L,J_L)\ge r_L`$ avec $`r_L\to\infty`$ et
$`r_L/L\to0`$, puis explorons symétriquement leurs deux boules de Kruskal
jusqu'à leur rencontre. Soit $`N_K=N_D^{\mathrm{geo}}(I_L,J_L)`$ le nombre
de cellules géométriquement admissibles sur $`K`$ échelles espacées.

La cible est : il existe $`c,\kappa>0`$ tels que

```math
\mathbb P_{\mathrm{pair}}
\left(N_K<cK\right)
\le
e^{-\kappa K}.
\qquad\text{(5.1)}
```

Une preuve plausible emploie :

- des annuli dyadiques séparés par des zones tampons ;
- RSW pour obtenir des probabilités uniformes de circuits et traversées ;
- une séparation de bras pour empêcher qu'une même interface ne porte toutes
  les cellules ;
- une extraction donnant une domination stochastique **minorante** de
  $`N_K`$, ou une majoration directe de l'événement mauvais ;
- une concentration après avoir prouvé ce découplage.

Seuls les labels bruts sur des annuli disjoints sont indépendants. Après
enracinement dans la paire et exploration de Kruskal, l'admissibilité est un
événement global : son indépendance ou sa domination conditionnelle est un
lemme à établir, pas une propriété disponible gratuitement.

Birkhoff ne suffit pas : il donnerait $`N_K/K\to c`$ presque sûrement sous
une loi stationnaire appropriée, mais pas nécessairement le coût
exponentiel des corridors pauvres exigé par le changement de mesure.

## 6. Pourquoi l'entropie intervient naturellement

Soit
$`W_{\mathrm{in}}=\|M_{\mathrm{in}}\|_{L^2(\pi_D)}^2\in[0,1]`$ l'énergie
de la corrélation de paire à l'entrée d'un groupe de blocs. Si
$`a=\mathbb E_{\mathbb P}[W_{\mathrm{in}}]=0`$, il n'y a plus rien à
prouver. Supposons $`a>0`$ et posons

```math
a=\mathbb E_{\mathbb P}[W_{\mathrm{in}}],
\qquad
\frac{d\mathbb Q}{d\mathbb P}=\frac{W_{\mathrm{in}}}{a},
\qquad
0\le W_{\mathrm{in}}\le1.
\qquad\text{(6.1)}
```

La loi $`\mathbb Q`$ est d'abord une **loi inclinée par l'énergie**. Elle ne
devient une Palm énergétique qu'après construction d'une mesure aléatoire
stationnaire sur des ancres, à déplacement et niveau fixés. Son coût
entropique est borné par

```math
D(\mathbb Q\Vert\mathbb P)
=
\mathbb E_{\mathbb Q}
\left[
\log\frac{W_{\mathrm{in}}}{a}
\right]
\le
\log\frac1a.
\qquad\text{(6.2)}
```

Soit $`A_K=\{N_K<cK\}`$. La contraction de l'entropie relative par
l'application $`\omega\mapsto\mathbf1_{A_K}(\omega)`$ donne

```math
D(\mathbb Q\Vert\mathbb P)
\ge
\mathbb Q(A_K)
\log\frac1{\mathbb P(A_K)}
-
\log2.
\qquad\text{(6.3)}
```

En combinant (5.1), (6.2) et (6.3),

```math
\mathbb Q(A_K)
\le
\frac{\log(1/a)+\log2}{\kappa K}.
\qquad\text{(6.4)}
```

Cette inégalité produit la bonne dichotomie.

- Si $`a\le e^{-\varepsilon\kappa K}`$, la corrélation est déjà
  exponentiellement petite.
- Sinon, le tilt énergétique ne dispose pas d'assez d'entropie pour placer
  presque toute sa masse sur les rares corridors pauvres.

Plus précisément, si $`a>e^{-\varepsilon\kappa K}`$,

```math
\mathbb Q(A_K)
\le
\varepsilon+\frac{\log2}{\kappa K}.
\qquad\text{(6.5)}
```

La dissipation reste donc une comparaison $`L^2`$. On lui ajoute un contrôle
KL, de type $`x\log x`$, pour transporter une grande déviation depuis la loi
géométrique ordinaire vers le **même tilt d'énergie d'entrée**. Ce n'est pas
un remplacement de $`x^2`$ par $`x\log x`$.

## 7. Lemme analytique de bloc

Pour une paire fixée, soit $`f_{ij}(\sigma)=\sigma_i\sigma_j`$ et soit
$`M_k=P_kf_{ij}`$ la martingale inverse associée à des projections collapsed
emboîtées. L'identité exacte est

```math
\|M_{k-1}\|_2^2-\|M_k\|_2^2
=
\|M_{k-1}-M_k\|_2^2.
\qquad\text{(7.1)}
```

Le lemme requis ne doit pas affirmer une marge pour tout potentiel de bord.
La cellule exacte $`L=4`$ a déjà réfuté cette uniformité. Le calcul local et
la cible de fermeture doivent être distingués.

### Formulation par potentiels atteints

Pour la loi réelle du potentiel entrant, une première cible diagnostique est

```math
\mathbb E
\left[
\|M_{k-1}-M_k\|_2^2
\,\middle|\,
\text{cellule admissible}
\right]
\ge
\eta\,
\mathbb E
\left[
\|M_{k-1}\|_2^2
\,\middle|\,
\text{cellule admissible}
\right].
\qquad\text{(7.2)}
```

Cette perte locale ne se compose pas automatiquement : après chaque update,
la loi inclinée deviendrait proportionnelle à $`M_k^2`$.

### Formulation canonique par macrobloc

Posons $`G_K=\{N_K\ge cK\}`$,
$`W_{\mathrm{in}}=\|M_{\mathrm{in}}\|_2^2`$ et
$`W_{\mathrm{out}}=\|M_{\mathrm{out}}\|_2^2`$. Regrouper les $`K`$ échelles
avant d'incliner, puis montrer qu'il existe $`\lambda>0`$ tel que

```math
\mathbb E_{\mathbb P}[W_{\mathrm{out}}]
\le
\mathbb E_{\mathbb P}
\left[W_{\mathrm{in}}\mathbf1_{G_K^c}\right]
+
e^{-\lambda K}
\mathbb E_{\mathbb P}
\left[W_{\mathrm{in}}\mathbf1_{G_K}\right].
\qquad\text{(7.3)}
```

Avec le tilt unique de (6.1), (7.3) devient

```math
\frac{\mathbb E[W_{\mathrm{out}}]}a
\le
\mathbb Q(G_K^c)+e^{-\lambda K}.
\qquad\text{(7.4)}
```

Cette formulation évite de changer de mesure après chaque micro-update et
s'accorde exactement avec le budget entropique de la section 6. Une autre
possibilité serait une itération démontrée sur un nombre divergent de
macroblocs avec un calendrier explicite des tilts. Un unique facteur fixe
$`1-\eta`$ ne suffit pas.

## 8. Ce que l'ergodicité peut et ne peut pas faire

| étape | apport légitime | insuffisance |
|---|---|---|
| choisir une paire typique | théorèmes ergodiques spatiaux et transport de masse, à déplacement fixé | ne donne pas de contrôle uniforme lorsque la distance de la paire diverge |
| compter des cellules | fréquence asymptotique si le processus indexé est stationnaire | ne donne pas (5.1) sans concentration |
| suivre un bras ancestral | possible après construction d'un environnement vu depuis la cellule | le point-shift parent est plusieurs-vers-un |
| gérer le biais du LCA | l'enracinement dans la paire absorbe naturellement le biais par les descendants | ne supprime pas le conditionnement à deux points |

Les point-shifts non bijectifs peuvent produire une loi invariante singulière
par rapport à la Palm ordinaire ; voir
[Baccelli–Haji-Mirsadeghi](https://arxiv.org/abs/1312.0287). Les transports
invariants et la mass-stationnarité sont développés par
[Last–Thorisson](https://arxiv.org/abs/0906.2062).

## 9. Ordre de travail

### Phase A — cohérence probabiliste

1. écrire l'espace de probabilité de la paire, des marques et du dendrogramme ;
2. vérifier que la définition des cellules n'enrichit pas la variable
   auxiliaire utilisée par le heat bath ;
3. dériver les lois enracinées pertinentes et identifier lesquelles sont
   réellement des Palm, sans facteur $`m`$ ajouté deux fois.

### Phase B — test géométrique

1. définir trois à cinq événements locaux d'annulus ;
2. tester leur stabilité par changement d'échelle sur petits tores ;
3. sélectionner un événement compatible avec une preuve RSW ;
4. tenter (5.1), d'abord par une domination minorante de $`N_K`$ ou une
   majoration directe de son événement de petite queue.

### Phase C — test analytique

1. dériver symboliquement la variance perdue par une cellule ;
2. échantillonner des potentiels réellement atteints ;
3. chercher la borne de macrobloc (7.3), plutôt qu'un infimum sur tout le
   simplexe ;
4. arrêter si la marge pondérée disparaît quand la taille augmente.

### Phase D — assemblage

Assembler (5.1), (6.5) et (7.3) sous le même tilt d'entrée, puis seulement
traduire les constantes en borne de weak recovery.

## 10. Critères d'arrêt

La piste doit être abandonnée ou reformulée si l'un des faits suivants est
établi :

- aucune famille mesurable de cellules ne possède une abondance d'ordre
  $`\log d(i,j)`$ ;
- toute famille abondante a une perte énergétique pondérée qui tend vers
  zéro trop vite ;
- le biais de paire détruit la grande déviation même avant le tilt
  énergétique ;
- la variable auxiliaire nécessaire à la distance modifie le heat bath et
  réintroduit une parité presque révélée ;
- le passage de la contraction pairwise au recouvrement global exige une
  hypothèse non satisfaite.

À l'inverse, obtenir simultanément (5.1) et (7.3), ou une itération
équivalente avec des constantes strictes, justifierait de lancer
l'optimisation quantitative du seuil.

## 11. Dépendances et lectures suivantes

- [statut canonique et priorités](../CURRENT_STATUS.md) ;
- [cellules critiques à deux projections](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) ;
- [dissipation du secteur impair](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) ;
- [loi des coupes conditionnées](../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) ;
- [bibliographie primaire et limites de transfert](../references/LITERATURE.md).
