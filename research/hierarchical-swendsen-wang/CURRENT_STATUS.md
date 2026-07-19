# Statut scientifique actuel

**Dernière mise à jour : 19 juillet 2026.** Cette page est la source de
vérité du projet. Les anciennes feuilles de route restent consultables dans
[`archive/roadmaps/`](archive/roadmaps/), mais ne fixent plus les priorités.

## 1. Réponse courte

Le dépôt contient déjà une amélioration rigoureuse de la borne de référence
du chapitre 11 :

```math
p_{\mathrm{WR}}
\ge
\frac{809439}{1000000}
=
0.809439.
\qquad\text{(1.1)}
```

Cette borne est obtenue par le
[canal triangulaire multi-état](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).
Elle ne valide pas encore la stratégie hiérarchique.

La voie hiérarchique reste plausible, mais seulement sous une forme
resserrée :

```math
\text{grande déviation géométrique sous une loi ordinaire enracinée dans une paire}
\quad+
\text{budget d'entropie du tilt énergétique}
\quad+
\text{contraction collapsed par blocs}.
\qquad\text{(1.2)}
```

## 2. Ce qui est définitivement acquis

### Voie non hiérarchique

- un canal rationnel à quatre états est less-noisy que le canal physique à
  $`p=0.809439`$ pour tout a priori ;
- les quatre contrôles polynomiaux sont certifiés par Sturm ;
- les hypothèses de sous-criticité du modèle de comparaison sont strictes ;
- information-percolation ferme le passage du canal local à l'absence de
  weak recovery.

La preuve canonique unique est le
[fichier 34](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).
Les certificats à $`0.805`$ et $`0.809`$ sont des jalons désormais subsumés,
classés dans [`archive/certificates/`](archive/certificates/).

### Voie hiérarchique

- la loi jointe observation–réplique–dendrogramme est définie exactement en
  volume fini ;
- chaque update est un heat bath exact de cette loi ;
- la weak recovery se ramène à une corrélation spin–spin pairwise ;
- les taux ancestraux et la loi résiduelle des marques de frontière sont
  calculés sous les conditionnements annoncés ;
- les projections collapsed donnent une identité pythagoricienne exacte de
  dissipation ;
- le corridor complet ne conserve pas plus d'énergie $`L^2`$ que le LCA seul ;
- une chaîne de cactus triangulaires fournit un certificat hiérarchique exact.

Ces résultats sont des briques de preuve. Aucun ne donne encore une borne
hiérarchique nouvelle sur la grille triangulaire.

## 3. Les deux no-go qui ont changé le programme

Le [fichier 29](diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) ferme deux
raccourcis importants.

1. **Criticalisation multiport.** Avancer une fusion réelle vers le niveau
   critique n'est pas une domination de Blackwell uniforme dès que plusieurs
   relations de bord réagissent différemment aux flips descendants.
2. **État local fidèle.** Si l'état de transfert contient assez
   d'information pour rendre le twist mesurable, le déficit Feynman--Kac
   local est exactement nul.

Conséquence : on ne cherche plus une chaîne Markovienne de bord de dimension
fixe ni une contraction uniforme fusion par fusion. On étudie la fonction de
paire effectivement propagée par des projections collapsed.

## 4. L'objet géométrique prioritaire

On tire $`I_L`$ et $`J_L`$ indépendamment et uniformément. Pour ne garder que
des paires lointaines sans perdre une fraction macroscopique de la moyenne,
on peut conditionner par $`d(I_L,J_L)\ge r_L`$, où $`r_L\to\infty`$ et
$`r_L/L\to0`$. On explore symétriquement les deux bras du dendrogramme
jusqu'à leur LCA.

Pour la paire $`\{i,j\}`$, une règle déterministe, symétrique et mesurable
depuis le dendrogramme non marqué sélectionne au plus une arête d'ancrage par
cellule **géométriquement admissible**. Si $`\mathcal A_D(i,j)`$ désigne ces
ancrages et $`\mathcal P_D(i,j)`$ le chemin unique entre les deux feuilles,
le compteur utile est

```math
N_D^{\mathrm{geo}}(i,j)
=
\sum_{e\in\mathcal P_D(i,j)}
\mathbf 1_{\{e\in\mathcal A_D(i,j)\}}.
\qquad\text{(4.1)}
```

Ce compteur est symétrique et mesurable depuis le dendrogramme non marqué. Il
n'est pas appelé distance, car la sélection peut dépendre de la paire. Il
compte seulement les occasions géométriques candidates. Leur activité
énergétique est l'objet d'un lemme séparé ; l'ultramétrique de coalescence
seule localise le LCA, mais ne mesure pas leur accumulation.

Le premier objectif géométrique crédible est une grande déviation de la
forme

```math
\mathbb P_{\mathrm{pair}}
\left(
N_D^{\mathrm{geo}}(i,j)<cK
\right)
\le
e^{-\kappa K},
\qquad
K\asymp\log d(i,j).
\qquad\text{(4.2)}
```

Il faut tenter de la démontrer avec RSW, séparation de bras, circuits et un
vrai lemme de découplage ou de domination conditionnelle. Les labels bruts
sur des annuli disjoints sont indépendants ; les événements d'admissibilité,
après enracinement dans une paire et exploration de Kruskal, ne le sont pas
automatiquement. Un théorème ergodique seul ne donne pas le taux exponentiel
de (4.2).

## 5. Le changement de mesure entropique

Soit
$`W_{\mathrm{in}}=\|M_{\mathrm{in}}\|_{L^2(\pi_D)}^2\in[0,1]`$ l'énergie
de paire à l'entrée du macrobloc. Si
$`a=\mathbb E_{\mathbb P}[W_{\mathrm{in}}]=0`$, la corrélation est déjà
nulle. Supposons donc $`a>0`$ et posons

```math
a=\mathbb E_{\mathbb P}[W_{\mathrm{in}}],
\qquad
\frac{d\mathbb Q}{d\mathbb P}=\frac{W_{\mathrm{in}}}{a}.
\qquad\text{(5.1)}
```

La mesure $`\mathbb Q`$ est exactement la loi inclinée par l'énergie que voit
le quotient de dissipation. Elle n'est appelée Palm énergétique qu'après la
construction d'une mesure aléatoire stationnaire appropriée. Comme
$`\log W_{\mathrm{in}}\le0`$,

```math
D(\mathbb Q\Vert\mathbb P)
\le
\log\frac1a.
\qquad\text{(5.2)}
```

Pour l'événement mauvais
$`A_K=\{N_D^{\mathrm{geo}}(i,j)<cK\}`$, la contraction de l'entropie relative vers
l'indicatrice de $`A_K`$ donne

```math
\mathbb Q(A_K)
\le
\frac{\log(1/a)+\log2}{\kappa K},
\qquad
\text{si }\mathbb P(A_K)\le e^{-\kappa K}.
\qquad\text{(5.3)}
```

En particulier, si $`a>e^{-\varepsilon\kappa K}`$,

```math
\mathbb Q(A_K)
\le
\varepsilon+\frac{\log2}{\kappa K}.
\qquad\text{(5.4)}
```

La preuve visée repose alors sur une dichotomie.

- Si $`a`$ est déjà exponentiellement petit en $`K`$, la corrélation est
  détruite sans autre argument.
- Sinon, le budget d'entropie du tilt est insuffisant pour concentrer toute
  l'énergie sur les corridors géométriquement exceptionnels.

Cette idée remplace la minoration directe, très fragile, de la fréquence des
bons blocs sous une loi inclinée par l'énergie d'entrée. Elle ne permet pas
d'assembler automatiquement des pertes micro-locales dont le tilt changerait
de $`M_k^2`$ à $`M_{k+1}^2`$.

## 6. Le lemme analytique encore manquant

La géométrie et l'entropie ne suffisent que si les $`K`$ échelles produisent
une contraction qui tend vers zéro, mesurée avec le **même tilt d'entrée**.
Posons $`G_K=A_K^c`$ et
$`W_{\mathrm{out}}=\|M_{\mathrm{out}}\|_{L^2(\pi_D)}^2`$. Une cible
minimale est : il existe $`\lambda>0`$ tel que

```math
\mathbb E_{\mathbb P}[W_{\mathrm{out}}]
\le
\mathbb E_{\mathbb P}
\left[W_{\mathrm{in}}\mathbf1_{A_K}\right]
+
e^{-\lambda K}
\mathbb E_{\mathbb P}
\left[W_{\mathrm{in}}\mathbf1_{G_K}\right].
\qquad\text{(6.1)}
```

Après division par $`a`$, (6.1) donne exactement

```math
\frac{\mathbb E[W_{\mathrm{out}}]}a
\le
\mathbb Q(A_K)+e^{-\lambda K}.
\qquad\text{(6.2)}
```

Cette formulation regroupe les updates avant le changement de mesure. Une
marge uniforme sur tous les potentiels extérieurs est fausse ; le lemme doit
porter sur les potentiels atteints ou moyenner leur queue polarisée. Une
autre option serait une suite explicite de macroblocs dont les tilts sont
contrôlés, mais un unique facteur fixe $`1-\eta`$ ne suffit pas.

## 7. Portes go/no-go

| porte | énoncé falsifiable | décision si échec |
|---|---|---|
| G0 — variable auxiliaire | le compteur symétrique et les cellules géométriques sont mesurables depuis le dendrogramme non marqué | redériver explicitement le heat bath marqué avant tout calcul |
| G1 — géométrie | obtenir (4.2) avec $`K\asymp\log d(i,j)`$ | abandonner la simple distance de comptage |
| G2 — bloc | obtenir (6.1), ou une itération équivalente, sous une classe atteignable de potentiels | abandonner cette famille de cellules |
| G3 — tilt | utiliser le même tilt d'entrée dans G1 et G2 | regrouper davantage les updates ou construire un calendrier explicite de tilts |
| G4 — fermeture | transformer la décroissance pairwise en absence de weak recovery | vérifier les moyennes de paire et les composantes distinctes avant toute annonce de seuil |

On ne lance un grand programme numérique qu'après G0. On ne lance une preuve
multiscalaire complète qu'après un témoin robuste pour G2.

## 8. Priorités immédiates

1. formaliser l'exploration symétrique enracinée dans une paire ;
2. fixer une définition mesurable des cellules critiques blindées ;
3. prouver une version pilote de G1 sur des annuli dyadiques espacés ;
4. dériver la formule de variance locale pour un bloc, puis chercher G2 ;
5. assembler le lemme entropique seulement après les deux validations
   précédentes ;
6. convertir le résultat en seuil numérique si et seulement si les constantes
   restent strictes.

La construction d'une mesure invariante complète « vue depuis l'ancêtre »
n'est pas prioritaire : le passage à l'ancêtre est plusieurs-vers-un et le
tilt énergétique varie avec le niveau. Cette construction risque de recréer
les mêmes difficultés sous un formalisme plus lourd.

## 9. Références de navigation

- [programme actif détaillé](active/35_DISTANCE_ENTROPIE_ERGODICITE.md) ;
- [socle de dissipation quadratique](active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) ;
- [cellules critiques à deux updates](active/33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) ;
- [audit multiport et état local](diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) ;
- [bibliographie commentée](references/LITERATURE.md) ;
- [catalogue de toutes les notes](INDEX.md).
