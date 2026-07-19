# Programme prioritaire : corridor hiérarchique aux rangs réels

> [!WARNING]
> **Document archivé.** Ce programme retrace une étape de la recherche ; il
> n'est plus prioritaire. Consulter le
> [statut scientifique actuel](../../CURRENT_STATUS.md).

> [!IMPORTANT]
> L'[audit à froid](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) remplace la
> criticalisation comme réduction principale. Celle-ci est fausse pour une
> fusion multiport, même à squelette et tailles fixés. Le programme actif
> traite le corridor réel à ses rangs réalisés. L'oracle critique reste un
> benchmark et le lemme de Blackwell reste valide pour un bucket mono-bit.
> Le no-go $`|U|=K`$ sur l'état fidèle rend en outre la fermeture locale
> bornée improbable. Le pivot prioritaire est la
> [dissipation quadratique du secteur impair](../../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md),
> désormais restreinte aux [cellules critiques
> consécutives](../../active/33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md).

> [!NOTE]
> Le premier gain quantitatif est maintenant établi par une route parallèle :
> le canal triangulaire rationnel du [fichier
> 34](../../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) donne
> $`p_{\mathrm{WR}}\ge0.809439`$. Cette route n'utilise pas le dendrogramme et ne
> doit pas être comptée comme un succès de la dynamique hiérarchique.

Cette page sépare les lemmes déjà prouvés, les réductions conditionnelles et
les verrous qui doivent encore être fermés pour obtenir une nouvelle borne de
weak recovery. Pour l'ordre actuel des travaux, partir de la
[feuille corrigée](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md), puis consulter la
[sous-feuille historique](27_SUBROADMAP_CORRIDOR_P0805.md) et les
[premiers résultats](../../diagnostics/finite_volume/28_FIRST_CORRIDOR_P0805_RESULTS.md).

## 1. Question cible

On considère le GSBM binaire homogène sur le tore triangulaire $`G_L`$. Les
observations s'écrivent, dans la jauge plantée,

```math
O_{xy}=\Sigma_x\Sigma_y Z_{xy},
\qquad
\mathbb P(Z_{xy}=+1)=p,
\qquad
u_p=\log\frac{p}{1-p}.
\qquad\text{(1.1)}
```

Une arête satisfaite par la réplique de référence reçoit une horloge
$`\mathrm{Exp}(u_p)`$ ; une arête insatisfaite reçoit l'horloge infinie. Les
arêtes sonnées avant $`\beta`$ engendrent une partition $`\Pi_\beta`$ et un
dendrogramme de Kruskal.

La meilleure borne rationnelle désormais établie dans ce dossier est :

```math
\boxed{
p\in\left[\frac12,\frac{809439}{10^6}\right]
\quad\Longrightarrow\quad
\text{absence de weak recovery sur le GSBM triangulaire},
}
\qquad\text{(1.2)}
```

donc $`p_{\mathrm{WR}}\ge0.809439>0.8`$. La preuve combine un certificat
less-noisy exact, information-percolation et la percolation triangulaire
corrélée sous-critique ; elle est indépendante de la hiérarchie. Le programme
ci-dessous cherche désormais à expliquer ou améliorer cette borne par la
dynamique hiérarchique. Voir la [feuille de route
resserrée](26_FEUILLE_DE_ROUTE_PSTAR.md).

## 2. Expérience favorable canonique

Pour deux sommets $`i,j`$, posons

```math
\beta_{ij}
:=
\inf\{\beta:i\leftrightarrow j\text{ dans }\Pi_\beta\}.
\qquad\text{(2.1)}
```

Si $`q_c=2\sin(\pi/18)`$ est le seuil de percolation par arêtes de la grille
triangulaire, le temps critique est défini par

```math
q_p(\beta_c)
:=
p(1-e^{-u_p\beta_c})
=q_c.
\qquad\text{(2.2)}
```

L'expression « $`i,j`$ lointains fusionnent à la percolation » signifie en
volume fini l'événement de fenêtre

```math
\mathcal F_{L,\rho,\varepsilon}
=
\left\{
d_L(i,j)\ge\rho L,
\ \beta_c-\varepsilon
\le\beta_{ij}\le\beta_c
\right\}.
\qquad\text{(2.3)}
```

L'égalité exacte $`\beta_{ij}=\beta_c`$ a probabilité nulle. On travaille
donc soit avec (2.3), soit avec la désintégration Palm du flux de fusions.

> [!NOTE]
> Au seuil bidimensionnel infini, il n'existe pas de composante infinie de
> densité positive. Dans ce dossier, « même composante géante à
> $`\beta_c`$ » désigne une composante critique macroscopique ou traversante
> sur le tore fini, avec le LCA de la paire localisé dans la fenêtre (2.3).

### Pourquoi ce cas reste un benchmark utile

Trois mécanismes expliquent son intérêt, avec une restriction essentielle.

1. Une paire macroscopiquement éloignée ne peut pas fusionner à une distance
   fixe sous le seuil avec une probabilité non négligeable.
2. Pour un bucket mono-bit à taille fixée, avancer un niveau postcritique vers
   $`\beta_c`$ améliore exactement le canal au sens de Blackwell. Cette
   affirmation est fausse pour une fusion multiport du corridor collapsed.
3. Si les deux sommets ne sont même pas connectés à $`\beta=1`$, les
   recolorations indépendantes de leurs racines effacent exactement leur
   parité.

Ainsi, l'oracle critique est un benchmark de connexion précoce et de qualité
locale des marques. Il n'est pas une enveloppe informationnelle du corridor
multiport réel.

## 3. Lemme scalaire et contre-exemple multiport

Il faut distinguer deux énoncés.

### Lemme F-scalaire — produit de buckets mono-bit, établi

Fixons une expérience artificielle dans laquelle chaque bucket $`r`$ dépend
d'un unique bit latent $`X_r`$, toutes ses arêtes se complémentant ensemble.
Les observations de buckets sont conditionnellement indépendantes sachant le
vecteur $`X`$, dont le prior peut être arbitrairement corrélé. Pour chaque
niveau postcritique $`t_r\ge\beta_c`$, posons

```math
t_r^{\mathrm{fav}}=\beta_c,
\qquad\text{(3.1)}
```

sans modifier le squelette. Alors, pour toute cible $`F`$, l'expérience
criticalisée est plus informative au sens de Blackwell et son second moment
postérieur est plus grand :

```math
\mathbb E\left[
\mathbb E(F\mid K^{\mathrm{réel}})^2
\right]
\le
\mathbb E\left[
\mathbb E(F\mid K^{\mathrm{fav}})^2
\right].
\qquad\text{(3.2)}
```

Ce résultat est prouvé par l'ordre de Blackwell à taille fixée, puis par
tensorisation conditionnelle au vecteur complet des parités. Il formalise
uniquement :

```math
\boxed{
\text{dans le surrogate produit mono-bit, le canal critique est le plus
informatif.}
}
\qquad\text{(3.3)}
```

Références : [19](../../foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md) et
[20](../../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md).

### Contre-théorème F-multiport — statut : certifié

Dans le corridor collapsed réel, un bucket ancestral peut porter plusieurs
relations variant séparément sous les flips descendants. Avec deux relations,
une gagnante marginalisée et la cible $`F(x)=x_1x_2`$, le canal tardif peut
avoir un second moment postérieur strictement supérieur au canal critique.
Il n'existe alors aucun noyau de dégradation de Blackwell.

Le certificat exact à $`p=0.805`$ et la cellule T2-Kruskal sont dans le
[fichier 29](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md). Une criticalisation ne
peut être réintroduite qu'après une domination cible-spécifique démontrée
sous la véritable loi de bord.

### Conjecture GF — domination de la géométrie critique, ouverte sur la grille

L'énoncé plus fort remplacerait le corridor réel par un corridor réellement
échantillonné sous la Palm de (2.3). Il doit comparer simultanément les
tailles de coupe, les formes, les ports latéraux et les messages ancestraux.

Cette domination est :

- établie exactement sur une chaîne de cactus triangulaires ;
- fausse si l'on compare seulement les temps en laissant varier les tailles
  arbitrairement ;
- ouverte sur le tore triangulaire bidimensionnel.

Le programme ne suppose jamais GF sans la marquer comme hypothèse. La voie
robuste cherche désormais une contraction annealed directement aux rangs
réalisés du squelette observé.

## 4. Les lemmes structurants

| code | lemme | statut | référence |
|---|---|---|---|
| L0 | $`\beta_{ij}>1`$ implique un effacement exact par racines indépendantes | établi | [19](../../foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md), [22](../../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md) |
| L1 | une paire lointaine fusionnant avant $`\beta_c-\delta`$ a une masse asymptotiquement nulle | établi à $`\delta`$ fixé | [12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md), [14](../../foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md) |
| L2 | les arêtes internes aux enfants ne votent pas dans le $`\Lambda_u`$ courant | établi exactement | [14](../../foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md) |
| L3 | conditionnellement à $`\Pi_\beta`$, les marques de frontière sont i.i.d. résiduelles | établi exactement | [14](../../foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md), [25](../../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L4 | la charge d'une coupe instantanée est gouvernée par $`m h_p(\beta)^2`$ | établi | [25](../../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L5 | une coupe fusionnante possède une correction Palm par l'arête gagnante | établi exactement | [09](../../diagnostics/09_CRITICAL_MERGER_ORACLE.md), [25](../../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L6 | l'intensité LCA pré-saut vaut $`mN_\rho`$ et un événement de fusion déjà réalisé se pondère par $`N_\rho`$ | établi en volume fini | [27](27_SUBROADMAP_CORRIDOR_P0805.md), [28](../../diagnostics/finite_volume/28_FIRST_CORRIDOR_P0805_RESULTS.md) |
| L7 | à taille fixée, le canal critique mono-bit Blackwell-domine le canal tardif | établi | [19](../../foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md) |
| L8 | les dégradations se tensorisent dans le surrogate produit mono-bit, même sous un prior corrélé | établi abstraitement ; faux comme description générale du corridor réel | [20](../../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md), [29](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) |
| L9 | le corridor collapsed est au plus persistant qu'un sweep des mêmes nœuds | établi | [20](../../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md), [22](../../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md) |
| L10 | une chaîne de cactus critique perd exponentiellement la corrélation avec la profondeur | établi exactement | [21](../../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md) |
| L11 | beaucoup de coupes screenées uniformément contractantes impliquent une obstruction pairwise | conditionnel à la géométrie | [23](23_OPTIMAL_WEAK_RECOVERY_OBSTRUCTION.md), [24](../../diagnostics/24_SIMPLE_RESIDUAL_BALANCE_OBSTRUCTION.md) |

## 5. Information portée par une coupe

Conditionnellement à la partition complète au temps $`\beta`$, une arête de
frontière est conforme avec probabilité

```math
s_p(\beta)
=
\frac{pe^{-u_p\beta}}{1-p+pe^{-u_p\beta}},
\qquad
h_p(\beta)
=
2s_p(\beta)-1
=
\tanh\left(\frac{u_p(1-\beta)}2\right).
\qquad\text{(5.1)}
```

### Coupe instantanée

Pour une coupe de taille $`m`$ non conditionnée à fusionner ensuite,

```math
K\mid X=+1\sim\mathrm{Bin}(m,s_p(\beta)),
\qquad
K\mid X=-1\sim\mathrm{Bin}(m,1-s_p(\beta)).
\qquad\text{(5.2)}
```

Le rapport signal sur bruit est

```math
\mathrm{SNR}^{\mathrm{snap}}
=
\frac{m h_p(\beta)^2}{1-h_p(\beta)^2}.
\qquad\text{(5.3)}
```

La variable géométrique prioritaire est donc

```math
\mathcal J_v=m_vh_p(\beta_v)^2.
\qquad\text{(5.4)}
```

### Coupe de fusion

Si cette coupe fusionne au niveau $`\beta`$, l'arête gagnante est conforme et

```math
K\mid X=+1
\sim
1+\mathrm{Bin}(m-1,s_p(\beta)).
\qquad\text{(5.5)}
```

Le log-rapport local devient

```math
\ell_{m,K}
=
\log\frac K{m-K}
+u_p(1-\beta)(2K-m).
\qquad\text{(5.6)}
```

À $`\beta=1`$, sa fiabilité locale vaut exactement $`1/m`$. Une coupe de
taille un est donc parfaite ; elle ne peut jamais servir de bloc
contractant.

### Sélection LCA-Palm

Une coupe $`E(A,B)`$ de taille $`m(A,B)`$ fusionne au taux

```math
m(A,B)u_ps_p(\beta).
\qquad\text{(5.7)}
```

Le LCA d'une paire lointaine ajoute le poids

```math
m(A,B)N_\rho(A,B),
\qquad\text{(5.8)}
```

où $`N_\rho(A,B)`$ compte les paires macroscopiquement éloignées séparées par
$`A\mid B`$. La géométrie LCA favorise donc les grandes interfaces et les
fusions portant beaucoup de paires.

## 6. Dynamique recommandée pour l'obstruction

| dynamique | usage | verdict |
|---|---|---|
| LCA seul | oracle local et borne supérieure de persistance | trop favorable pour exploiter la distance |
| sweep top-down | dynamique séquentielle naturelle | utile comme contre-audit ; feedback ancestral difficile |
| sweep bottom-up | parcours finissant au LCA | au plus persistant que le LCA seul |
| corridor collapsed | heat bath conjoint des deux bras | dynamique prioritaire pour la preuve |

Si $`P_{\mathcal C}`$ est le heat bath collapsed du corridor et $`K`$ un
sweep des mêmes nœuds,

```math
\|Kg\|_2^2
=
\|P_{\mathcal C}g\|_2^2
+
\|K(I-P_{\mathcal C})g\|_2^2.
\qquad\text{(6.1)}
```

Le corridor collapsed minimise donc la persistance $`L^2`$ parmi ces
rééchantillonnages. Il est le meilleur choix pour construire une obstruction
pairwise calculable.

## 7. Le verrou des $`\Lambda_v`$ ancestraux

Au nœud $`u`$ qui sépare les deux branches de $`i,j`$, les quatre poids du
heat bath contiennent tous les ancêtres :

```math
q_u^{ab}
\propto
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
\qquad\text{(7.1)}
```

Pour $`v\succ u`$, sa frontière doit être divisée selon son incidence avec
les deux descendants de $`u`$. Le groupe invariant ne peut pas être supprimé
avant d'appliquer

```math
F_v(x)=xe^{(1-\beta_v)x},
\qquad\text{(7.2)}
```

car $`F_v`$ est non linéaire. Une majorité scalaire sur la coupe ne détermine
donc pas le signe du message ancestral.

La variable minimale à contrôler sous la Palm d'événement aux rangs réels est

```math
\left(
m_v,\beta_v,
m_{v,0},m_{v,1},m_{v,2},
Z_v,B_v
\right)_{v\in\mathcal C_{ij}},
\qquad\text{(7.3)}
```

où $`Z_v`$ encode les ports latéraux et $`B_v`$ le message extérieur.

## 8. Théorème conditionnel à fermer

Soient $`I_L,J_L`$ uniformes et lointains. Le critère maître est

```math
\mathbb E\left[
H_{\mathcal C}(I_L,J_L)^2
\right]
\longrightarrow0.
\qquad\text{(8.1)}
```

Pour la branche hiérarchique, après les contre-audits à $`p=0.805`$, la
chaîne conditionnelle est la suivante.

1. **Décomposition exacte.** Les paires précoces ont masse $`o(1)`$, les
   racines distinctes ont persistance nulle et les autres corridors gardent
   leurs rangs effectivement réalisés.
2. **Filtration collapsed.** Intégrer des ensembles croissants
   d'orientations du corridor et poser
   $`M_k=\mathbb E[f_{ij}\mid\mathcal F_k]`$.
3. **Dissipation exacte.** Utiliser l'identité pythagoricienne
   $`\|M_{k-1}\|_2^2-\|M_k\|_2^2=\|M_{k-1}-M_k\|_2^2`$.
4. **Cellule cible-spécifique.** Certifier une perte au second update sur la
   fonction réellement produite, sans norme uniforme sur les constantes.
5. **Lemme géométrique annealed.** Sous la loi marquée du corridor final,
   minorer la perte en proportion de l'énergie entrante sur un nombre
   croissant de cellules consécutives dans des fenêtres near-critical. Les
   annuli génériques ne sont plus une cible.
6. **Clôture.** La composition donne

```math
\mathbb E[H_{\mathcal C}^2]
\le
\mathbb E\left[
\exp\!\left(-\sum_{k=1}^{K_L}\alpha_{k,L}\right)
\right]+o(1)
\longrightarrow0.
\qquad\text{(8.2)}
```

Le théorème pairwise des fichiers [03](../../foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) et
[20](../../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md) transforme alors (8.1) en absence de
weak recovery.

## 9. Ordre de travail

### Priorité 1 — borne rationnelle, fermée

À $`p=809439/10^6`$, la PSD uniforme du canal
$`(a,s,e)=(166642280,55571811,166642287)/(5\times10^8)`$ est certifiée
exactement, avec la marge $`Q_E-Q_Y\ge\mathrm{Var}/(5\times10^7)`$. Le
relèvement aux facteurs, la
tensorisation et le passage au tore sont fermés dans le [fichier
34](../../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md). Le point tangent
$`0.809909\ldots`$ reste une optimisation ouverte.

### Priorité 2 — D1 et D2, effectués

D1 exhibe une seconde perte réelle après propagation sur un witness, mais la
marge uniforme s'annule sur les potentiels de bord. D2, pour des paires à
distance maximale du tore $`L=4`$, montre une dissipation dominée par un petit
nombre de paquets et une queue rare.

### Priorité 3 — population effectuée, pivot critique

L'audit non sélectionné de 302 cellules confirme la queue rare, mais les 14
cellules dans $`|q-q_c|\le0.02`$ portent $`34.1\%`$ de la perte pour
$`4.13\%`$ de l'énergie entrante. La seule suite justifiée est la
[sous-feuille critique](../../active/33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) :
prouver une inégalité locale pondérée sur un compact de potentiels intérieurs,
puis une occupation énergétique répétée de ces cellules.

### Priorité 4 — outlets critiques seulement après le lemme local

Ne construire des outlets épaissis que si le lemme local critique donne une
marge énergétique robuste. Les annuli collapsed à rang arbitraire sont
abandonnés. Si la mesure énergétiquement inclinée se polarise encore sous un
blindage critique, arrêter la branche hiérarchique.

### Priorité 5 — dernier contre-test T2

Conserver la jauge locale comme test d'une éventuelle compression spéciale
des potentiels atteignables. Si elle n'est pas Markov-fermée, ne pas relancer
R2--R4 du fichier 29.

## 10. Calibration à $`p=0.8`$

Au seuil critique,

```math
\beta_c=0.410716539196\ldots,
\qquad
s_c=0.693582222752\ldots,
\qquad
h_c=0.387164445505\ldots.
\qquad\text{(10.1)}
```

Ainsi,

```math
\mathcal J_{m,\beta_c}
=
0.149896307863\ldots\,m.
\qquad\text{(10.2)}
```

Une grande coupe au LCA critique est donc très informative. L'obstruction ne
peut pas reposer sur ce seul bucket. Elle doit exploiter la profondeur du
corridor, des coupes plus tardives ou des blocs ambigus répétés.

Pour un bucket $`m=2`$ neutre,

```math
\Gamma_2(\beta_c;0.8)
=
s_c
=
0.693582222752\ldots.
\qquad\text{(10.3)}
```

Dans le corridor factorisé, $`N`$ tels blocs donnent $`s_c^N`$. Le verrou est
donc leur abondance sous Palm, leur screening et la validité de la
composition sur la grille.

## 11. Contre-audits obligatoires

| raccourci tentant | verdict correct |
|---|---|
| « même composante à $`\beta_c`$ » suffit à prouver la weak recovery | faux : il faut contrôler la loi d'événement réelle et son transfert multiport |
| le LCA critique est globalement le plus favorable sur la grille | faux en général même à squelette, taille et incidences fixés ; vrai dans le surrogate mono-bit et sur le cactus d'articulations |
| les fausses arêtes internes aux clusters votent au nœud courant | faux : elles s'annulent dans le flip relatif |
| une grande coupe critique perd son information | faux à $`p=0.8`$ : sa charge croît comme $`0.149896m`$ |
| une coupe tardive est automatiquement peu informative | faux sans connaître sa taille |
| un bucket $`m=1`$ contracte | faux : l'arête gagnante révèle parfaitement la parité |
| le LCA seul exploite la distance entre $`i,j`$ | faux : il ignore la profondeur des deux bras |
| un message local nul annule la corrélation globale | faux sans screening des ancêtres et routes latérales |
| le certificat cactus prouve le résultat sur la grille | faux : les cycles chevauchants et l'état de bord restent ouverts |
| un nœud de Kruskal réalisé doit encore être pondéré par $`mN_\rho`$ | faux : il faut seulement $`N_\rho`$ ; sinon on crée $`m^2N_\rho`$ |
| le benchmark snapshot à $`q_c`$ Blackwell-domine le corridor final | faux : il change le squelette |
| la marge E1+ inférieure à $`0.3`$ est uniforme en potentiel extérieur | faux : le second moment brut tend vers un sous un champ polarisant |
| conserver davantage de micro-état aide à produire un déficit local | faux : si l'état révèle le twist, $`\lvert U\rvert=K`$ et $`d=0`$ |

## 12. Parcours de lecture recommandé

1. [README](../../README.md) — intuition, statuts et navigation.
2. [01 — cadre mathématique](../../foundations/01_MATHEMATICAL_FRAMEWORK.md) — mesure jointe et
   heat baths.
3. [03 — weak recovery](../../foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) — critère pairwise.
4. [25 — information des coupes](../../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md)
   — géométrie conditionnée et loi Palm.
5. [08 — chaîne ancestrale](../../foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md) — calcul des
   $`\Lambda_v^{ab}`$.
6. [20 — corridor collapsed](../../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md) — projection
   exacte et tensorisation scalaire restreinte.
7. [21 — certificat cactus](../../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md) — premier
   modèle exact.
8. [05 — feuille de route](05_PROOF_ROADMAP.md) — dépendances techniques.
9. [27 — sous-feuille P0805](27_SUBROADMAP_CORRIDOR_P0805.md) — ordre
   falsifiable et portes go/no-go.
10. [28 — premiers résultats](../../diagnostics/finite_volume/28_FIRST_CORRIDOR_P0805_RESULTS.md) — audits
    Palm, cellule E1+ et prochaine cellule T2-Kruskal.
11. [29 — audit à froid](../../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) —
    contre-exemple multiport, T2 réel et feuille de route corrigée.
12. [30 — dissipation du secteur impair](../../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) —
    pivot opératoriel, critère annealed et diagnostic exact.

Les calculs et commandes de validation sont documentés dans
[computations/README.md](../../computations/README.md).
