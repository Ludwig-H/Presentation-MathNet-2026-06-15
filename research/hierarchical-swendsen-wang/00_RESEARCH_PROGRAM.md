# Programme prioritaire : oracle critique et corridor hiérarchique

> [!IMPORTANT]
> La [feuille de route resserrée](26_FEUILLE_DE_ROUTE_PSTAR.md) remplace
> l'oracle « LCA ponctuel critique » comme réduction principale. Elle traite
> le corridor réel, criticalise seulement les canaux tardifs à squelette et
> tailles fixés, et conserve les attaches postcritiques. L'oracle critique de
> cette page reste un test favorable et un catalogue de lemmes, pas une
> hypothèse typique sur toutes les paires.

Cette page sépare les lemmes déjà prouvés, les réductions conditionnelles et
les verrous qui doivent encore être fermés pour obtenir une nouvelle borne de
weak recovery. Pour l'ordre actuel des travaux, partir de la
[sous-feuille de route](27_SUBROADMAP_CORRIDOR_P0805.md), puis des
[premiers résultats à 0,805](28_FIRST_CORRIDOR_P0805_RESULTS.md).

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

Le premier pré-certificat quantitatif est :

```math
\boxed{
p_0=\frac45
\quad\Longrightarrow\quad
\text{absence de weak recovery sur le GSBM triangulaire.}
}
\qquad\text{(1.2)}
```

Par dégradation BSC, une preuve à $`p_0`$ donnerait aussi l'impossibilité pour
tout $`p\le p_0`$. Ce résultat améliorerait la borne de référence
$`0.794659\ldots`$, mais il n'est pas encore démontré.

Une preuve au seul point $`p_0=4/5`$ ne donnerait toutefois pas encore un
seuil strictement supérieur à $`0.8`$. Le premier objectif publiable renforcé
est de prolonger le certificat, par continuité avec marge ou directement, au
point rationnel $`p_1=0.805`$. Voir la
[feuille de route resserrée](26_FEUILLE_DE_ROUTE_PSTAR.md).

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

### Pourquoi ce cas est favorable

Trois mécanismes complémentaires le justifient.

1. Une paire macroscopiquement éloignée ne peut pas fusionner à une distance
   fixe sous le seuil avec une probabilité non négligeable.
2. À taille et squelette fixés, avancer un niveau postcritique vers
   $`\beta_c`$ améliore exactement le canal au sens de Blackwell.
3. Si les deux sommets ne sont même pas connectés à $`\beta=1`$, les
   recolorations indépendantes de leurs racines effacent exactement leur
   parité.

Ainsi, l'oracle critique concentre les situations les plus favorables à la
conservation de la corrélation : connexion aussi précoce que la géométrie
macroscopique le permet, marques résiduelles de meilleure qualité et absence
d'effacement par racines séparées.

## 3. Théorème favorable réellement établi

Il faut distinguer deux énoncés.

### Théorème F — criticalisation à squelette fixé, établi

Fixons le corridor d'une paire, toutes ses incidences, ses tailles de coupe
$`m_r`$, ses états de bord et une loi arbitrairement corrélée de ses parités
latentes. Pour chaque niveau postcritique $`t_r\ge\beta_c`$, posons

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
rigoureusement :

```math
\boxed{
\text{à géométrie fixée, la fusion critique est le cas postcritique
le plus favorable.}
}
\qquad\text{(3.3)}
```

Références : [19](19_FAVORABLE_SWEEP_PROJECTIONS.md) et
[20](20_COLLAPSED_CORRIDOR_BLACKWELL.md).

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
robuste cherche plutôt une contraction uniforme après criticalisation du
squelette réellement observé.

## 4. Les lemmes structurants

| code | lemme | statut | référence |
|---|---|---|---|
| L0 | $`\beta_{ij}>1`$ implique un effacement exact par racines indépendantes | établi | [19](19_FAVORABLE_SWEEP_PROJECTIONS.md), [22](22_LCA_VS_FULL_HIERARCHY.md) |
| L1 | une paire lointaine fusionnant avant $`\beta_c-\delta`$ a une masse asymptotiquement nulle | établi à $`\delta`$ fixé | [12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md), [14](14_CRITICAL_COMPONENT_BOUNDARY.md) |
| L2 | les arêtes internes aux enfants ne votent pas dans le $`\Lambda_u`$ courant | établi exactement | [14](14_CRITICAL_COMPONENT_BOUNDARY.md) |
| L3 | conditionnellement à $`\Pi_\beta`$, les marques de frontière sont i.i.d. résiduelles | établi exactement | [14](14_CRITICAL_COMPONENT_BOUNDARY.md), [25](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L4 | la charge d'une coupe instantanée est gouvernée par $`m h_p(\beta)^2`$ | établi | [25](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L5 | une coupe fusionnante possède une correction Palm par l'arête gagnante | établi exactement | [09](09_CRITICAL_MERGER_ORACLE.md), [25](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L6 | la loi LCA-Palm repondère une coupe par $`mN_\rho`$ | établi en volume fini | [25](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |
| L7 | à taille fixée, le canal critique Blackwell-domine le canal tardif | établi | [19](19_FAVORABLE_SWEEP_PROJECTIONS.md) |
| L8 | les dégradations se tensorisent sur un corridor fixé, même sous un prior corrélé | établi | [20](20_COLLAPSED_CORRIDOR_BLACKWELL.md) |
| L9 | le corridor collapsed est au plus persistant qu'un sweep des mêmes nœuds | établi | [20](20_COLLAPSED_CORRIDOR_BLACKWELL.md), [22](22_LCA_VS_FULL_HIERARCHY.md) |
| L10 | une chaîne de cactus critique perd exponentiellement la corrélation avec la profondeur | établi exactement | [21](21_CACTUS_COLLAPSED_CERTIFICATE.md) |
| L11 | beaucoup de coupes screenées uniformément contractantes impliquent une obstruction pairwise | conditionnel à la géométrie | [23](23_OPTIMAL_WEAK_RECOVERY_OBSTRUCTION.md), [24](24_SIMPLE_RESIDUAL_BALANCE_OBSTRUCTION.md) |

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

La variable minimale à contrôler sous Palm critique est

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

Après le premier audit à $`p=0.805`$, une preuve suffisante suivrait la
chaîne suivante.

1. **Réduction favorable.** Les paires précoces ont masse $`o(1)`$, les
   racines distinctes ont persistance nulle et les corridors postcritiques
   sont criticalisés à squelette fixé.
2. **Transfert complet.** Une cellule T2-Kruskal conserve la fusion, les
   ports, les attaches, le potentiel extérieur et les
   $`\Lambda_v^{ab}`$ ancestraux.
3. **Composition tordue.** Après une normalisation commune, chaque transition
   porte un déficit $`d_r(z,z')\ge0`$ dans une représentation de
   Feynman--Kac.
4. **Lemme géométrique annealed.** Sous la loi marquée du corridor final, le
   déficit cumulé diverge, sans supposer un screening uniforme.
5. **Clôture.** La composition donne

```math
\mathbb E[H_{\mathcal C}^2]
\le
\mathbb E\left[
\exp\left(-\sum_{r=1}^{N_L}d_r(Z_{r-1},Z_r)\right)
\right]+o(1)
\longrightarrow0.
\qquad\text{(8.2)}
```

Le théorème pairwise des fichiers [03](03_HIERARCHICAL_WEAK_RECOVERY.md) et
[20](20_COLLAPSED_CORRIDOR_BLACKWELL.md) transforme alors (8.1) en absence de
weak recovery.

## 9. Ordre de travail

### Priorité 1 — composition finie tordue

Pour un transfert positif levé, séparer le noyau de masse $`K_r`$ du secteur
signé $`U_r`$, puis exploiter $`|U_r|\le K_r`$. Le but est une formule de
Feynman--Kac avec déficit dépendant de la transition, sans demander une
marge uniforme sur tous les messages.

### Priorité 2 — cellule T2-Kruskal

Étendre la cellule E1+ certifiée au cas d'une fusion réelle avec partition
ouverte, trois ports au moins, une attache en peigne, potentiel extérieur et
quatre $`\Lambda_v^{ab}`$. Deux implémentations indépendantes doivent
coïncider avant toute optimisation.

### Priorité 3 — diagnostic Palm enrichi

Le corridor final réalisé utilise le poids $`N_\rho`$ sur ses nœuds ; le
facteur $`m`$ est déjà dans la course de Kruskal. Le poids $`mN_\rho`$
s'applique aux coupes candidates pré-saut, dans la même fenêtre de rang.
Ajouter aux petites coupes déjà comptées les ports, contournements et boîtes
de potentiel requis par T2.

### Priorité 4 — abondance du déficit

Prouver seulement que l'espérance exponentielle de (8.2) s'annule sous la
loi marquée réelle. Une loi limite complète du dendrogramme et un gap
spectral global ne sont pas requis.

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
| « même composante à $`\beta_c`$ » suffit à prouver la weak recovery | faux : il faut contrôler la masse de l'événement ou établir une domination favorable |
| le LCA critique est globalement le plus favorable sur la grille | ouvert : établi seulement à squelette fixé et sur le cactus |
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

## 12. Parcours de lecture recommandé

1. [README](README.md) — intuition, statuts et navigation.
2. [01 — cadre mathématique](01_MATHEMATICAL_FRAMEWORK.md) — mesure jointe et
   heat baths.
3. [03 — weak recovery](03_HIERARCHICAL_WEAK_RECOVERY.md) — critère pairwise.
4. [25 — information des coupes](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md)
   — géométrie conditionnée et loi Palm.
5. [08 — chaîne ancestrale](08_ANCESTRAL_LAMBDA_CHAIN.md) — calcul des
   $`\Lambda_v^{ab}`$.
6. [20 — corridor collapsed](20_COLLAPSED_CORRIDOR_BLACKWELL.md) — dynamique
   prioritaire et tensorisation.
7. [21 — certificat cactus](21_CACTUS_COLLAPSED_CERTIFICATE.md) — premier
   modèle exact.
8. [05 — feuille de route](05_PROOF_ROADMAP.md) — dépendances techniques.
9. [27 — sous-feuille P0805](27_SUBROADMAP_CORRIDOR_P0805.md) — ordre
   falsifiable et portes go/no-go.
10. [28 — premiers résultats](28_FIRST_CORRIDOR_P0805_RESULTS.md) — audits
    Palm, cellule E1+ et prochaine cellule T2-Kruskal.

Les calculs et commandes de validation sont documentés dans
[computations/README.md](computations/README.md).
