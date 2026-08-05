# Statut scientifique, no-go et vérifications

## 1. Résultat central, avec son périmètre

Le dossier établit la chaîne exacte suivante sur le benchmark de broadcast :

```math
\text{deux répliques postérieures}
\ \longrightarrow\
\text{deux dendrogrammes indépendants}
\ \longrightarrow\
\text{deux marginalisations}
\ \longrightarrow\
\theta^2\text{ par branche}
\ \longrightarrow\
\beta_\chi=\beta_c^{\mathrm{geom}}
\ \Longleftrightarrow\
d\theta^2=1.
\qquad\text{(1.1)}
```

Il aligne cette représentation avec le théorème connu du SBM fini. Il ne
ferme pas encore une nouvelle preuve du théorème fini fondée uniquement sur
le mélange d'un sweep hiérarchique.

## 2. Tableau des énoncés

| énoncé | statut | justification |
|---|---|---|
| postérieure Ising finie avec port de magnétisation | **identité exacte** | calcul de Bayes |
| vérité remplaçable par une réplique postérieure | **théorème** | Nishimori, chapitre 11 |
| critère $Q_n$ à deux répliques pour un avantage positif | **proposition établie** | preuve spectrale courte |
| $\Pi_1$ égale la partition Swendsen--Wang | **identité exacte** | couplage exponentiel |
| couper à tout $\beta$ puis sommer tous les facteurs conserve le Gibbs | **identité exacte** | associativité du sum-product |
| deux coupes indépendantes donnent $\theta^2$ | **identité exacte sur une arête** | produit de deux moyennes |
| $\beta_\chi=\beta_c^{\mathrm{geom}}$ si et seulement si $d\theta^2=1$ | **identité exacte de calibration** | $q_{\beta_\chi}=\theta^2$ et $q_{\beta_c^{\mathrm{geom}}}=1/d$ |
| une coupe partagée donne une valeur strictement plus grande | **identité exacte sur une arête** | moyenne d'un carré |
| le double broadcast a le seuil $d\theta^2=1$ | **théorème classique** | sandwich global / reconstruction sur arbre |
| le SBM binaire fini a le même seuil de weak recovery | **théorème classique** | Mossel--Neeman--Sly, Massoulié |
| le port global fini s'élimine par convolution | **identité exacte en volume fini** | messages de magnétisation |
| remplacer le port par des racines indépendantes au seuil | **réfuté pour la route full-$D$ étudiée** | les deux grandes orientations sont asymptotiquement opposées ; le port doit rester |
| un nombre explicite de sweeps hiérarchiques contracte si $\lambda\le1$ | **ouvert** | mélange/spectral gap non contrôlé |
| $\mathbb E[e^{-\Delta_v/2}]=\rho_n^{n-1}$ pour le sweep à $\beta=0$ | **identité exacte** | [07, prop. 3.1](07_SEUILS_PAR_LA_DYNAMIQUE.md) |
| seuil de l'arbre : reconstruction ssi $\lambda>1$ (hors cas critique) | **démontré depuis zéro** | [08, th. I.17–I.18](08_PREUVES_COMPLETES_SEUILS.md) ; cas $\lambda=1$ cité (E2) |
| almost exact impossible si $\lambda_n\not\to\infty$ | **démontré depuis zéro** | [08, th. II.6](08_PREUVES_COMPLETES_SEUILS.md) (Le Cam + oracle par site) |
| la vérité est stable sous un sweep de Glauber ssi $(\sqrt A-\sqrt B)^2>2$ | **établi** (stabilité) / **établi sous références** (instabilité) | [07, prop. 4.1 et 4.3](07_SEUILS_PAR_LA_DYNAMIQUE.md) |
| la hiérarchie atteint l'exposant almost exact optimal | **programme de recherche** | lift Hellinger écrit, fermeture absente ; stabilité par sommet établie ([07 §5](07_SEUILS_PAR_LA_DYNAMIQUE.md)) |
| la hiérarchie atteint le seuil exact recovery | **programme de recherche** | queue $o(1/n)$ non prouvée ; Glauber depuis un départ froid : ouvert ([07, SBM-DYN0](07_SEUILS_PAR_LA_DYNAMIQUE.md)) |

## 3. Trois réplications à ne pas confondre

### Réplication postérieure

```math
(\sigma^{(1)},D^{(1)}),
(\sigma^{(2)},D^{(2)})
\overset{\mathrm{i.i.d.}}{\sim}
\nu_A
\quad\text{conditionnellement à }A.
\qquad\text{(3.1)}
```

Elle calcule le carré d'overlap de weak recovery.

### Réplication dynamique

Le même environnement $(A,\sigma,D)$ est conservé et deux aléas de sweep
sont tirés indépendamment. Elle calcule le second moment conditionnel d'un
noyau, pas le carré postérieur.

### Réplication Hellinger

Deux hypothèses $X_v=+1$ et $X_v=-1$ possèdent deux fonctions de partition
hiérarchiques. Chacune est marginalisée séparément avant de prendre la
moyenne géométrique. Ce ne sont pas deux répliques de la même postérieure.

## 4. No-go démontrés ou esquissés

(A est démontré par calcul exact ; B est un mécanisme correct mais
esquissé, non rédigé en détail ; C–F sont des garde-fous immédiats au
niveau des définitions.)

### No-go A — dendrogramme partagé

Pour $q=q_\beta$,

```math
\sum_b\pi_bc_b^2
=
\theta^2+\frac{q(1-\theta)^2}{1-q}
>
\theta^2.
\qquad\text{(4.1)}
```

Une preuve fondée sur cette quantité ne vise pas le carré postérieur exact.

### No-go B — dendrogramme complet figé

Sur un arbre clairsemé, les fusions à une arête fixent les parités. Le Gibbs
conditionnel à $D$ voit alors $d\theta$, pas $d\theta^2$.

### No-go C — blocs rendus indépendants après la coupe

Les facteurs postcritiques et les ports recouplent les blocs. Les supprimer
revient à remplacer la postérieure par un canal dégradé non justifié.

### No-go D — $\beta$ interprété comme température

Le temps des horloges ne tempère pas la vraisemblance. Une vraie postérieure
tempérée hors du point bayésien ne satisfait plus l'identité vérité--réplique
du chapitre 11.

### No-go E — $\beta=0$ déclaré égal à Glauber

Sous la projection de coupe, $B_0$ est vide : après marginalisation du
reste de $D$, une mise à jour séquentielle par les conditionnelles de
$\mu_A$ est bien Glauber. En revanche, couper à zéro un dendrogramme complet
déjà conditionné conserve ses temps et facteurs ancêtres ; ses feuilles ne
sont pas des mises à jour de Glauber pour $\mu_A$. Sous une bisection
exacte, il faut en outre des mises à jour équilibrées par paires.

### No-go F — overlap global utilisé pour l'exact recovery

$Q_n\to1$ ne voit pas un nombre sous-linéaire d'erreurs. L'exact recovery
demande une queue locale assez forte pour éliminer simultanément toutes les
erreurs.

## 5. Test falsifiable minimal

Le point scalaire le plus important se vérifie sur
$d=3,\theta=1/2$. On doit obtenir

```math
\lambda_{\mathrm{KS}}=d\theta^2=0.75,
\qquad
d\eta_{\mathrm{partagée}}=1.125.
\qquad\text{(5.1)}
```

Toute implémentation qui renvoie $1.125$ comme transfert postérieur utilise
une coupe partagée ou conserve une information auxiliaire commune.

Le script existant est
[`sbm_critical_cut_replica_diagnostic.py`](../hierarchical-swendsen-wang/computations/sbm_critical_cut_replica_diagnostic.py).

Commande :

```bash
python3 research/hierarchical-swendsen-wang/computations/sbm_critical_cut_replica_diagnostic.py \
  --degree 3 \
  --theta 0.5
```

## 6. Vérifications reproductibles

Les identités SBM existantes sont couvertes par les tests suivants :

```bash
PYTHONDONTWRITEBYTECODE=1 \
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_sbm*.py' \
  -v
```

Ces tests vérifient notamment :

- les formules de $\beta_c^{\mathrm{geom}}$, de corrélation résiduelle et
  d'inflation par
  coupe partagée ;
- la densité d'évolution scalaire du broadcast ;
- les affinités de Bhattacharyya et les trois régimes de recovery ;
- la convolution exacte du port global ;
- les identités finies des diagnostics à deux dendrogrammes.

Ils ne certifient ni le mélange du noyau, ni le passage arbre--graphe, ni un
nouveau seuil hiérarchique.

Les fichiers principaux sont :

- [`sbm_tree_threshold_proofs.py`](../hierarchical-swendsen-wang/computations/sbm_tree_threshold_proofs.py) ;
- [`sbm_broadcast_density_evolution.py`](../hierarchical-swendsen-wang/computations/sbm_broadcast_density_evolution.py) ;
- [`sbm_critical_cut_replica_diagnostic.py`](../hierarchical-swendsen-wang/computations/sbm_critical_cut_replica_diagnostic.py) ;
- [`sbm_global_port_convolution.py`](../hierarchical-swendsen-wang/computations/sbm_global_port_convolution.py) ;
- [`sbm_recovery_regimes_diagnostic.py`](../hierarchical-swendsen-wang/computations/sbm_recovery_regimes_diagnostic.py).

Pour vérifier le Markdown et les formules du dépôt :

```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
```

## 7. Portes avant toute revendication nouvelle

Une preuve dynamique complète sur le SBM fini devra fermer, dans cet ordre :

1. **Invariance :** tous les facteurs et le port global sont présents.
2. **Réplication :** les deux augmentations sont indépendantes
   conditionnellement à $A$.
3. **Géométrie :** les cycles et les intersections de partitions sont
   contrôlés au-delà de la limite locale fixe.
4. **Linéarisation :** l'opérateur overlap a bien rayon $d\theta^2$.
5. **Non-linéaire :** les grands messages sont dominés avec la même
   frontière.
6. **Dynamique :** un nombre explicite de sweeps efface l'état initial.
7. **Fermeture :** cette contraction implique $Q_n\to0$ sous le seuil.
8. **Queues :** pour almost/exact, le fonctionnel Hellinger donne le taux
   local requis.

Avant la huitième porte, il est légitime de parler de représentation exacte,
de calibration ou de programme de preuve, mais pas d'un nouveau théorème
d'achievability hiérarchique.
