# SBM — dendrogrammes répliqués et seuil de weak recovery

Ce dossier transpose au **SBM symétrique binaire** la démarche du
[chapitre 11 du manuscrit](../../Manuscrit_de_thèse.pdf#page=133) :
postérieure bayésienne, rééchantillonnage de Nishimori, dynamique invariante,
puis obstruction informationnelle. L'objet nouveau est un dendrogramme
d'horloges. La coupe appelée $\beta_c$ dans le programme de recherche sera
notée ici $\beta_c^{\mathrm{geom}}$ (ou, plus brièvement,
$\beta_{\mathrm{geom}}$) pour la distinguer du temps informationnel
$\beta_\chi$ introduit ci-dessous.

> **Réponse centrale.** La représentation sans biais du carré d'overlap
> utilise deux répliques postérieures et **deux dendrogrammes indépendants
> conditionnellement à la seule observation**. On coupe chacun à
> $\beta_c^{\mathrm{geom}}$, on conserve tous les facteurs situés au-dessus
> de la coupe, puis on marginalise chaque dendrogramme séparément. Le mode
> d'overlap transmis par une arête est alors exactement $\theta^2$. Un
> dendrogramme partagé représente un autre couplage, plus informatif, et
> donne en général un faux seuil.

Pour

```math
\mathbb P(A_{ij}=1\mid X_iX_j=+1)=\frac an,
\qquad
\mathbb P(A_{ij}=1\mid X_iX_j=-1)=\frac bn,
\qquad a>b>0,
```

on pose

```math
d=\frac{a+b}{2},
\qquad
\theta=\frac{a-b}{a+b},
\qquad
\lambda=d\theta^2=\frac{(a-b)^2}{2(a+b)}.
```

Le seuil à retrouver est

```math
\boxed{\lambda_c=1.}
```

Autrement dit, la weak recovery est impossible pour $\lambda\le1$ et
possible pour $\lambda>1$. La coupe géométrique et le temps informationnel
sont définis sur la **même horloge** par

```math
q_\beta
=
\frac{1+\theta}{2}
\left(
1-\exp\left[
-\beta\log\frac{1+\theta}{1-\theta}
\right]
\right),
\qquad
q_{\beta_c^{\mathrm{geom}}}=\frac1d,
\qquad
q_{\beta_\chi}=\theta^2.
```

La coupe $\beta_c^{\mathrm{geom}}$ rend les blocs géométriquement critiques ;
**elle ne crée pas le carré**. Le carré vient du produit de deux transferts
Gibbs indépendants :

```math
\theta\times\theta=\theta^2,
\qquad
d\text{ branches}\ \Longrightarrow\ d\theta^2.
```

On replace ensuite ce carré sur l'horloge du dendrogramme. Comme $q_\beta$
est croissante,

```math
\boxed{
\begin{aligned}
\beta_\chi<\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2<1,\\
\beta_\chi=\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2=1,\\
\beta_\chi>\beta_c^{\mathrm{geom}}
&\Longleftrightarrow d\theta^2>1.
\end{aligned}
}
```

Ainsi, au seuil de Kesten--Stigum, le temps de l'information rencontre
**exactement** la coupe géométrique. Cette égalité est la lecture demandée du
seuil au niveau $\beta_c$. Elle est distributionnelle : la percolation
d'information de rétention $\theta^2$ a la même loi que la coupe au temps
$\beta_\chi$. Elle ne permet ni de partager un dendrogramme entre les deux
répliques, ni de supprimer les facteurs postcritiques.

Le dossier coupe chaque arbre à $\beta_c^{\mathrm{geom}}$ pour organiser le
calcul, mais conserve tout le corridor postcritique. Le coefficient
$\theta^2$ reste le résultat de deux marginalisations exactes et vaut quel
que soit l'ordre de coupe.

La solution finie définissant $\beta_c^{\mathrm{geom}}$ exige $a>2$ ; elle
appartient à $[0,1]$ exactement lorsque $a-b\ge2$. Cette condition est
satisfaite au seuil de Kesten--Stigum non dégénéré.

![Deux dendrogrammes indépendants, coupés séparément au même beta critique géométrique, puis marginalisés avant de former l'overlap.](figures/deux_dendrogrammes_beta_c.svg)

## Parcours conseillé

1. [Du chapitre 11 au critère à deux répliques](01_DU_CHAPITRE_11_AU_SBM.md)
   fixe le modèle, l'overlap, Nishimori et la quantité quadratique exacte.
2. [Deux dendrogrammes coupés à $\beta_c^{\mathrm{geom}}$](02_DEUX_DENDROGRAMMES_A_BETA_C.md)
   construit la coupe et démontre pourquoi les deux copies doivent être
   indépendantes.
3. [Calibration hiérarchique et preuve sur le broadcast](03_PREUVE_DU_SEUIL_WEAK_RECOVERY.md)
   ferme le calcul sur le broadcast $\mathrm{PGW}(d)$ et précise le passage
   au théorème du SBM fini.
4. [Dynamique hiérarchique](04_DYNAMIQUE_HIERARCHIQUE.md) décrit le noyau
   invariant, les messages de blocs et distingue le full-$D$ de la famille
   de projections reliant Glauber à Swendsen--Wang.
5. [Almost exact et exact recovery](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md)
   explique pourquoi le fonctionnel doit devenir
   Hellinger--Chernoff et pourquoi
   $\beta_{c,n}^{\mathrm{geom}}\to0$ dans les régimes d'almost/exact
   recovery.
6. [Statut scientifique et tests de réfutation](06_STATUT_SCIENTIFIQUE.md)
   sépare les identités établies de la nouvelle preuve dynamique qui reste
   à construire.
7. [Références](REFERENCES.md) regroupe le manuscrit, les articles primaires
   et les calculs reproductibles déjà présents dans le dépôt.

## Trois précisions indispensables

### $\beta$ est un temps d'horloge, pas une température

Une arête satisfaite reçoit une horloge exponentielle. La partition
$\Pi_\beta$ contient les arêtes dont l'horloge a sonné avant $\beta$.
Modifier $\beta$ modifie donc la **résolution du dendrogramme**, pas la
postérieure. Tempérer la postérieure ferait perdre l'identité de Nishimori,
sauf à changer aussi le modèle génératif.

### Couper ne signifie pas tronquer

Les facteurs postcritiques ne sont jamais supprimés. La coupe à
$\beta_c^{\mathrm{geom}}$
sert uniquement à grouper les premières sommes du calcul de Gibbs. Les états
de ports des blocs coupés sont encore tirés conjointement sous tous les
facteurs supérieurs.

### Le résultat exact et la preuve dynamique ne sont pas identiques

Le dossier donne une représentation hiérarchique exacte du seuil sur la
limite locale de broadcast et l'aligne avec le théorème connu du SBM fini.
Il ne prétend pas encore qu'un nombre fixé de sweeps du nouveau noyau
hiérarchique prouve à lui seul le seuil sur le graphe fini. Les non-arêtes
ou la contrainte de balance y forment un port global qui doit rester dans le
Gibbs.

## Résumé en une ligne

```math
\boxed{
\text{même observation}
\ +\
\text{deux Gibbs postérieurs}
\ +\
\text{deux dendrogrammes indépendants}
\ +\
\text{deux marginalisations}
\ =\
d\theta^2.
}
```
