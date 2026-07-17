# Bilan d'information ancestral au LCA critique

Ce fichier remplace une tentative trop rapide de définir une « libre entropie
autoduale hiérarchique » en tronquant le produit des facteurs
$`\Lambda_v e^{(1-\beta_v)\Lambda_v}`$. Une telle troncature n'est pas une
marginalisation et ne définit pas, à elle seule, une expérience statistique.

La formulation correcte est un **bilan d'information**. Pour le produit de
spins d'une paire, on révèle successivement :

1. l'environnement de la mise à jour LCA, notamment les orientations non
   rafraîchies ;
2. les données du nœud $u$ ;
3. les données de chaque ancêtre $`v\succ u`$.

Les magnétisations conditionnelles forment alors une martingale. Leur énergie
$L^2$ et leur entropie se décomposent exactement, étage par étage. Les quatre
$`\Lambda_v^{ab}`$ servent à calculer les incréments de cette filtration,
après marginalisation des étages encore cachés.

Le résultat principal de recul est le suivant.

> **Dichotomie exacte.** Si l'on oublie complètement le dendrogramme $D$, sa
> loi jointe redonne la postérieure et tous les facteurs hiérarchiques
> disparaissent. Si l'on conserve $D$, la règle LCA utilise bien tous les
> $`\Lambda_v`$, mais bénéficie d'une information oracle qu'il faut mesurer.

Cette dichotomie n'empêche pas d'utiliser le cas favorable demandé. Elle
précise sa logique : pour une borne d'impossibilité, une expérience plus
informative est légitime ; pour identifier un seuil exact, son avantage doit
être comptabilisé.

## 1. La filtration exacte d'une parité LCA

Fixons un volume fini et une paire déterministe $`i,j`$. Sous la loi
bayésienne jointe, $`\Sigma\mid O`$ a pour loi $`\mu_O`$ ; après cela, tirons
$`D\mid(\Sigma,O)`$ par les horloges. La variable cible est

```math
Y_{ij}:=\Sigma_i\Sigma_j\in\{-1,+1\}.
```

Cette convention évite une ambiguïté de signe : le premier étage du bilan est
alors exactement la corrélation postérieure usuelle. Supposons, pour
l'instant, que $`i,j`$ appartiennent au même arbre et notons
$`u=\mathrm{LCA}_D(i,j)`$.

Le heat bath de $u$ conserve la tribu orbitale
$`\mathscr G_{ij}^{\mathrm{LCA}}`$ du fichier 06 : le dendrogramme non marqué,
les spins hors des deux fils et les relations internes à chacun des deux
fils, mais pas leurs deux orientations globales. On choisit une factorisation
mesurable

```math
\mathscr G_{ij}^{\mathrm{LCA}}
=
\sigma(O,\mathscr E_{ij},A_0,\ldots,A_H),
```

où $`A_0`$ contient la donnée de $u$, $`A_k`$ celle de son $k$-ième ancêtre
strict, et $`\mathscr E_{ij}`$ contient le reste de l'environnement ainsi que
les variables de sélection nécessaires pour décoder cette chaîne. On peut
prendre $`H\le n-2`$ déterministe en complétant la chaîne par un symbole
cimetière. Cette écriture est une **factorisation de l'information**, pas une
hypothèse d'indépendance. En particulier, toute information révélée par le
choix aléatoire de $u$, sa profondeur ou son squelette est comptée dans
$`\mathscr E_{ij}`$ ou dans un $`A_k`$ ; elle ne disparaît pas du bilan.

L'identité de l'arête gagnante reste marginalisée, conformément à la
convention du fichier 08. Posons

```math
\mathscr F_{-1}=\sigma(O),
\qquad
\mathscr F_0=\sigma(O,\mathscr E_{ij}),
```

et, pour $`0\le k\le H`$,

```math
\mathscr F_{k+1}
=
\mathscr F_k\vee\sigma(A_k).
```

Définissons

```math
m_k:=\mathbb E[Y_{ij}\mid\mathscr F_k],
\qquad
L_k:=\log\frac{1+m_k}{1-m_k},
```

avec $`L_k=\pm\infty`$ lorsque $`m_k=\pm1`$. Au dernier étage,

```math
m_{H+1}
=
\mathbb E[Y_{ij}\mid\mathscr G_{ij}^{\mathrm{LCA}}].
```

Après sommation des deux orientations absolues, son signe dépend du choix de
la configuration de référence, mais son carré n'en dépend pas et vaut
exactement

```math
m_{H+1}^2
=
\tanh^2(L_u/2)
=
\eta_u,
```

où $`L_u`$ est le log-rapport complet des fichiers 08 et 12.

## 2. Théorème du bilan ancestral

Posons

```math
\mathsf h(m)
:=
h_2\left(\frac{1+m}{2}\right).
```

### Théorème 2.1 — bilan $L^2$ et entropique, statut : établi

La suite $`(m_k)_{-1\le k\le H+1}`$ est une martingale bornée. Pour tout
$`-1\le k\le H`$,

```math
\boxed{
\mathbb E[m_{k+1}^2]-\mathbb E[m_k^2]
=
\mathbb E[(m_{k+1}-m_k)^2].
}
\tag{L2-k}
```

De plus,

```math
\boxed{
\mathbb E[\mathsf h(m_k)]
-\mathbb E[\mathsf h(m_{k+1})]
=
I(Y_{ij};\mathscr F_{k+1}\mid\mathscr F_k),
}
\tag{H-k}
```

où l'information mutuelle est mesurée en bits. Par sommation,

```math
\boxed{
\mathbb E[m_{H+1}^2]
=
\mathbb E[m_{-1}^2]
+\Delta_u^{\mathrm{env}}
+\sum_{k=0}^{H}\Delta_{u,k}^{(2)},
}
\tag{LEDGER-L2}
```

avec

```math
\Delta_u^{\mathrm{env}}
:=
\mathbb E[(m_0-m_{-1})^2],
\qquad
\Delta_{u,k}^{(2)}
:=
\mathbb E[(m_{k+1}-m_k)^2].
```

L'analogue entropique est

```math
\boxed{
H(Y_{ij}\mid O)
-H(Y_{ij}\mid\mathscr F_{H+1})
=
I(Y_{ij};\mathscr E_{ij}\mid O)
+\sum_{k=0}^{H}I(Y_{ij};A_k\mid\mathscr F_k).
}
\tag{LEDGER-H}
```

### Preuve

La propriété de tour donne

```math
\mathbb E[m_{k+1}\mid\mathscr F_k]=m_k.
```

Ainsi $`m_{k+1}-m_k`$ est orthogonal dans $L^2$ à toute variable
$`\mathscr F_k`$-mesurable, en particulier à $`m_k`$. En développant
$`m_{k+1}^2`$, on obtient (L2-k). L'identité (H-k) est la définition de
l'information conditionnelle écrite comme diminution moyenne de l'entropie
postérieure. Les deux formules télescopiques suivent.

### Interprétation pour la weak recovery

Sous la loi jointe non conditionnée, le premier terme de (LEDGER-L2) est la
vraie quantité postérieure de la paire :

```math
\mathbb E[m_{-1}^2]
=
\mathbb E[c_{ij}(O)^2].
```

Le score LCA terminal contient trois avantages supplémentaires :

1. $`\Delta_u^{\mathrm{env}}`$, dû aux autres orientations de la réplique et
   au dendrogramme hors chaîne ;
2. $`\Delta_{u,0}^{(2)}`$, dû à la donnée du LCA lui-même ;
3. $`\sum_{k\ge1}\Delta_{u,k}^{(2)}`$, dû aux données de tous ses ancêtres
   stricts.

On obtient donc exactement

```math
\boxed{
\text{score LCA local}
=
\text{corrélation postérieure au carré}
+\text{fuite environnementale}
+\text{fuite au LCA et ancestrale}.
}
```

Pour une preuve d'impossibilité, majorer la vraie corrélation par le score LCA
reste valide. Pour une condition nécessaire et suffisante ou une dérivation
du seuil de Nishimori, les fuites ne peuvent pas être ignorées.

Si l'on travaille sous une loi favorable $`\mathbb P^\star`$ conditionnée par
« même arbre et LCA critique », toutes les identités restent vraies sous
$`\mathbb P^\star`$, mais $`m_{-1}`$ devient déjà une magnétisation **oracle**.
On ne peut la comparer à $`c_{ij}(O)`$ sous la loi originale qu'au moyen du
lemme de domination HF ou d'un changement de mesure explicite.

### Corollaire 2.2 — pas d'annulation ancestrale en moyenne

Pour tout préfixe correctement collapsed $`0\le k\le H`$,

```math
\boxed{
\mathbb E[\eta_u]
-\mathbb E[m_k^2]
=
\sum_{\ell=k}^{H}
\mathbb E[(m_{\ell+1}-m_\ell)^2]
\ge0.
}
\tag{no-mean-cancel}
```

Un message ancestral peut annuler ponctuellement le logit local dans une
réalisation donnée, mais révéler des ancêtres ne peut jamais diminuer le score
$L^2$ **en moyenne** lorsqu'on compare deux étages de la même filtration. Le
score obtenu en posant artificiellement $`B_u=0`$ n'est pas un tel étage : il
change l'expérience au lieu de marginaliser les ancêtres cachés.

Cette observation réoriente la preuve. Il ne faut pas espérer une annulation
moyenne des $`\Lambda_v`$ ; il faut montrer que le premier score collapsed est
petit et que la somme des incréments géométriques et exponentiels est
sommable.

## 3. Incréments collapsed et bornes par les logits

Supposons $`|m_k|,|m_{k+1}|<1`$ et posons

```math
\delta_k:=L_{k+1}-L_k.
```

### Proposition 3.1 — contrôle universel d'un étage, statut : établi

En nats,

```math
I(Y_{ij};\mathscr F_{k+1}\mid\mathscr F_k)
=
\mathbb E\left[
D_{\mathrm{KL}}
\left(
\mathrm{Bern}\left(\frac{1+m_{k+1}}2\right)
\middle\|
\mathrm{Bern}\left(\frac{1+m_k}2\right)
\right)
\right].
```

On a les bornes

```math
\boxed{
I(Y_{ij};\mathscr F_{k+1}\mid\mathscr F_k)
\le
\frac18\mathbb E[\delta_k^2],
}
\tag{KL-logit}
```

et

```math
\boxed{
\Delta_{u,k}^{(2)}
\le
\min\left\{
\frac14\mathbb E[\delta_k^2],
2I(Y_{ij};\mathscr F_{k+1}\mid\mathscr F_k)
\right\}.
}
\tag{L2-logit}
```

Dans le second terme du minimum, l'information est en nats.

### Preuve

Pour $`A(\ell)=\log(1+e^\ell)`$, la divergence de Bernoulli est la divergence
de Bregman associée à $A$. Comme

```math
0\le A''(\ell)\le\frac14,
```

Taylor donne la borne $`\delta_k^2/8`$. Ensuite

```math
\left|
\tanh\left(\frac{L_{k+1}}2\right)
-\tanh\left(\frac{L_k}2\right)
\right|
\le\frac12|\delta_k|.
```

La seconde borne vient de Pinsker : pour deux lois de Bernoulli, la variation
totale vaut $`|m_{k+1}-m_k|/2`$.

### Pourquoi les quatre $`\Lambda_v^{ab}`$ restent indispensables

Le logit terminal est calculé par les quatre poids

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v^{ab}e^{(1-\beta_v)\Lambda_v^{ab}}.
```

En revanche, l'incrément $`\delta_k`$ n'est généralement **pas** la seule
différence des termes $`\log\Lambda_{v_k}^{ab}+(1-\beta_{v_k})
\Lambda_{v_k}^{ab}`$. Avant la révélation de $`A_k`$, les ancêtres non révélés
doivent être marginalisés. Le bon calcul est donc :

1. former les quatre poids avec les $`\Lambda_v^{ab}`$ ;
2. sommer exactement les orientations et données encore cachées ;
3. prendre ensuite le log-rapport collapsed ;
4. soustraire les deux logits successifs.

Sur un cactus, ces sommes relèvent d'une récurrence finie. Sur une bande, elles
relèvent d'une matrice de transfert. Remplacer un facteur caché par $1$ ou un
taux par sa moyenne ne produit pas $`\delta_k`$.

### Proposition 3.2 — coût exact d'une horloge exponentielle censurée

Cette proposition isole la brique qui permet de transformer des estimations
des $`\Lambda_v`$ en une borne d'information. Soit
$`R_{\lambda,\tau}`$ l'observation d'une horloge de taux $`\lambda`$ jusqu'à
l'horizon $`\tau>0`$ : on observe son temps si elle sonne avant $`\tau`$, et
un symbole de censure sinon. Pour $`\lambda,\mu>0`$,

```math
\boxed{
D_{\mathrm{KL}}
\left(P_{\lambda,\tau}\middle\|P_{\mu,\tau}\right)
=
(1-e^{-\lambda\tau})
\left(
\log\frac{\lambda}{\mu}
+\frac{\mu}{\lambda}-1
\right).
}
\tag{Exp-KL}
```

Si un bit $`X`$ a une probabilité a priori $`\pi`$ d'être positif et si les
deux taux sont $`\lambda_+,\lambda_-`$, alors, en nats,

```math
\boxed{
I(X;R)
\le
\pi(1-\pi)\,
\mathsf J_\tau(\lambda_+,\lambda_-),
}
\tag{Exp-J}
```

où

```math
\mathsf J_\tau(\lambda,\mu)
:=
D_{\mathrm{KL}}(P_{\lambda,\tau}\|P_{\mu,\tau})
+D_{\mathrm{KL}}(P_{\mu,\tau}\|P_{\lambda,\tau}).
```

La version à quatre taux est tout aussi explicite. Conditionnellement aux
deux parités, supposons que la nuisance d'orientation choisisse les mélanges

```math
P_+
=
w_0^+P_{\Lambda^{00},\tau}
+w_1^+P_{\Lambda^{11},\tau},
```

```math
P_-
=
w_0^-P_{\Lambda^{01},\tau}
+w_1^-P_{\Lambda^{10},\tau}.
```

Dans (4-rate), on abrège

```math
(\Lambda_+^0,\Lambda_+^1)
=(\Lambda^{00},\Lambda^{11}),
\qquad
(\Lambda_-^0,\Lambda_-^1)
=(\Lambda^{01},\Lambda^{10}).
```

Pour tout couplage $`\gamma^+`$ des poids $`w^+,w^-`$ et tout couplage
$`\gamma^-`$ des poids $`w^-,w^+`$,

```math
\begin{aligned}
I(X;R)
\le \pi(1-\pi)\Bigg[&
\sum_{a,b}\gamma^+_{ab}
D_{\mathrm{KL}}
(P_{\Lambda_+^a,\tau}\|P_{\Lambda_-^b,\tau})\\
&+
\sum_{b,a}\gamma^-_{ba}
D_{\mathrm{KL}}
(P_{\Lambda_-^b,\tau}\|P_{\Lambda_+^a,\tau})
\Bigg].
\end{aligned}
\tag{4-rate}
```

On peut minimiser séparément les deux membres droits sur les couplages
$`2\times2`$. Cette borne utilise donc bien les **quatre**
$`\Lambda_v^{ab}`$, sans remplacer un mélange par son taux moyen.

### Preuve

Posons $`\widetilde R=R`$ si l'horloge sonne et
$`\widetilde R=\tau`$ en cas de censure. Sous $`P_{\lambda,\tau}`$, le
log-rapport vaut

```math
\mathbf1_{\{R\ne\dagger\}}\log\frac\lambda\mu
-(\lambda-\mu)\widetilde R.
```

Or

```math
\mathbb P_\lambda(R\ne\dagger)=1-e^{-\lambda\tau},
\qquad
\mathbb E_\lambda[\widetilde R]
=
\frac{1-e^{-\lambda\tau}}\lambda,
```

ce qui donne (Exp-KL). Si $`M=\pi P+(1-\pi)Q`$, la convexité de la divergence
en son second argument donne

```math
D_{\mathrm{KL}}(P\|M)
\le
(1-\pi)D_{\mathrm{KL}}(P\|Q),
```

et l'inégalité symétrique donne (Exp-J). Enfin, la convexité jointe de la
divergence appliquée à un couplage des deux mélanges donne (4-rate).

Si un taux est nul, la formule se lit par limite. En particulier,
$`D_{\mathrm{KL}}(P_{0,\tau}\|P_{\mu,\tau})=\mu\tau`$, tandis que la divergence
inverse est infinie pour $`\mu>0`$. Les coins nuls doivent donc être isolés ;
une borne uniforme fondée seulement sur un contraste relatif ne peut pas les
ignorer.

### Portée exacte pour un ancêtre de Kruskal

Conditionnellement à une coupe candidate fixée et au fait qu'elle a survécu
jusqu'au temps d'entrée, la mémoire sans vieillissement donne précisément le
canal censuré précédent, avec le temps restant pour horizon et
$`\lambda=\Lambda_v^{ab}`$. En revanche, sélectionner la coupe parce qu'elle
est le prochain ancêtre de la paire est une information géométrique
supplémentaire. Elle appartient à $`\mathscr E_{ij}`$ ou au squelette de
$`A_k`$ et doit être comptée séparément. La proposition ne prétend donc pas
que les étages réels sont indépendants.

Enfin, le facteur de vraisemblance **ajouté par le dendrogramme** est
$`\Lambda_v e^{-\beta_v\Lambda_v}`$. Le facteur compensé
$`\Lambda_v e^{(1-\beta_v)\Lambda_v}`$ des $`q_u^{ab}`$ inclut déjà
l'information de $O$ présente à l'étage $`\mathscr F_{-1}`$. Employer ce
dernier facteur comme si toute son exponentielle était une nouvelle fuite de
$D$ compterait deux fois l'observation.

### Corollaire 3.3 — réduction quantitative aux $`\Lambda_v`$

Écrivons chaque donnée ancestrale sous la forme $`A_k=(S_k,R_k)`$, où $`S_k`$
est le sélecteur géométrique de la coupe et $`R_k`$ son horloge censurée.
Conditionnellement à $`(\mathscr F_k,S_k)`$, notons $`\mathcal C_k`$ le membre
droit optimisé de (4-rate), calculé avec les quatre taux et leurs poids
collapsed. Alors

```math
\boxed{
\mathbb E[\eta_u]
\le
\mathbb E[m_0^2]
+2\sum_{k=0}^{H}
\left\{
I(Y_{ij};S_k\mid\mathscr F_k)
+\mathbb E[\mathcal C_k]
\right\}.
}
\tag{Lambda-ledger}
```

Toutes les informations de (Lambda-ledger) sont en nats. La preuve combine la
règle de chaîne, (L2-logit), (4-rate) et le télescopage de (LEDGER-L2).

Cette inégalité sépare exactement les deux verrous :

1. la loi du **squelette sélectionné** le long d'une paire critique, mesurée
   par les termes $`I(Y_{ij};S_k\mid\mathscr F_k)`$ ;
2. les **quatre taux** à coupe fixée, mesurés par $`\mathcal C_k`$.

Une sommabilité des seuls contrastes de taux ne suffit donc pas si la
géométrie sélectionnée révèle elle-même une information macroscopique. En
revanche, une borne sommable sur ces deux familles, jointe à
$`\mathbb E[m_0^2]\to0`$, force le score favorable à s'annuler et ferme la
borne d'impossibilité via HF.

## 4. Théorème de marginalisation : pourquoi $`\widehat\Psi_K`$ n'était pas canonique

La mesure jointe exacte s'écrit

```math
\nu(d\sigma,dD\mid O)
=
\mu(d\sigma\mid O)\,P(dD\mid\sigma,O).
```

### Proposition 4.1 — dichotomie retenir ou oublier, statut : établi

On a

```math
\int\nu(d\sigma,dD\mid O)=\mu(d\sigma\mid O).
```

Par conséquent :

- une dualité appliquée **après marginalisation complète de $D$** est la
  dualité de la postérieure initiale ; la hiérarchie n'ajoute aucun nouveau
  facteur thermodynamique ;
- une dualité ou un heat bath appliqué **conditionnellement à $D$** utilise
  les facteurs $`\Lambda_v`$, mais traite une expérience plus informative ;
- le produit tronqué des $K$ premiers facteurs n'est la densité d'aucune
  marginale canonique tant que l'intégration des facteurs omis et la variable
  conservée ne sont pas spécifiées.

### Preuve

Pour chaque $`(\sigma,O)`$, $`P(\cdot\mid\sigma,O)`$ est une probabilité. Son
intégrale vaut $`1`$. Les trois conséquences sont immédiates.

Cette proposition retire la notation $`\widehat\Psi_K`$ de la feuille de
route. Le défaut de face $`\Psi_0^{\mathrm{face}}`$ du fichier 13 reste exact,
mais son extension hiérarchique doit être indexée par une véritable variable
retenue et accompagnée de (LEDGER-L2) ou (LEDGER-H).

## 5. Contre-audit exact sur une face triangulaire

Considérons les trois bruits $`Z_1,Z_2,Z_3`$ et posons

```math
K:=\#\{e:Z_e=+1\},
\qquad
u_p:=\log\frac p{1-p}.
```

Une arête satisfaite porte une horloge $`\mathrm{Exp}(u_p)`$. À l'instant
$t$, sa probabilité d'ouverture est

```math
a_p(t):=1-e^{-tu_p}.
```

Notons $`C_t`$ l'événement où les trois sommets appartiennent au même arbre à
l'instant $t$. Deux arêtes ouvertes suffisent et sont nécessaires. On obtient
exactement :

| $K$ | syndrome | $`\mathbb P(C_t\mid K)`$ |
|---:|:---:|---:|
| 3 | $+$ | $`c_3(t)=3a_p(t)^2-2a_p(t)^3`$ |
| 2 | $-$ | $`c_2(t)=a_p(t)^2`$ |
| 1 | $+$ | $0$ |
| 0 | $-$ | $0$ |

Posons

```math
r_p=p^3+3p(1-p)^2,
```

```math
\alpha_+
:=\mathbb P(K=3\mid S=+1)=\frac{p^3}{r_p},
```

et

```math
\alpha_-
:=\mathbb P(K=2\mid S=-1)
=\frac{3p^2(1-p)}{1-r_p}.
```

### Proposition 5.1 — information du cas « même arbre », statut : établi

L'information révélée par la seule connexion vaut

```math
\begin{aligned}
I(Z;C_t\mid S)
={}&r_p\left[
h_2(\alpha_+c_3)-\alpha_+h_2(c_3)
\right]\\
&+(1-r_p)\left[
h_2(\alpha_-c_2)-\alpha_-h_2(c_2)
\right].
\end{aligned}
\tag{IC}
```

Définissons $`Y_t`$ comme suit : $`Y_t=\varnothing`$ si $`C_t`$ échoue ; si
$`C_t`$ a lieu, $`Y_t`$ indique la première paire de singletons fusionnée.
Alors

```math
\begin{aligned}
I(Z;Y_t\mid S)
={}&r_p\left[
h_2(\alpha_+c_3)-\alpha_+h_2(c_3)
\right]\\
&+(1-r_p)\left[
h_2(\alpha_-c_2)+\alpha_-c_2\log_2 3
-\alpha_-(h_2(c_2)+c_2)
\right].
\end{aligned}
\tag{IY}
```

Enfin, $`Y_t`$ est une fonction du dendrogramme non marqué, donc

```math
I(Z;D_t\mid S)\ge I(Z;Y_t\mid S).
```

### Preuve

Conditionnellement à $`S=+1`$, seules les classes $`K=3`$ et $`K=1`$ sont
possibles ; seule la première peut réaliser $`C_t`$. Cela donne le premier
crochet de (IC). Le raisonnement est identique pour $`S=-1`$ avec $`K=2`$ et
$`K=0`$.

Lorsque $`K=3`$, la première paire est uniforme parmi trois mais ne distingue
pas davantage l'unique mot $`(+,+,+)`$. Lorsque $`K=2`$, elle est uniforme
parmi les deux arêtes positives. Marginalement sous $`S=-1`$, les trois
paires sont symétriques. Les entropies conditionnelles donnent (IY). La
dernière inégalité est le traitement des données.

### Corollaire 5.2 — fuite d'une parité de paire, statut : établi

Fixons l'arête $`e_1`$. Sur la face isolée et pour une observation fixée,
$`\Sigma_x\Sigma_y=O_{e_1}Z_1`$ ; l'information sur $`Z_1`$ est donc
exactement l'information sur cette parité de spins. Posons

```math
x_+=\frac{1+2\alpha_+}{3},
\qquad
r_+^0
=
\frac{\alpha_+(1-c_3)+(1-\alpha_+)/3}
{1-\alpha_+c_3},
```

et

```math
x_-=\frac{2\alpha_-}{3},
\qquad
r_-^0
=
\frac{2\alpha_-(1-c_2)/3}
{1-\alpha_-c_2}.
```

Alors

```math
\begin{aligned}
I(Z_1;C_t\mid S)
={}&r_p\left[
h_2(x_+)-(1-\alpha_+c_3)h_2(r_+^0)
\right]\\
&+(1-r_p)\left[
h_2(x_-)-\alpha_-c_2h_2(2/3)
-(1-\alpha_-c_2)h_2(r_-^0)
\right],
\end{aligned}
\tag{IC-edge}
```

tandis que

```math
\begin{aligned}
I(Z_1;Y_t\mid S)
={}&r_p\left[
h_2(x_+)-(1-\alpha_+c_3)h_2(r_+^0)
\right]\\
&+(1-r_p)\left[
h_2(x_-)-\frac23\alpha_-c_2
-(1-\alpha_-c_2)h_2(r_-^0)
\right].
\end{aligned}
\tag{IY-edge}
```

En effet, dans le secteur $`S=-1,K=2`$, la connexion seule laisse une
probabilité $`2/3`$ que $`Z_1=+1`$, d'où $`h_2(2/3)`$. Si la première paire
est aussi révélée, l'entropie conditionnelle moyenne sur les trois paires vaut
$`2/3`$ bit : elle vaut $0$ lorsque la première paire est $`e_1`$ et $1$ bit
pour chacune des deux autres. Le secteur plus et l'issue nulle donnent les
autres termes.

Par traitement des données,

```math
I(\Sigma_x\Sigma_y;D_t\mid O)
\ge
I(Z_1;Y_t\mid S)
```

sur la face isolée, après la bijection de jauge du fichier 13.

Le même calcul donne directement le bilan $L^2$. Avant révélation de la
connexion,

```math
M_0
=
r_p(2x_+-1)^2+(1-r_p)(2x_--1)^2.
```

Après révélation de $`C_t`$,

```math
\begin{aligned}
M_C={}&r_p\left[
\alpha_+c_3+(1-\alpha_+c_3)(2r_+^0-1)^2
\right]\\
&+(1-r_p)\left[
\frac{\alpha_-c_2}{9}
+(1-\alpha_-c_2)(2r_-^0-1)^2
\right].
\end{aligned}
```

Après révélation de $`Y_t`$, il suffit de remplacer
$`\alpha_-c_2/9`$ par $`\alpha_-c_2/3`$ ; notons le résultat $`M_Y`$. Les
gains martingale sont exactement

```math
\Delta_C^{(2)}=M_C-M_0,
\qquad
\Delta_Y^{(2)}=M_Y-M_0.
\tag{face-L2}
```

### Valeurs au point de Nishimori et au temps de percolation

Soit

```math
q_c=2\sin(\pi/18),
\qquad
p_{\mathrm N}^{(0)}=0.835805792367\ldots.
```

Le temps critique vérifie

```math
p\,a_p(\beta_c)=q_c.
```

À $`p=p_{\mathrm N}^{(0)}`$ et $`t=\beta_c(p)`$, le calcul exact donne

```math
I(Z;C_{\beta_c}\mid S)
=0.043883918779\ldots\ \text{bit},
```

et

```math
\boxed{
I(Z;D_{\beta_c}\mid S)
\ge
I(Z;Y_{\beta_c}\mid S)
=0.078638140273\ldots\ \text{bit}.
}
\tag{LEAK-face}
```

Pour une paire fixée, le corollaire donne plus directement

```math
I(Z_1;C_{\beta_c}\mid S)
=0.027809400607\ldots\ \text{bit},
```

et

```math
\boxed{
I(\Sigma_x\Sigma_y;D_{\beta_c}\mid O)
\ge
I(Z_1;Y_{\beta_c}\mid S)
=0.042759377412\ldots\ \text{bit}.
}
\tag{LEAK-pair}
```

Les gains $L^2$ correspondants sont

```math
\Delta_C^{(2)}
=0.006320258880\ldots,
\qquad
\Delta_Y^{(2)}
=0.019523088673\ldots.
\tag{LEAK-pair-L2}
```

La comparaison avec la meilleure borne d'impossibilité actuelle du dossier
est instructive :

| paramètre | $`I(Z_1;Y_{\beta_c}\mid S)`$ | $`\Delta_Y^{(2)}`$ | entropie de Palm |
|:---|---:|---:|---:|
| $`p_{\mathrm{info}}=(1+\sqrt{q_c})/2=0.794659275831\ldots`$ | $`0.061716309354`$ | $`0.031834915394`$ | $`0.314600241391`$ |
| $`p_{\mathrm N}^{(0)}=0.835805792367\ldots`$ | $`0.042759377412`$ | $`0.019523088673`$ | $`0.251560120699`$ |

L'oracle de face reste donc quantitativement informatif aux deux repères. Ce
tableau est un **contre-audit local**, pas un minorant pour une paire à
distance divergente : seul le cactus croissant puis la bande peuvent décider
si ce coût se propage ou se dissipe avec l'échelle.

Ce minorant ne compte même pas l'information des forêts partielles ni celle
des temps exacts.

## 6. Loi de Palm d'une fusion exactement critique

L'événement $`\beta_u=t`$ a probabilité nulle. La formulation correcte fixe
la densité du second temps de fusion. Supposons que la première paire fusionne
à l'instant $`s`$ et que le triangle devienne connexe à l'instant $`t`$, avec
$`0<s<t`$.

Pour une première paire fixée, les densités conditionnelles aux mots de bruit
sont :

| syndrome | états compatibles | densité par état |
|:---:|---:|---:|
| $+$ | un état $`K=3`$ | $`2u_p^2e^{-u_ps-2u_pt}`$ |
| $-$ | deux états $`K=2`$ | $`u_p^2e^{-u_p(s+t)}`$ |

### Théorème 6.1 — entropie de l'oracle de Palm, statut : établi

Sous la loi de Palm du second merge à l'instant $t$,

```math
\boxed{
\mathbb P(S=-1\mid\beta_u=t)
=
\frac{1-p}{1-p+p e^{-tu_p}}.
}
\tag{Palm-S}
```

Conditionnellement à l'observation et au dendrogramme non marqué complet,
l'entropie du mot de bruit vaut $0$ bit lorsque $`S=+1`$ et $1$ bit lorsque
$`S=-1`$. Son entropie moyenne sous cette loi de Palm est donc

```math
\boxed{
H_{\mathrm{Palm},t}(Z\mid O,D)
=
\frac{1-p}{1-p+p e^{-tu_p}}\ \text{bit}.
}
\tag{Palm-H}
```

Au temps de percolation triangulaire,

```math
p e^{-\beta_cu_p}=p-q_c,
```

donc

```math
\boxed{
H_{\mathrm{Palm},\beta_c}(Z\mid O,D)
=
\frac{1-p}{1-q_c}\ \text{bit}.
}
```

Au point $`p=p_{\mathrm N}^{(0)}`$,

```math
H_{\mathrm{Palm},\beta_c}(Z\mid O,D)
=0.251560120699\ldots\ \text{bit}.
```

### Preuve

Après multiplication par les probabilités a priori, la masse de l'état plus
est

```math
p^3\,2u_p^2e^{-u_ps-2u_pt},
```

tandis que la masse totale des deux états moins est

```math
2p^2(1-p)u_p^2e^{-u_p(s+t)}.
```

Leur rapport donne (Palm-S), indépendamment de $s$. Le syndrome est observé.
Dans le secteur plus, le mot est unique. Dans le secteur moins, la première
paire est positive et il reste deux choix équiprobables pour la seconde arête
positive. Les entropies sont donc $0$ et $1$ bit. La relation critique donne
la dernière simplification.

### Contre-audit de portée

Au même paramètre, l'entropie de face non conditionnée du fichier 13 vaut
exactement $1$ bit, alors que l'entropie de l'expérience de Palm favorable
vaut $`0.251560\ldots`$ bit. Ces quantités vivent sous deux lois différentes ;
on ne doit pas les soustraire comme une information mutuelle. Leur écart
prouve néanmoins qu'elles ne peuvent pas être identifiées.

Le cas « même arbre, fusion au seuil » est donc réellement beaucoup plus
favorable, comme souhaité. Cette propriété renforce une preuve d'impossibilité
qui réussirait sous cet oracle, mais interdit de réutiliser directement
l'équation d'entropie de Nishimori dans l'expérience conditionnée.

## 7. Nouvelle cible sur cactus

Le premier calcul non ambigu n'est plus un zéro $`\widehat\Psi_1`$. C'est le
bilan (LEDGER-L2) sur un cactus de deux triangles partageant un sommet, avec
une paire dont le LCA possède un ancêtre strict retenu.

Il faut calculer simultanément :

1. $`\mathbb E[m_{-1}^2]`$, la corrélation postérieure réelle ;
2. $`\Delta_u^{\mathrm{env}}`$, l'avantage du sélecteur LCA, du squelette hors
   chaîne et des orientations non rafraîchies ;
3. $`\Delta_{u,0}^{(2)}`$, l'information du merge critique ;
4. $`\Delta_{u,1}^{(2)}`$, l'information du premier ancêtre et donc des quatre
   $`\Lambda_{v_1}^{ab}`$ ;
5. le score terminal $`\mathbb E[\tanh^2(L_u/2)]`$.

Le cactus à deux triangles possède cinq sommets et six arêtes. Il suffit
d'énumérer $`2^6`$ mots de bruit, les types d'ordre de Kruskal compatibles et
les intégrales d'ordre exponentielles. Chaque ligne du bilan doit être
vérifiée de deux manières : par conditionnements successifs et par calcul
direct de la postérieure.

Une fois ce calcul fermé, deux résultats deviennent possibles.

- Si les incréments oracle restent grands alors que la corrélation réelle est
  petite, la dynamique LCA à un pas ne peut pas identifier le seuil exact sans
  débiaisage.
- Si les incréments ancestraux sont sommables uniformément et si le score
  terminal contracte même sous l'expérience favorable, (LEDGER-L2) et HF
  donnent une borne d'impossibilité rigoureuse.

## 8. Statut et ordre de travail

| Élément | Statut |
|---|---|
| Bilan martingale (LEDGER-L2) | **Établi, volume fini** |
| Bilan entropique (LEDGER-H) | **Établi, volume fini** |
| Bornes par incréments de logit | **Établies** |
| KL d'une horloge censurée et borne à quatre taux | **Établis** |
| Dichotomie marginaliser / conditionner par $D$ | **Établie** |
| Information de connexion d'une face | **Établie exactement** |
| Loi de Palm du second merge d'une face | **Établie exactement** |
| Bilan complet sur le cactus à deux triangles | **Prochain calcul** |
| Sommabilité des incréments sous le biais d'une paire critique | **À prouver** |
| Identification du seuil de weak recovery à $`p_{\mathrm N}^{(0)}`$ | **Conjecture** |

Les contre-audits reproductibles sont dans
[`computations/ancestral_information_ledger.py`](computations/ancestral_information_ledger.py)
et
[`computations/test_ancestral_information_ledger.py`](computations/test_ancestral_information_ledger.py).
