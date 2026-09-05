# La hiérarchie et ses lois

## 1. Le modèle et l'horloge

On se limite ici à deux communautés et à un prior i.i.d. uniforme.
L'observation fixe le graphe et ses poids signés finis. Posons
$`a_e=|W_e|`$ et $`s_e(\sigma)=\mathbf1_{\{W_e\sigma_i\sigma_j>0\}}`$ ;
les arêtes de poids nul sont ignorées. La postérieure est

```math
\pi_O(\sigma)\propto\exp\left(\sum_e a_e s_e(\sigma)\right)
=\text{constante}\times\exp\left(\frac12\sum_{e=\{i,j\}}W_e\sigma_i\sigma_j\right).
```

Sur chaque arête satisfaite, tirer indépendamment
$`\xi_e\sim\mathrm{Exp}(a_e)`$ ; sinon poser $`\xi_e=\infty`$.
Les composantes de $`A_t=\{e:\xi_e\le t\}`$ donnent les partitions emboîtées
$`\Pi_t`$, pour $`0\le t\le1`$. C'est le dendrogramme des diapositives.
**Le niveau $`t`$ est un temps de filtration ; la postérieure ne change pas de température.**

Conditionnellement aux spins, entre deux fusions, deux composantes $`A,B`$
fusionnent au taux $`\sum_{e\in E(A,B)}a_es_e(\sigma)`$.
Le temps jusqu'à la prochaine fusion est exponentiel de taux égal à la
somme de ces taux ; la paire choisie l'est proportionnellement à son taux.
C'est la propriété sans mémoire des horloges. Les fusions successives
ne sont donc pas des variables indépendantes de lois identiques.

Kruskal calcule les partitions. Il ne remplace pas le graphe dans l'énergie :
**toutes les arêtes entre deux amas comptent.**

## 2. La loi exacte à une coupe

Fixer $`t`$ avant de tirer les horloges. Ne conserver que $`A_t`$ : les temps
exacts et les événements futurs sont intégrés. L'identité élémentaire
de gel partiel donne

```math
\pi_O(\sigma\mid A_t)
\propto
\mathbf1_{\{s_e(\sigma)=1\ \forall e\in A_t\}}
\exp\left((1-t)\sum_{e\notin A_t}a_es_e(\sigma)\right).
```

**Preuve.** Une arête ouverte contribue
$`e^{a_es_e}(1-e^{-ta_e})s_e=(e^{a_e}-e^{(1-t)a_e})s_e`$ ;
une arête fermée contribue $`e^{a_es_e}e^{-ta_es_e}=e^{(1-t)a_es_e}`$.
Multiplier ces facteurs donne la formule. C'est une application de la
[représentation d'Edwards–Sokal](https://doi.org/10.1103/PhysRevD.38.2009).

Choisir dans chaque composante $`C`$ une configuration compatible $`g_i`$,
déterminée par les signes des arêtes ouvertes. Toute configuration admissible
s'écrit $`\sigma_i=g_i z_C`$, avec $`z_C\in\{-1,1\}`$.
La loi des orientations est un Ising sur les amas :

```math
\pi_O(z\mid A_t)\propto\exp\left(\sum_{C\lt D}J_{CD}z_Cz_D\right),
\qquad
J_{CD}=\frac{1-t}{2}\sum_{\substack{e=\{i,j\}\\i\in C,\ j\in D}}W_e g_i g_j.
```

À une coupe intermédiaire, **les amas restent couplés**.
Les signes positifs et négatifs de cette somme peuvent se compenser.

## 3. Un pas de dynamique

1. Choisir $`t`$, déterministe ou aléatoire indépendamment des spins à observation fixée.
2. Tirer $`A_t`$ depuis les spins courants.
3. Parcourir ses composantes dans un ordre fixé par leurs sommets, indépendant des spins.
4. Pour chaque composante $`C`$, retourner tous ses spins avec la probabilité ci-dessous, recalculée après chaque mise à jour.
5. Oublier $`A_t`$ et recommencer avec de nouveaux tirages.

```math
\mathbb P(\text{retourner }C\mid\sigma,A_t)
=\frac{1}{1+\exp\left((1-t)\sum_{e=\{i,j\}\in\partial C}W_e\sigma_i\sigma_j\right)}.
```

Ici $`\partial C`$ contient les arêtes ayant exactement une extrémité dans $`C`$.
Cette probabilité est le poids de la configuration retournée divisé par
la somme des deux poids : chaque mise à jour préserve donc la conditionnelle.
Leur composition, précédée du tirage auxiliaire, préserve $`\pi_O`$.
Un balayage ordonné est invariant, sans être nécessairement réversible.

À $`t=0`$, les composantes sont les singletons : c'est exactement un balayage
de Glauber. À $`t=1`$, chaque composante est retournée avec probabilité $`1/2`$,
indépendamment des autres : c'est exactement Swendsen–Wang signé.
À poids finis, l'événement $`A_t=\varnothing`$ a une probabilité positive et
permet d'atteindre toute configuration en un balayage ; la chaîne est donc
irréductible et apériodique. Cela ne donne aucune borne utile sur son temps de mélange.

On peut composer plusieurs niveaux, **en rafraîchissant l'auxiliaire à chaque niveau**.
Choisir une coupe d'après les fusions observées exige une nouvelle preuve de sa loi.

## 4. Si l'on veut garder le dendrogramme non marqué

Notons $`D_{\le t}`$ les partitions et dates de fusion jusqu'à $`t`$, sans
l'identité de l'arête gagnante. Pour une fusion $`u`$ à la date $`b_u`$,
$`\Lambda_u`$ somme les poids satisfaits entre ses deux fils.
$`\Lambda_\infty`$ somme ceux entre les racines de $`\Pi_t`$.
Sur une topologie donnée et des dates compatibles, la densité exacte est

```math
\mathbb P(dD_{\le t}\mid\sigma)
=e^{-t\Lambda_\infty(\sigma)}
\prod_u\Lambda_u(\sigma)e^{-b_u\Lambda_u(\sigma)}\,d\mathbf b.
```

Chaque arête appartient à une unique coupe de fusion, ou relie deux racines.
Dans une coupe, aucune horloge satisfaite ne sonne avant $`b_u`$ ; l'une
sonne à cette date. Cela explique chacun des facteurs, y compris la survie
entre racines. Bayes donne alors

```math
\pi_O(\sigma\mid D_{\le t})\propto
e^{(1-t)\Lambda_\infty(\sigma)}
\prod_u\Lambda_u(\sigma)e^{(1-b_u)\Lambda_u(\sigma)}.
```

Retourner une racine conserve tous les facteurs internes : on retrouve
la probabilité de la section 3. Pour retourner un sous-amas interne,
il faut aussi recalculer les facteurs de ses ancêtres affectés.
Le dendrogramme non marqué et les liens ouverts sont deux auxiliaires
distincts ; ils donnent ici le même noyau de retournement des racines.
