# Échantillonnage exact et direct du random-cluster signé \(q=2\) à partir du graphe observé

## 0. Objet de la note

Cette note décrit un échantillonneur **exact**, **direct** et **sans variable de spin initiale** du random-cluster associé à un modèle d’Ising signé et frustré sur un graphe fini pondéré

\[
G=(V,E,W),
\qquad W_e\in\mathbb R.
\]

L’algorithme utilise uniquement :

- les signes et modules des poids observés ;
- une structure union-find avec parités ;
- des priorités exponentielles indépendantes ;
- des pièces de Bernoulli ;
- un oracle stochastique de masse de complétion, lui-même construit par un passage glouton.

Il n’utilise :

- ni ground truth ;
- ni configuration de spins courante ;
- ni chaîne de Markov ;
- ni burn-in ;
- ni calcul explicite d’une fonction de partition.

Une exécution terminée fournit exactement :

1. un random-cluster signé \(A\) ;
2. un dendrogramme marqué par des temps exponentiels ;
3. si désiré, une configuration de spins distribuée exactement selon la mesure de Gibbs.

L’échantillonneur est de type **Las Vegas** : sa sortie est toujours exacte et il termine presque sûrement sur tout graphe fini. Son coût moyen peut néanmoins être élevé dans les régimes très cycliques ou fortement frustrés.

---

# 1. Convention pour le modèle d’Ising signé

Pour chaque arête \(e=\{u,v\}\), posons

\[
s_e:=\operatorname{sign}(W_e)\in\{-1,+1\},
\qquad
a_e:=|W_e|\ge 0.
\]

Les arêtes de poids nul peuvent être supprimées. On supposera donc \(a_e>0\).

Nous utilisons la convention énergétique

\[
\boxed{
\mu_W(\sigma)
=
\frac1{Z_W}
\exp\left[-\sum_{e=\{u,v\}\in E}
a_e\,
\mathbf 1_{\{\sigma_u\sigma_v\neq s_e\}}
\right],
}
\tag{1.1}
\]

pour \(\sigma\in\{-1,+1\}^V\).

Ainsi :

- \(s_e=+1\) favorise \(\sigma_u=\sigma_v\) ;
- \(s_e=-1\) favorise \(\sigma_u=-\sigma_v\) ;
- le rapport entre le poids d’une arête satisfaite et celui d’une arête non satisfaite vaut \(e^{a_e}\).

Dans la convention standard

\[
\mu(\sigma)
\propto
\exp\left(\beta\sum_{e=\{u,v\}}J_e\sigma_u\sigma_v\right),
\]

il faut prendre

\[
s_e=\operatorname{sign}(J_e),
\qquad
a_e=2\beta|J_e|.
\]

Définissons enfin

\[
\boxed{
p_e:=1-e^{-a_e}\in(0,1).
}
\tag{1.2}
\]

---

# 2. Augmentation d’Edwards--Sokal et loi cible

Introduisons une configuration de liens

\[
A\subseteq E.
\]

La mesure jointe spins-liens est

\[
\boxed{
\Psi_W(\sigma,A)
\propto
\prod_{e=\{u,v\}\in E}
\left[
(1-p_e)\mathbf 1_{\{e\notin A\}}
+
p_e\mathbf 1_{\{e\in A\}}
\mathbf 1_{\{\sigma_u\sigma_v=s_e\}}
\right].
}
\tag{2.1}
\]

Une arête ouverte doit donc être satisfaite. Une arête fermée n’impose aucune contrainte.

## Proposition 2.1 — marginale en spins

La marginale de \(\Psi_W\) en \(\sigma\) est exactement \(\mu_W\).

### Preuve

Fixons \(\sigma\). Pour une arête satisfaite,

\[
(1-p_e)+p_e=1.
\]

Pour une arête non satisfaite, seule l’option fermée est possible et fournit

\[
1-p_e=e^{-a_e}.
\]

Par conséquent,

\[
\sum_{A\subseteq E}\Psi_W(\sigma,A)
\propto
\prod_{e:\,\sigma_u\sigma_v\neq s_e}e^{-a_e}
=
\exp\left[-\sum_ea_e\mathbf 1_{\{\sigma_u\sigma_v\neq s_e\}}\right].
\]

C’est bien la mesure (1.1). \(\square\)

---

# 3. Sous-graphes équilibrés et facteur \(2^{k(A)}\)

## Définition 3.1 — équilibre signé

Un sous-graphe \(A\subseteq E\) est dit **équilibré**, **compatible** ou **non frustré** lorsque le système

\[
\sigma_u\sigma_v=s_e,
\qquad e=\{u,v\}\in A,
\tag{3.1}
\]

admet au moins une solution.

Les propriétés suivantes sont équivalentes :

1. \(A\) est équilibré ;
2. tout cycle \(\gamma\subseteq A\) vérifie
   \[
   \prod_{e\in\gamma}s_e=+1 ;
   \]
3. tout cycle de \(A\) contient un nombre pair d’arêtes négatives ;
4. il existe une jauge \(g:V\to\{-1,+1\}\) telle que
   \[
   s_e=g_ug_v
   \qquad\forall e=\{u,v\}\in A.
   \]

Notons \(k(A)\) le nombre de composantes connexes de \((V,A)\), sommets isolés compris.

## Proposition 3.2 — nombre de colorations compatibles

Si \(A\) est équilibré, le nombre de configurations de spins satisfaisant toutes les contraintes ouvertes vaut exactement

\[
\boxed{2^{k(A)}.}
\tag{3.2}
\]

### Preuve

Dans chaque composante connexe, les signes relatifs sont déterminés par les contraintes. Il reste un unique flip global libre par composante. \(\square\)

En sommant \(\sigma\) dans (2.1), on obtient la loi cible.

## Théorème 3.3 — random-cluster signé \(q=2\)

La marginale en liens est

\[
\boxed{
\phi_W(A)
=
\frac1{Z_{\mathrm{RC}}}
\mathbf 1_{\{A\text{ équilibré}\}}
2^{k(A)}
\prod_{e\in A}p_e
\prod_{e\notin A}(1-p_e).
}
\tag{3.3}
\]

Tout l’objet de la suite est d’échantillonner exactement cette loi à partir de \((G,W)\), sans introduire préalablement une configuration \(\sigma\).

---

# 4. Pourquoi le glouton local n’est pas exact

La règle

> « proposer indépendamment les arêtes, puis rejeter seulement l’arête qui crée un cycle frustré »

produit toujours un sous-graphe équilibré, mais pas la loi (3.3).

Deux effets manquent :

1. chaque fusion de deux composantes doit payer un facteur \(1/2\), car elle diminue \(k(A)\) de un ;
2. l’effet d’une fusion dépend des cycles et chemins qui pourront encore apparaître dans le graphe restant.

Le premier effet est local. Le second est intrinsèquement futur et empêche une règle purement myope d’être exacte sur un graphe cyclique général.

L’algorithme ci-dessous conserve l’esprit glouton, mais corrige exactement les complétions futures par des pièces de Bernoulli.

---

# 5. Ordre de traitement et union-find avec parités

Fixons un ordre

\[
\pi=(e_1,\ldots,e_m),
\qquad m=|E|,
\]

indépendant du random-cluster à échantillonner.

L’ordre peut être déterministe. Pour obtenir une construction hiérarchique liée aux poids, on peut tirer des **priorités auxiliaires**

\[
R_e\sim\operatorname{Exp}(a_e)
\]

indépendamment et trier les arêtes par \(R_e\) croissant.

> Ces priorités servent seulement à choisir l’ordre de décision. Les temps physiques du dendrogramme exact seront tirés séparément après l’échantillonnage de \(A\), à la section 12.

Après avoir traité les \(j\) premières arêtes, l’état courant est un sous-graphe équilibré \(A_j\), représenté par une union-find avec parités.

Pour deux sommets déjà connectés, la structure fournit la relation imposée

\[
\pi_j(u,v)\in\{-1,+1\}.
\]

Pour l’arête suivante \(e=\{u,v\}\), trois cas sont possibles :

1. **fusion** : \(u\) et \(v\) sont dans deux composantes différentes ;
2. **cycle compatible** : ils sont connectés et \(\pi_j(u,v)=s_e\) ;
3. **cycle frustré** : ils sont connectés et \(\pi_j(u,v)\neq s_e\).

---

# 6. Masse normalisée des complétions

Soit

\[
R_j:=\{e_{j+1},\ldots,e_m\}
\]

l’ensemble des arêtes restant à traiter.

Pour un état équilibré \(S_j\), de sous-graphe courant \(A_j\), définissons

\[
\boxed{
H_j(S_j)
:=
\sum_{\substack{B\subseteq R_j\\A_j\cup B\text{ équilibré}}}
2^{k(A_j\cup B)-k(A_j)}
\prod_{e\in B}p_e
\prod_{e\in R_j\setminus B}(1-p_e).
}
\tag{6.1}
\]

La quantité \(H_j(S_j)\) est la masse totale normalisée de toutes les complétions futures.

Comme toute fusion future diminue le nombre de composantes,

\[
2^{k(A_j\cup B)-k(A_j)}\le1.
\]

Par conséquent,

\[
0<H_j(S_j)\le1.
\tag{6.2}
\]

La positivité vient de la complétion où toutes les arêtes restantes sont fermées :

\[
H_j(S_j)
\ge
\prod_{e\in R_j}(1-p_e)>0.
\tag{6.3}
\]

À la fin,

\[
H_m(S_m)=1.
\tag{6.4}
\]

---

# 7. Récursions exactes arête par arête

Considérons l’arête suivante \(e=e_{j+1}=\{u,v\}\), de paramètre \(p_e\).

## 7.1 Cas d’une fusion

Notons :

- \(S^0\) l’état où \(e\) est fermée ;
- \(S^1\) l’état où \(e\) est ouverte et où les deux composantes sont fusionnées avec la contrainte \(g_ug_v=s_e\).

La branche ouverte apporte le facteur \(p_e\), mais aussi le facteur \(1/2\) correspondant à la perte d’une composante. Donc

\[
\boxed{
H_j(S)
=
(1-p_e)H_{j+1}(S^0)
+
\frac{p_e}{2}H_{j+1}(S^1).
}
\tag{7.1}
\]

La probabilité conditionnelle exacte d’ouvrir \(e\) est

\[
\boxed{
q_e(S)
=
\frac{\frac{p_e}{2}H_{j+1}(S^1)}
{(1-p_e)H_{j+1}(S^0)+\frac{p_e}{2}H_{j+1}(S^1)}.
}
\tag{7.2}
\]

## 7.2 Cas d’un cycle compatible

Ouvrir ou fermer \(e\) ne change ni la partition, ni les parités, ni le problème futur. Les deux branches ont donc la même masse de complétion.

Ainsi

\[
\boxed{\mathbb P(e\in A\mid S)=p_e.}
\tag{7.3}
\]

## 7.3 Cas d’un cycle frustré

L’ouverture rendrait le sous-graphe incompatible. Elle a donc une masse nulle :

\[
\boxed{\mathbb P(e\in A\mid S)=0.}
\tag{7.4}
\]

La seule difficulté est donc le cas d’une fusion, où apparaissent les deux masses futures inconnues \(H_{j+1}(S^0)\) et \(H_{j+1}(S^1)\).

---

# 8. Le coin de complétion exact

La quantité \(H_j(S_j)\) n’est jamais calculée numériquement. On construit à la place une variable de Bernoulli dont la probabilité de succès est exactement \(H_j(S_j)\).

## Algorithme `CompletionCoin`

Entrées :

- un état équilibré \(S_j\) ;
- l’indice de la première arête restante ;
- l’ordre \(\pi\).

On travaille sur une copie de la union-find avec parités.

Pour chaque arête restante \(e=\{u,v\}\), dans l’ordre :

### A. Les extrémités sont dans deux composantes différentes

1. Tirer
   \[
   B_e\sim\operatorname{Bernoulli}(p_e).
   \]
2. Si \(B_e=0\), laisser l’arête fermée et continuer.
3. Si \(B_e=1\), tirer
   \[
   C_e\sim\operatorname{Bernoulli}(1/2).
   \]
4. Si \(C_e=0\), retourner immédiatement `0`.
5. Si \(C_e=1\), fusionner les composantes avec la parité \(s_e\), puis continuer.

### B. Les extrémités sont connectées avec la bonne parité

Ne rien tirer et continuer.

En effet, les options ouverte et fermée ont une masse totale

\[
(1-p_e)+p_e=1
\]

et conduisent au même état futur.

### C. Les extrémités sont connectées avec la mauvaise parité

1. Tirer
   \[
   C_e\sim\operatorname{Bernoulli}(1-p_e).
   \]
2. Si \(C_e=0\), retourner immédiatement `0`.
3. Sinon, continuer avec l’arête fermée.

### Fin

Si toutes les arêtes ont été parcourues sans échec, retourner `1`.

## Théorème 8.1 — exactitude du coin

\[
\boxed{
\mathbb P\bigl(\texttt{CompletionCoin}(S_j)=1\bigr)
=
H_j(S_j).
}
\tag{8.1}
\]

### Preuve

On raisonne par récurrence sur le nombre d’arêtes restantes.

- Sans arête restante, le coin retourne \(1\), comme \(H_m=1\).
- Pour une fusion, sa probabilité de succès vaut
  \[
  (1-p_e)H(S^0)+p_e\frac12H(S^1),
  \]
  qui est exactement la récursion (7.1).
- Pour un cycle compatible, le coin passe directement à l’état suivant, ce qui correspond au facteur total \(1\).
- Pour un cycle frustré, le coin survit avec probabilité \(1-p_e\), exactement le poids de l’arête nécessairement fermée.

Les récursions et la condition terminale coïncident. \(\square\)

---

# 9. Course de Bernoulli pour une arête fusionnante

Pour une fusion, posons

\[
c_0:=1-p_e,
\qquad
c_1:=\frac{p_e}{2}.
\]

Nous voulons choisir \(I\in\{0,1\}\) avec

\[
\mathbb P(I=i)
=
\frac{c_iH_i}{c_0H_0+c_1H_1},
\tag{9.1}
\]

avec

\[
H_i:=H_{j+1}(S^i).
\]

## Algorithme `FusionRace`

Répéter :

1. proposer la branche ouverte \(I=1\) avec probabilité
   \[
   \frac{c_1}{c_0+c_1}
   =
   \frac{p_e}{2-p_e}
   =
   \tanh\left(\frac{a_e}{2}\right),
   \tag{9.2}
   \]
   et la branche fermée sinon ;
2. lancer un `CompletionCoin` indépendant dans l’état \(S^I\) ;
3. si le coin retourne \(1\), accepter définitivement la branche \(I\) ;
4. sinon, recommencer la course avec de nouveaux aléas.

## Théorème 9.1 — exactitude de la course

`FusionRace` retourne \(I=i\) avec la probabilité exacte (9.1).

### Preuve

À une tentative donnée,

\[
\mathbb P(\text{succès avec la branche }i)
=
\frac{c_i}{c_0+c_1}H_i.
\]

La probabilité totale de succès d’une tentative vaut

\[
s
=
\frac{c_0H_0+c_1H_1}{c_0+c_1}>0.
\]

Conditionnellement au premier succès,

\[
\mathbb P(I=i)
=
\frac{\frac{c_i}{c_0+c_1}H_i}{s}
=
\frac{c_iH_i}{c_0H_0+c_1H_1}.
\]

C’est exactement la loi voulue. \(\square\)

La course termine presque sûrement, car chaque tentative réussit avec une probabilité strictement positive.

---

# 10. Algorithme complet d’échantillonnage du random-cluster

## Algorithme `ExactSignedRandomCluster`

### Entrée

Un graphe fini signé et pondéré

\[
G=(V,E,W).
\]

### Préparation

Pour chaque arête :

\[
s_e=\operatorname{sign}(W_e),
\qquad
a_e=|W_e|,
\qquad
p_e=1-e^{-a_e}.
\]

Tirer éventuellement des priorités auxiliaires indépendantes

\[
R_e\sim\operatorname{Exp}(a_e)
\]

et trier les arêtes selon \(R_e\). Tout autre ordre indépendant convient également.

### État

- \(A\leftarrow\varnothing\) ;
- union-find avec parités, initialement composée des singletons.

### Boucle

Pour \(e=\{u,v\}\) dans l’ordre :

1. **Si \(u\) et \(v\) sont dans deux composantes différentes :**
   - construire virtuellement \(S^0\) et \(S^1\) ;
   - appeler `FusionRace` ;
   - si la branche ouverte gagne, ajouter \(e\) à \(A\) et fusionner les deux composantes avec la contrainte \(g_ug_v=s_e\).

2. **Si \(u,v\) sont déjà connectés avec la bonne parité :**
   - tirer \(B_e\sim\operatorname{Bernoulli}(p_e)\) ;
   - ajouter \(e\) à \(A\) si \(B_e=1\).

3. **Si \(u,v\) sont déjà connectés avec la mauvaise parité :**
   - forcer \(e\notin A\).

Retourner \(A\).

---

# 11. Preuve globale d’exactitude

## Théorème 11.1

La sortie de `ExactSignedRandomCluster` suit exactement la loi

\[
\boxed{A\sim\phi_W.}
\tag{11.1}
\]

### Preuve

Fixons l’ordre \(\pi\).

À chaque étape :

- un cycle compatible est ouvert avec sa probabilité conditionnelle exacte \(p_e\) ;
- un cycle frustré est fermé avec probabilité \(1\) ;
- une fusion est ouverte ou fermée selon la probabilité conditionnelle exacte (7.2), grâce aux théorèmes 8.1 et 9.1.

L’algorithme tire donc successivement chaque indicatrice \(\mathbf 1_{\{e_j\in A\}}\) selon sa loi conditionnelle sous \(\phi_W\), sachant toutes les décisions précédentes.

Par la règle de la chaîne,

\[
\mathbb P(A\mid\pi)=\phi_W(A).
\]

Cette identité vaut pour tout ordre \(\pi\). Si l’ordre est lui-même aléatoire et indépendant, sa marginalisation ne change donc pas la loi de \(A\). \(\square\)

---

# 12. Construction du dendrogramme exponentiel exact

L’ordre auxiliaire utilisé par l’algorithme précédent ne doit pas être confondu avec les temps physiques du dendrogramme d’Edwards--Sokal.

Une fois \(A\sim\phi_W\) obtenu, tirer indépendamment, pour chaque arête ouverte,

\[
\boxed{
T_e\mid(e\in A)
\sim
\operatorname{Exp}(a_e)\ \text{conditionnée par }T_e\le1.
}
\tag{12.1}
\]

Sa densité est

\[
\boxed{
f_e(t)
=
\frac{a_ee^{-a_et}}{1-e^{-a_e}}
\mathbf 1_{\{0<t\le1\}}.
}
\tag{12.2}
\]

Par inversion, avec \(U_e\sim\operatorname{Unif}[0,1]\),

\[
\boxed{
T_e
=
-\frac1{a_e}
\log\left[1-U_e(1-e^{-a_e})\right].
}
\tag{12.3}
\]

Pour une arête fermée, poser \(T_e=+\infty\).

La loi jointe exacte des liens et des marques est

\[
\boxed{
\begin{aligned}
\overline\phi_W(A,dT)
\propto{}&
\mathbf 1_{\{A\text{ équilibré}\}}
2^{k(A)}
\prod_{e\in A}
\left[a_ee^{-a_eT_e}\mathbf 1_{\{0<T_e\le1\}}dT_e\right]
\\
&\times
\prod_{e\notin A}e^{-a_e}.
\end{aligned}
}
\tag{12.4}
\]

L’intégration des temps ouverts donne \(1-e^{-a_e}=p_e\), donc la marginale en \(A\) est bien (3.3).

## Dendrogramme de Kruskal

Trier les arêtes ouvertes par \(T_e\) croissant et appliquer Kruskal :

- une arête reliant deux composantes crée un nœud de fusion ;
- une arête fermant un cycle compatible ne change pas la partition.

La forêt de Kruskal \(F\subseteq A\) vérifie, pour tout \(t\in[0,1]\),

\[
\boxed{
\operatorname{CC}\{e\in A:T_e\le t\}
=
\operatorname{CC}\{e\in F:T_e\le t\}.
}
\tag{12.5}
\]

Le dendrogramme signé exact peut donc être stocké sous la forme

\[
D=(F,(T_e)_{e\in F},(s_e)_{e\in F}).
\tag{12.6}
\]

Les arêtes cycliques compatibles sont nécessaires pour représenter le random-cluster complet, mais pas pour représenter sa hiérarchie de composantes.

---

# 13. Échantillonnage exact d’une configuration de Gibbs

Soit \(F\) une forêt couvrante des composantes de \(A\), par exemple la forêt de Kruskal.

Dans chaque composante \(C\) :

1. choisir une racine \(r_C\) ;
2. définir, pour \(v\in C\),
   \[
   g_v
   :=
   \prod_{e\in P_F(r_C,v)}s_e ;
   \tag{13.1}
   \]
3. tirer
   \[
   \xi_C\sim\operatorname{Unif}\{-1,+1\}
   \]
   indépendamment ;
4. poser
   \[
   \boxed{\sigma_v=\xi_Cg_v.}
   \tag{13.2}
   \]

Comme \(A\) est équilibré, toute arête ouverte hors forêt est automatiquement satisfaite par cette jauge : elle forme avec le chemin de forêt un cycle de produit signé \(+1\).

Ainsi, conditionnellement à \(A\), la procédure choisit uniformément l’une des \(2^{k(A)}\) configurations compatibles.

## Théorème 13.1

La configuration obtenue vérifie

\[
\boxed{\sigma\sim\mu_W.}
\tag{13.3}
\]

### Preuve

La procédure reconstruit exactement la conditionnelle \(\Psi_W(\sigma\mid A)\). Comme \(A\sim\phi_W\), le couple \((\sigma,A)\) suit la mesure jointe (2.1). La marginale en spins est \(\mu_W\) par la proposition 2.1. \(\square\)

---

# 14. Pseudocode

## 14.1 Coin de complétion

```text
CompletionCoin(state, first_remaining_index, order):
    D <- copy(state.parity_union_find)

    for position from first_remaining_index to |E|:
        e = order[position]
        (u,v) = endpoints(e)
        p = 1 - exp(-abs(W[e]))
        s = sign(W[e])

        relation = D.relative_parity(u,v)

        if relation == DISCONNECTED:
            if Bernoulli(p) == 1:
                if Bernoulli(1/2) == 0:
                    return 0
                D.union(u,v,constraint=s)

        else if relation == s:
            # Cycle compatible : ouverture/fermeture intégrées,
            # masse totale (1-p)+p = 1.
            continue

        else:
            # Cycle frustré : l'arête doit rester fermée.
            if Bernoulli(1-p) == 0:
                return 0

    return 1
```

## 14.2 Course de fusion

```text
FusionRace(state_closed, state_open, edge, next_index, order):
    p = 1 - exp(-abs(W[edge]))
    proposal_open = p / (2 - p)

    repeat:
        if Bernoulli(proposal_open) == 1:
            I = OPEN
            S = state_open
        else:
            I = CLOSED
            S = state_closed

        if CompletionCoin(S, next_index, order) == 1:
            return I
```

## 14.3 Échantillonneur complet

```text
ExactSignedRandomCluster(G,W):
    for each edge e:
        priority[e] ~ Exp(abs(W[e]))

    order <- edges sorted by increasing priority

    A <- empty set
    D <- parity union-find on V

    for position from 1 to |E|:
        e = order[position]
        (u,v) = endpoints(e)
        p = 1 - exp(-abs(W[e]))
        s = sign(W[e])

        relation = D.relative_parity(u,v)

        if relation == DISCONNECTED:
            S0 <- copy(D)
            S1 <- copy(D)
            S1.union(u,v,constraint=s)

            decision = FusionRace(S0,S1,e,position+1,order)

            if decision == OPEN:
                A.add(e)
                D <- S1

        else if relation == s:
            if Bernoulli(p) == 1:
                A.add(e)

        else:
            # Cycle frustré : fermeture forcée.
            continue

    return A
```

Après cet appel :

1. tirer les temps tronqués (12.3) pour les arêtes de \(A\) ;
2. appliquer Kruskal pour produire le dendrogramme ;
3. colorer les composantes par (13.2) pour obtenir un spin Gibbs exact.

---

# 15. Invariants d’implémentation

Une implémentation doit vérifier systématiquement les invariants suivants.

## 15.1 Équilibre du sous-graphe courant

À tout instant, le sous-graphe effectivement retenu est équilibré.

## 15.2 Convention de parité

Une convention simple est de stocker

\[
\mathrm{parity}[v]=g_vg_{\mathrm{parent}(v)}.
\]

Après compression de chemin, `find(v)` renvoie :

- la racine de \(v\) ;
- le produit \(g_vg_{\mathrm{racine}(v)}\).

Pour deux sommets connectés,

\[
\pi(u,v)
=
(g_ug_r)(g_vg_r)
=
g_ug_v.
\]

## 15.3 Aléas frais

Chaque appel à `CompletionCoin` et chaque tentative de `FusionRace` doivent utiliser de nouveaux aléas indépendants.

## 15.4 Copie de l’état

Les rollouts de complétion travaillent sur des copies de l’union-find. Ils ne doivent jamais modifier l’état principal tant que la course n’a pas accepté une branche.

## 15.5 Séparation des deux familles de temps

- les priorités auxiliaires choisissent l’ordre de décision ;
- les temps tronqués de la section 12 définissent le dendrogramme probabiliste exact.

Les confondre produirait une mauvaise loi jointe des liens et des temps.

---

# 16. Terminaison et complexité

Pour une fusion, la probabilité de succès d’une tentative de course est

\[
\frac{(1-p_e)H_0+\frac{p_e}{2}H_1}
{(1-p_e)+\frac{p_e}{2}}.
\tag{16.1}
\]

Elle est strictement positive puisque \(H_0,H_1>0\). Le nombre de tentatives est donc géométrique et fini presque sûrement.

L’algorithme complet termine presque sûrement sur tout graphe fini.

En revanche, la borne

\[
H_j(S_j)
\ge
\prod_{e\in R_j}(1-p_e)
\tag{16.2}
\]

peut être très petite. Le coût moyen peut donc être exponentiel dans le pire cas.

Cette note établit l’exactitude et la construction directe. Les accélérations possibles — factorisation par composantes, résolution immédiate des parties arborescentes, tables de petits blocs, mémoïsation sur les états de bord, coins hiérarchiques — peuvent être ajoutées ultérieurement sans modifier la preuve, dès lors que le nouvel oracle reste un coin de probabilité exactement égale à \(H_j(S_j)\).

---

# 17. Portée et limites

La construction précédente suppose :

- des spins binaires ;
- un a priori uniforme sans champ externe ;
- une énergie constituée uniquement d’interactions à deux corps signées.

Avec un a priori non uniforme ou un champ externe, le facteur \(2^{k(A)}\) est remplacé par une fonction de partition de composantes. L’algorithme doit alors être modifié.

La construction accepte :

- des poids positifs ou négatifs ;
- une frustration arbitraire ;
- des poids hétérogènes ;
- un graphe non planaire ;
- un graphe non connexe.

---

# 18. Résultat final

À partir du seul graphe observé \((G,W)\), on obtient la chaîne exacte

\[
\boxed{
(G,W)
\longrightarrow
A\sim\phi_W
\longrightarrow
D\text{ dendrogramme exponentiel exact}
\longrightarrow
\sigma\sim\mu_W.
}
\tag{18.1}
\]

La loi du random-cluster est

\[
\boxed{
\phi_W(A)
\propto
\mathbf 1_{\{A\text{ équilibré}\}}
2^{k(A)}
\prod_{e\in A}(1-e^{-|W_e|})
\prod_{e\notin A}e^{-|W_e|}.
}
\tag{18.2}
\]

Le point central de l’échantillonneur est la formule de fusion

\[
\boxed{
\mathbb P(e\text{ ouverte}\mid\text{passé})
=
\frac{\frac{p_e}{2}H_{\mathrm{ouvert}}}
{(1-p_e)H_{\mathrm{fermé}}+\frac{p_e}{2}H_{\mathrm{ouvert}}},
}
\tag{18.3}
\]

sans jamais calculer les masses \(H\) :

- `CompletionCoin` génère exactement des Bernoulli de paramètres \(H\) ;
- `FusionRace` choisit exactement entre les deux branches ;
- la règle de la chaîne donne le random-cluster cible.

---

# Références

1. R. G. Edwards and A. D. Sokal, **Generalization of the Fortuin–Kasteleyn–Swendsen–Wang representation and Monte Carlo algorithm**, *Physical Review D* 38, 2009–2012, 1988. DOI: `10.1103/PhysRevD.38.2009`.
2. S. M. Schmon, A. Doucet and G. Deligiannidis, **Bernoulli Race Particle Filters**, *Proceedings of AISTATS 2019*, PMLR 89, 2350–2358. La course utilisée ici est le cas élémentaire à deux branches ; sa preuve complète est donnée dans la présente note.

