# Community detection on signed graphs

[![Research checks](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/actions/workflows/research-checks.yml/badge.svg)](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/actions/workflows/research-checks.yml)

Ce dépôt réunit trois présentations, le chapitre de thèse qui les motive et
un cahier de recherche reproductible sur la weak recovery dans les graphes
signés.

> [!IMPORTANT]
> **Résultat rigoureux actuel.** Sur le tore triangulaire, le dépôt établit
> l'absence de weak recovery pour tout
> $`p\in[1/2,0.809439]`$, soit
> $`p_{\mathrm{WR}}\ge0.809439`$. La
> [preuve canonique](research/hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md)
> utilise un canal triangulaire multi-état. Elle est rigoureuse, mais elle
> **ne provient pas** encore de la dynamique hiérarchique.

## Par où commencer ?

Choisissez le parcours qui correspond à votre objectif.

| Je veux… | Première lecture | Puis… |
|---|---|---|
| comprendre le résultat en cinq minutes | [statut scientifique actuel](research/hierarchical-swendsen-wang/CURRENT_STATUS.md) | [certificat à p = 0,809439](research/hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) |
| comprendre la dynamique hiérarchique | [cible prioritaire répliquée](research/hierarchical-swendsen-wang/active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) | [pilote SBM et port global fini](research/hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md) |
| retrouver une note précise | [index exhaustif des notes](research/hierarchical-swendsen-wang/INDEX.md) | dossier indiqué par son statut |
| reproduire les preuves assistées par calcul | [guide des calculs](research/hierarchical-swendsen-wang/computations/README.md) | commandes de validation ci-dessous |
| consulter les supports de séminaire | [table des présentations](#présentations) | source et PDF dans chaque dossier |

## Les deux voies scientifiques

Le dépôt sépare désormais explicitement deux axes qui ne doivent pas être
confondus.

| voie | acquis | question ouverte |
|---|---|---|
| canal triangulaire multi-état | borne rigoureuse $`p_{\mathrm{WR}}\ge0.809439`$ | pousser le certificat vers le point tangent candidat |
| dynamique hiérarchique | mesure jointe exacte, Gibbs d'arbre entier et calibration broadcast | annuler le reste signé inter-cellules de deux hiérarchies |

La cible triangulaire utilise chaque arbre associé à une composante géante,
coupé au niveau $`\beta_c(p)=q_p^{-1}(q_c)`$, mais conserve le Gibbs exact
de l'arbre entier. Pour calculer le carré postérieur, il faut deux
dendrogrammes indépendants conditionnellement à l'observation : partager un
dendrogramme gonfle artificiellement l'overlap. Le pilote dérive le jacobien
$`d\theta^2`$ et retrouve le seuil du broadcast après marginalisation
exacte, pour toute coupe ; il ne prouve ni une contraction dynamique ni le
transfert au SBM fini, où balance ou non-arêtes recouplent les racines. Sur
le graphe fini, ce recouplage est désormais écrit exactement comme un port
scalaire de magnétisation, sans que sa comparaison au broadcast soit encore
prouvée. Sur la grille, les contributions hors double géante et la diagonale
critique sont désormais éliminées. La cible exacte est le reste signé entre
cellules distinctes de l'intersection des deux géantes. Son statut est
**programme de preuve**, pas théorème. Ce programme ne répète pas la borne
par recoloriage du chapitre 11 : il remplace son unique objet gelé par deux
Gibbs conditionnels exacts sur deux dendrogrammes complets.

Le premier audit à $`L=4,p=0.81`$ trouve l'enveloppe à un dendrogramme
macroscopique ($`0.9507\ldots`$) et le reste signé à deux dendrogrammes encore
positif en moyenne. Il valide la décomposition, pas une tendance en volume.

## Carte du dépôt

| chemin | rôle |
|---|---|
| [`research/`](research/) | porte d'entrée de la recherche |
| [`research/hierarchical-swendsen-wang/`](research/hierarchical-swendsen-wang/) | projet principal, organisé par statut scientifique |
| [`beamer-presentation/`](beamer-presentation/) | séminaire MathNet du 15 juin 2026 |
| [`beamer-presentation-neo/`](beamer-presentation-neo/) | séminaire NEO du 25 juin 2026 |
| [`beamer-presentation-reunion-2026-07-16/`](beamer-presentation-reunion-2026-07-16/) | réunion de recherche du 16 juillet 2026 |
| [`ChapII.tex`](ChapII.tex) | source du chapitre 11 utilisé comme point de départ |
| [`Manuscrit_de_these.pdf`](Manuscrit_de_thèse.pdf) | manuscrit de thèse complet |

Dans le projet de recherche, les dossiers ont un sens précis :

```text
hierarchical-swendsen-wang/
├── CURRENT_STATUS.md   état de l'art interne et prochaine étape
├── INDEX.md            catalogue complet des notes
├── foundations/        identités et cadres réutilisables
├── results/            théorèmes et certificats établis
├── active/             programme hiérarchique prioritaire et outils associés
├── diagnostics/        expériences finies, benchmarks et no-go
├── archive/            anciennes feuilles de route et jalons subsumés
├── computations/       scripts et tests reproductibles
└── references/         littérature et bibliographie
```

Les numéros `00` à `40` gardent la chronologie du cahier. Ils ne définissent
plus l'ordre de lecture.

## Présentations

| date | support | sources | PDF |
|---|---|---|---|
| 15 juin 2026 | Séminaire MathNet | [`beamer-presentation/`](beamer-presentation/) | [`main.pdf`](beamer-presentation/main.pdf) |
| 25 juin 2026 | Séminaire NEO | [`beamer-presentation-neo/`](beamer-presentation-neo/) | [`PresentationNIM_2026-06-25.pdf`](beamer-presentation-neo/PresentationNIM_2026-06-25_LouisHauseux_ABayesianFrameworkForCommunityDetectionOnSignedGraphs.pdf) |
| 16 juillet 2026 | Réunion de recherche | [`beamer-presentation-reunion-2026-07-16/`](beamer-presentation-reunion-2026-07-16/) | [`Presentation_2026-07-16.pdf`](beamer-presentation-reunion-2026-07-16/Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf) |

Chaque présentation est un instantané autonome : son thème, ses images et sa
bibliographie restent dans son propre dossier. Pour la compiler :

```bash
cd beamer-presentation
make
```

## Reproduire et valider

Depuis la racine du dépôt :

```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_*.py'
python3 -m compileall -q \
  research/hierarchical-swendsen-wang/computations
```

Les calculs n'ont pas de dépendance scientifique externe. Une valeur
flottante n'est traitée comme une preuve que lorsqu'elle est accompagnée
d'un certificat exact, rationnel ou par intervalles.

## Lire les statuts sans ambiguïté

- **établi** : preuve complète dans les hypothèses annoncées ;
- **conditionnel** : implication prouvée sous un lemme explicitement nommé ;
- **diagnostic** : calcul fini ou simulation, sans extrapolation en volume ;
- **no-go** : raccourci réfuté ou fermeture démontrée inadéquate ;
- **actif** : objectif de recherche non encore démontré ;
- **archivé** : document conservé pour l'historique, mais non directeur.

Quelques fichiers historiques restent volontairement à la racine afin de
préserver les liens existants : les figures `LargestComp*` et `Overlap_*`, le
gabarit Beamer Inria, le PDF Sankararaman–Baccelli, la bibliographie de thèse
et les sources du chapitre. Les présentations n'en dépendent pas : elles
contiennent leurs propres copies compilables.
