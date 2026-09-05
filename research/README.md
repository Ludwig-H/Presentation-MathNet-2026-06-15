# Une dynamique de clusters vraiment hiérarchique

**L'idée des horloges est bonne. Le point à corriger est le conditionnement.**
Couper un dendrogramme complet tout en conservant ses fusions supérieures ne
redonne pas Glauber aux feuilles. Il faut intégrer les informations situées
au-dessus de la coupe et garder les interactions restantes entre amas.

La construction retenue conserve exactement la postérieure du chapitre 11 :

| Coupe | Mise à jour |
|---|---|
| zéro : sommets isolés | balayage de Glauber |
| intermédiaire | retournements d'amas avec interactions résiduelles |
| un : racines Swendsen–Wang | retournements indépendants et uniformes |

Trois notes, dans l'ordre :

1. [La hiérarchie et ses lois](01_HIERARCHIE.md) : construction, probabilités, preuve d'invariance.
2. [L'audit](02_AUDIT.md) : ce qui tient et ce qui doit être corrigé.
3. [La coupe critique et la recovery](03_RECOVERY.md) : la quantité précise qu'il reste à contrôler.

**Acquis :** une interpolation exacte en volume fini.
**Ouvert :** une meilleure borne de weak recovery grâce à cette dynamique.
La borne antérieure à 0,809439 reste accessible par son
[certificat historique](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/blob/b89adcc07bce6b13d7732233c183a0eb63654d99/research/hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) ; elle ne vient pas de la dynamique.

Vérification sur de petits graphes, depuis la racine du dépôt :

```bash
python3 research/check_hierarchy.py
```

Le script énumère les états et les transitions ; ses contrôles numériques
complètent les preuves, sans établir de seuil asymptotique.
