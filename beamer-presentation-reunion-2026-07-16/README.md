# Réunion de recherche — 16 juillet 2026

Support de la réunion Louis Hauseux, Nahuel Soprano-Loto et Konstantin
Avrachenkov consacrée à la dynamique hiérarchique de Swendsen--Wang.

> [!IMPORTANT]
> Les slides décrivant la dynamique, les coupes de fusion et la chaîne des
> $`\Lambda_v`$ sont regroupées dans
> [`hierarchical_sw_frames.tex`](hierarchical_sw_frames.tex). Leur programme
> mathématique actualisé se trouve dans le
> [dossier de recherche](../research/).

- [PDF compilé](Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf)
- [Source principale](main.tex)
- [Frames hiérarchiques](hierarchical_sw_frames.tex)
- [Bibliographie](referencesThesis.bib)

## Compilation

```bash
make
```

Pour supprimer les auxiliaires :

```bash
make clean
```

Le thème Beamer Inria est inclus dans [`theme/`](theme/).

## Parcours conseillé

1. lire les slides 31–33 dans le PDF ;
2. consulter l'[audit des affirmations](../research/02_AUDIT.md) ;
3. poursuivre avec la [construction corrigée](../research/01_HIERARCHIE.md)
   et la [question de recovery](../research/03_RECOVERY.md).
