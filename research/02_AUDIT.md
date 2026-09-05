# Audit : trois corrections décisives

## Une feuille conditionnée n'est pas Glauber

Deux sommets, une arête de poids $`W>0`$. Si le dendrogramme complet contient
leur fusion, les deux spins doivent être égaux. Retourner une seule feuille
a donc une probabilité nulle. Le Glauber de la postérieure autorise ce
retournement avec probabilité $`1/(1+e^W)>0`$.

**Correction :** pour retrouver Glauber à la coupe zéro, oublier les fusions
supérieures. À une coupe intermédiaire, cela restitue les interactions
résiduelles de la [construction](01_HIERARCHIE.md).
La formule conditionnelle du dendrogramme complet dans la soutenance est
correcte ; c'est son identification aux extrêmes qu'il faut préciser.

## Une coupe ne rend pas ses amas indépendants

L'argument du chapitre 11 utilise le recoloriage indépendant et uniforme
des clusters. L'invariance de la dynamique ne suffit pas à cet argument.
À $`t<1`$, les orientations ont des interactions $`(1-t)W`$ entre amas.
Des amas petits peuvent donc encore transmettre collectivement de l'information.

**Correction :** conserver ces interactions dans une preuve de recovery.
Une fusion n'est pas un canal indépendant de ses ancêtres ; multiplier
des fiabilités locales exige une preuve supplémentaire.
Deux répliques partageant le même dendrogramme donnent en général une
majoration par convexité, pas le carré postérieur exact.

## Deux détails de preuve dans le chapitre

Pour les petits clusters de taille inférieure à $`\delta n`$,

```math
\sum_{C:\,|C|<\delta n}\frac{|C|^2}{n^2}\le\delta.
```

La fluctuation est donc $`O_{\mathbb P}(\sqrt\delta)`$, pas automatiquement
$`o_{\mathbb P}(1)`$ à $`\delta`$ fixé. La preuve se répare en prenant d'abord
la limite supérieure en $`n`$, puis $`\delta\downarrow0`$.

Pour un prior général, « préserver le prior » ne prouve pas la balance
détaillée annoncée dans le résultat sur les recoloriages. Il faut écrire
une conditionnelle exacte ou une balance explicite. Le présent dossier
se place sous prior i.i.d. uniforme, où ces difficultés disparaissent.

## Sources et portée

Lecture du [chapitre 11, source actuelle](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/main/Manuscrit_de_these/Manuscrit%20these%20Louis%20Hauseux/PartIII/ChapII.tex),
notamment « Dynamiques de clusters » et « Borne sur la fraction recouvrable via percolation » ;
du [chapitre local](../ChapII.tex) ; des présentations
[MathNet](../beamer-presentation/main.tex), [NEO](../beamer-presentation-neo/main.tex),
de [juillet](../beamer-presentation-reunion-2026-07-16/hierarchical_sw_frames.tex)
et de la [soutenance](https://github.com/Ludwig-H/Manuscrit-de-th-se/blob/main/Soutenance/soutenance/Soutenance_These_2026-09-08_LouisHauseux.pdf),
page PDF 62 et compléments 93–98 sur les horloges.
Consultation : 5 septembre 2026. Les supports eux-mêmes ne sont pas modifiés.

Le certificat antérieur à 0,809439 a été contrôlé séparément : calcul rationnel
reproduit et hypothèses des théorèmes utilisés vérifiées, sans défaut identifié.
Il relève d'un canal triangulaire, pas de l'interpolation proposée ici.
Les anciennes notes sont remplacées par ces fondations courtes ; leur
[version antérieure](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/tree/b89adcc07bce6b13d7732233c183a0eb63654d99/research) reste dans l'historique Git.
