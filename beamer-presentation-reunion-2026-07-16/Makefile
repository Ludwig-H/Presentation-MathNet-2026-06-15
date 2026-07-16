PDFNAME = Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf

all: $(PDFNAME)

$(PDFNAME): main.tex hierarchical_sw_frames.tex
	latexmk -jobname=Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin main.tex

clean:
	latexmk -jobname=Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin -C
	@rm -f Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.nav Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.snm
