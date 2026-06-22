all: main.pdf

main.pdf: main.tex
	latexmk main.tex

clean:
	latexmk -C
	@rm -f main.nav main.snm

