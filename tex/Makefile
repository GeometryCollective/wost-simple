FILE = WoSt-tutorial

all:
	pdflatex -shell-escape $(FILE)
	open $(FILE).pdf

bib:
	pdflatex -shell-escape $(FILE)
	bibtex $(FILE)
	pdflatex -shell-escape $(FILE)
	open $(FILE).pdf

im:
	roger trios media/*.trio.json --tex-labels -o images/

clean:
	rm -f $(FILE).pdf
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -rf _minted*
	rm -rf svg-inkscape

