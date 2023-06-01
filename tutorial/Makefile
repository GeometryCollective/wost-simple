FILE = WoSt-tutorial

all:
	pdflatex -shell-escape $(FILE)
	bibtex $(FILE)
	pdflatex -shell-escape $(FILE)
	open $(FILE).pdf

im:
	roger trio media/WalkOnStars.trio.json -o images/WalkOnStars.svg
	roger trio media/AtomicClosestPoint.trio.json -o images/AtomicClosestPoint.svg
	roger trio media/AtomicClosestRay.trio.json -o images/AtomicClosestRay.svg

clean:
	rm -f $(FILE).pdf
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -rf _minted*
	rm -rf svg-inkscape

