FILE = WoSt-tutorial

all:
	pdflatex -shell-escape $(FILE)
	bibtex $(FILE)
	pdflatex -shell-escape $(FILE)
	open $(FILE).pdf

im:
	roger trio media/AtomicClosestPoint.trio.json --tex-labels -o images/AtomicClosestPoint.svg
	roger trio media/AtomicClosestRay.trio.json --tex-labels -o images/AtomicClosestRay.svg
	roger trio media/AtomicClosestSilhouette.trio.json --tex-labels -o images/AtomicClosestSilhouette.svg
	roger trio media/ClosestPoint.trio.json --tex-labels -o images/ClosestPoint.svg
	roger trio media/FirstRay.trio.json --tex-labels -o images/FirstRay.svg
	roger trio media/WalkOnSpheres.trio.json --tex-labels -o images/WalkOnSpheres.svg
	roger trio media/WalkOnStars.trio.json --tex-labels -o images/WalkOnStars.svg
	roger trio media/WalkOnStarsTrajectories.trio.json --tex-labels -o images/WalkOnStarsTrajectories.svg

clean:
	rm -f $(FILE).pdf
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -rf _minted*
	rm -rf svg-inkscape

