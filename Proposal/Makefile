all: proposal

proposal:
	pdflatex proposal
	bibtex proposal
	pdflatex proposal
	pdflatex proposal


clean:
	rm *~ *.blg *.log *.dvi *.aux *.brf *.out *.bbl $(INCLUDED)
