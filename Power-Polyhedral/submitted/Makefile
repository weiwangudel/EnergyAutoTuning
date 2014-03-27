ALL = EAT.pdf
PS  = EAT.ps
PDF = EAT.pdf

TEX = \
	EAT.tex \
	abstract.tex \
	benchmarks.tex \
	conclusion.tex \
	introduction.tex \
	related.tex \
	results.tex \
	setup.tex \
	tools.tex 

FIG = \


BIBS = \
	EAT.bib 

all:	$(ALL)


EAT.pdf: $(TEX) $(FIG) $(BOILERPLATE) Makefile $(BIBS)
	yes | tr 'y' ' ' | pdflatex EAT
	yes | tr 'y' ' ' | bibtex EAT
	yes | tr 'y' ' ' | pdflatex EAT
	yes | tr 'y' ' ' | pdflatex EAT
	yes | tr 'y' ' ' | evince EAT.pdf


# we should have dependence rules for the bibs too.

.PHONY:	ps pdf 

pdf:	$(PDF)

clean:
	yes | rm -f $(PS) $(PDF) *.dvi *.log *.aux *.bbl *.blg EAT.out *~
