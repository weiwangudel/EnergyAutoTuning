ALL = SEAK.pdf
PS  = SEAK.ps
PDF = SEAK.pdf

TEX = \
	SEAK.tex \
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
	SEAK.bib 

all:	$(ALL)


SEAK.pdf: $(TEX) $(FIG) $(BOILERPLATE) Makefile $(BIBS)
	yes | tr 'y' ' ' | pdflatex SEAK
	yes | tr 'y' ' ' | bibtex SEAK
	yes | tr 'y' ' ' | pdflatex SEAK
	yes | tr 'y' ' ' | pdflatex SEAK
	yes | tr 'y' ' ' | evince SEAK.pdf


# we should have dependence rules for the bibs too.

.PHONY:	ps pdf 

pdf:	$(PDF)

clean:
	yes | rm -f $(PS) $(PDF) *.dvi *.log *.aux *.bbl *.blg SEAK.out *~
