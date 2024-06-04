RM = rm -f
CP = cp
TOUCH = touch
REXE = R --vanilla
RCMD = $(REXE) CMD
RSCRIPT = Rscript --no-save --no-restore

CC = gcc
CXX = g++
CFLAGS = -O4 -Wall
CXXFLAGS = -O4 -Wall
LDLIBS =
FC = gfortran
FFLAGS = -O4
INSTALL = install

LATEX = latex -shell-escape
PDFLATEX = pdflatex -shell-escape
BIBTEX = bibtex
MAKEIDX = makeindex

DVIPS = dvips -tletter -Ppdf
PSTOPDF = ps2pdf
EPSTOPDF = epstopdf
TTH = htlatex # tth -E2 -e2
TTM = ttm -E2 -e2

default: glmer.pdf ## glmerAppendix.pdf
## in other dir:  lmer.pdf

.PHONY:

%.tex: %.Rnw
	$(RSCRIPT) -e "library(knitr); knit(\"$*.Rnw\")"
##	$(RSCRIPT) -e "Sweave(\"$*.Rnw\")"

%.tex: %.txt
	pandoc -s -o $*.tex $*.txt

%.R: %.Rnw
	$(RCMD) Stangle $*

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.ps
	$(PSTOPDF) $*.ps $@

%.ps: %.dvi
	$(DVIPS) $*.dvi -o $@

%.html: %.tex
	$(TTH) $*

%.xml: %.tex
	$(TTM) $*

%.pdf: %.fig
	$(F2D) -Lpdf $*.fig $@

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

%.eps: %.fig
	$(F2D) -Leps $*.fig $@

%.png: %.fig
	$(F2D) -Lpng $*.fig $@

%.bbl: %.tex
	-$(LATEX) $*
	$(BIBTEX) $*

%.idx: %.tex
	-$(LATEX) $*

%.ind: %.idx
	$(MAKEIDX) $*

%: %.cc
	$(CXX) $(CXXFLAGS) $*.cc $(LDLIBS) -o $*

%: %.c
	$(CC) $(CFLAGS) $*.c $(LDLIBS) -o $*

%: %.f
	$(FC) $(FFLAGS) $*.f $(LDLIBS) -o $*

clean:
	$(RM) *.o *.log *.aux *.toc *.blg *.bbl *.out *-???.pdf Rplots.ps

.SECONDARY:
