###--- Makefile (adapted from pcalg package)

DRIVERS = lmrob_simulation

## "Writing R Extensions" requires: 1st target must also clean up
doc: all
	${MAKE} clean

all: $(DRIVERS:=.tex) $(DRIVERS:=.pdf)

## render pdf and compress it using qpdf
%.pdf: %.tex *.bib
	texi2pdf $<
	qpdf --stream-data=compress $@ temp$@
	mv temp$@ $@

purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl) $(DRIVERS:=-*.pdf)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
