# To add a new set of functions --
#   Add any ".o" files that it uses to the OBJ list
#   Add a line for each function to the FUN list
#   Add a line for each documented function (should be all of them!) to DOC
# If the group is large, you may want to separate it out just to help keep
#   this file more readable, as is done with "SURV" below.
#
#   Use  "make" or "make all" to make a local copy of S.  You may want to
#     run S in this directory just to make sure all of your changes work.
#   make fun - functions only
#   make doc - documentation only
#

OBJ   = bdsmatrix_index1.o \
	bdsmatrix_index2.o \
	bdsmatrix_index3.o \
	bdsmatrix_prod.o  \
	bdsmatrix_prod2.o  \
	bdsmatrix_prod3.o  \
	bdsmatrix_prod4.o  \
	chinv4.o \
	chinv5.o \
	cholesky4.o \
	cholesky5.o \
	chsolve4.o \
	chsolve5.o \
	dmatrix.o \
	gcback.o \
	gchol.o \
	gchol_bds.o

HFILE=  bdsS.h \
	bdsmatrix.h 

FUN =	bdsmatrix \
	gchol \
	as.matrix.bdsmatrix \
	backsolve \
	bdsBlock \
	bdsI \
	bdsmatrix.ibd \
	bdsmatrix.reconcile \
	diag.bdsmatrix \
	gchol.bdsmatrix \
	listbdsmatrix \
	multiply.bdsmatrix \
	solve.bdsmatrix \
	solve.gchol \
	solve.gchol.bdsmatrix

NOWEB= 

# Define command and locations
Scmd    =  R --save --quiet --no-site-file --no-init-file

Ssrc    = ../R
Sccode  = ../src
NWsrc   = ../noweb

# Nothing from here on down should ever have to change
all:  fun bdsmatrix.so

# These are the implicit rules that drive the makefile.  If you can
#   read them, consider yourself an expert!

touchdir/%.s: $(Ssrc)/%.R
	$(Scmd) -f $< > /dev/null
	touch $@

touchdir/%.R: $(NWsrc)/all.nw
	notangle -R$* $(NWsrc)/all.nw > $@
	$(Scmd) -f $@ > /dev/null

$(NWsrc)/all.nw: $(NWsrc)
	$(MAKE) -C $(NWsrc) all.nw

%.o: %.c
	R CMD COMPILE $<

#%.h: $(Sccode)/%.h
#	cp $< $@

LFUN = $(FUN:%=touchdir/%.s)
NFUN = $(NOWEB:%=touchdir/%.R)
SrcOBJ= $(OBJ:%=../src/%)

fun: $(LFUN) $(NFUN)


bdsmatrix.so:  $(SrcOBJ) 
	R CMD SHLIB -o bdsmatrix.so $(SrcOBJ) 

