IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -L$(LDIR)

ODIR=obj
LDIR =../lib

LIBS=-llevmar -llapack -lblas -lm

_DEPS = oofit.h levmar.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ooread.o fit.o main.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

oofit: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

