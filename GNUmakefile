
FC = gfortran
FCOPTS =

LN = $(FC)
LNOPTS = 

NAME = gausseli

FILAC = accuracy.f90
FILIO = io.f90
FILEQ = eqsolver.f90
FILGS = gausseli.f90

OBJAC = accuracy.o
OBJIO = io.o
OBJEQ = eqsolver.o
OBJGS = gausseli.o
ALOBJS = $(OBJGS) $(OBJAC) $(OBJIO) $(OBJEQ)


$(NAME): $(ALOBJS)
	$(LN) $(LNOPTS) -o $@ $^

$(OBJAC): $(FILAC)
	$(FC) $(FCOPTS) -c $<

$(OBJIO): $(FILIO) $(OBJAC)
	$(FC) $(FCOPTS) -c $<

$(OBJEQ): $(FILEQ) $(OBJAC)
	$(FC) $(FCOPTS) -c $<

$(OBJGS): $(FILGS) $(OBJAC) $(OBJIO) $(OBJEQ)
	$(FC) $(FCOPTS) -c $<

.PHONY: clean realclean

clean:
	rm -f *.o *.mod

realclean: clean
	rm -f gausseli
