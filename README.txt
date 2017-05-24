This program will solve a system like Ax=b

To compile the program use GNUmakefile by typing "make"
For cleaning unneccessary files type "make clean"
To clean the whole program type "make realclean"


Required input file named "gauss.inp".
Formatted like

3           !Dimension of the Matrix or number of variables
1 2 3       !
1 2 3       ! rows of the Matrix A
1 2 3       !
2           ! amount of solution vectors
2 3 4       ! solution vector b 1
3 4 5       ! solution vector b 2
simpscrn    ! Defines output options
5
3 2 1 5 3
.
.
.
.
0           ! If dimension <= 0 the programm will be stopped

Outputoptions:
      simpfile: only solutionvector on screen
      compfile: full output on screen
      simpscrn: only solutionvektor in file
      compscrn: full output in file

Error: "underterminated system" will occur, if the input matrix is underterminated
       the program will continue with the next matrices.




