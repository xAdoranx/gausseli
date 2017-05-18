This program will solve a system like Ax=b

The files "accuracy.f90", "io.f90", "eqsolver.f90", "gausseli.f90"
must be compiled in this order.


Required input file named "gauss.inp".
Formatted like

3           !Dimension of the Matrix or number of variables
1 2 3       !
1 2 3       ! rows of the Matrix A
1 2 3       !
2 3 4       ! solution vector b
5
3 2 1 5 3
.
.
.
0           ! with dimension <= 0 the programm will be stopped

Solution is given in file "output.dat"

Error: "underterminated system" will occur, if the input matrix is underterminated
       the program will continue with the next matrices.




