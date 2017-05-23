program gausseli
  use accuracy
  use io
  use eqsolver
  Implicit none

  integer :: ii
  real(dp), allocatable :: res(:), bvec(:), array(:,:), Rarray(:,:), Larray(:,:), Parray(:,:)
  logical :: cancel
  !character(len=26) :: inputfile  = "gauss.inp"
  !character(len=26) :: outputfile = "output.dat"

  !Inputfile:
  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  !Outputfile:
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    call readinput(array, bvec, cancel)
    if (cancel) then
      exit mainloop
    end if
      
    call ludecompose(array, Rarray, Larray, Parray)

    call substituteback(array, Rarray, Larray, Parray, bvec, res)

    call writeoutput(Rarray, res)

    deallocate(array)
    deallocate(Rarray)
    deallocate(Larray)
    deallocate(Parray)
    deallocate(bvec)
    
  end do mainloop
  close(11)
  close(21)

end program gausseli
