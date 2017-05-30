!> program to solve a given system like Ax=b.
!!
!! Requires modules accuracy, io, eqsolver to work.
program gausseli
  use accuracy
  use io
  use eqsolver
  Implicit none

  integer :: ii
  real(dp), allocatable :: res(:,:), bvec(:,:), array(:,:), RLarray(:,:), Parray(:,:)
  logical :: cancel
  character(len=8) :: outform
  !character(len=26) :: inputfile  = "gauss.inp"
  !character(len=26) :: outputfile = "output.dat"

  !> Inputfile is defined with the ID 21
  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  !> Outputfile is defined with the ID 11
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    call readinput(array, bvec, cancel, outform)
    if (cancel) then
      exit mainloop
    end if
      
    call ludecompose(array, RLarray, Parray)

    call substituteback(RLarray, Parray, bvec, res)

    select case (outform)
    case ("simpfile", "compfile")
      call writetofile(RLarray, res, outform)
    case ("simpscrn", "compscrn")
      call writetoscreen(RLarray, res, outform)
    end select
      
    

    deallocate(array)
    deallocate(RLarray)
    deallocate(Parray)
    deallocate(bvec)
    
  end do mainloop
  close(11)
  close(21)

end program gausseli
