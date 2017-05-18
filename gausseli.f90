program gausseli
  use accuracy
  use io
  use eqsolver
  Implicit none

  integer :: dim
  real(dp), allocatable :: Rarray(:,:), Larray(:,:), Parray(:,:), res(:), bvec(:),&
      & array(:,:)
  logical :: cancel

  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    call readinput(dim, array, bvec, cancel)
    if (cancel) then
      exit mainloop
    end if
      
    call ludecompose(array, Rarray, Larray, Parray)

    deallocate(array)

    call substituteback(Rarray, Larray, Parray, bvec, res)

    deallocate(Parray)
    deallocate(Larray)
    deallocate(bvec)

    call writeoutput(Rarray, res)
        
  end do mainloop
  close(11)
  close(21)

end program gausseli
