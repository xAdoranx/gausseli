program gausseli
  use accuracy
  use io
  use eqsolver
  Implicit none

  integer :: dim, ii
  real(dp), allocatable :: res(:), bvec(:), array(:,:)
  logical :: cancel
  !character(len=26) :: outputform = '(1000000F10.2)'

  
  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    call readinput(dim, array, bvec, cancel)
    if (cancel) then
      exit mainloop
    end if
      
    call ludecompose(dim, array)

    !do ii=1,4*dim
    !  write(*,outputform) array(ii,:)
    !end do    

    call substituteback(dim, array, bvec, res)

    call writeoutput(dim, array, res)
        
  end do mainloop
  close(11)
  close(21)

end program gausseli
