module io
  use accuracy
  implicit none

contains
  
  subroutine readinput(array,bvec,cancel)

    real(dp), allocatable, intent(out) :: array(:,:), bvec(:)
    integer :: dim
    logical, intent(out) :: cancel

    read(21,*) dim
    if (dim <= 0) then
      cancel = .true.
      return
    end if

    allocate(array(dim,dim))
    allocate(bvec(dim))
    
    read(21,*) array
    
    read(21,*) bvec
    array = transpose(array)
    
    
  end subroutine readinput

  subroutine writeoutput(Rarray, res)
    real(dp), intent(in) :: Rarray(:,:), res(:)
    integer :: ii, dim
    character(len=26) :: outputform = '(1000000F10.2)'

    dim = size(Rarray, dim=1)
    write(11,*) "Upper triangle matrix:"
    do ii = 1,dim
      write(11,outputform) Rarray(ii,:)
    end do

    output: do ii = 1,dim
      write(11,"(A,I0,A2,F10.6)") "x_", ii, "=", res(ii)
    end do output
    
  end subroutine writeoutput
  

end module io
