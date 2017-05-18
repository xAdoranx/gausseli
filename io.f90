module io
  use accuracy
  implicit none

contains
  
  subroutine readinput(dim,array,bvec,cancel)

    real(dp), allocatable, intent(out) :: array(:,:), bvec(:)
    integer, intent(out) :: dim
    logical, intent(out) :: cancel

    read(21,*) dim
    if (dim <= 0) then
      cancel = .true.
      return
    end if

    allocate(array(4*dim,dim))
    allocate(bvec(dim))
    
    read(21,*) array(1:dim,:)
    
    read(21,*) bvec
    array(1:dim,:) = transpose(array(1:dim,:))
    
    
  end subroutine readinput

  subroutine writeoutput(dim, array, res)
    real(dp), allocatable, intent(inout) :: array(:,:), res(:)
    integer, intent(in) :: dim
    integer :: ii
    character(len=26) :: outputform = '(1000000F10.2)'

    write(11,*) "Upper triangle matrix:"
    do ii = 1,dim
      write(11,outputform) array(dim+ii,:)
    end do

    output: do ii = 1,dim
      write(11,"(A,I0,A2,F10.6)") "x_", ii, "=", res(ii)
    end do output

    deallocate(array)
    deallocate(res)
    
  end subroutine writeoutput
  

end module io
