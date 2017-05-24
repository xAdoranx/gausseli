module io
  use accuracy
  implicit none

contains
  
  subroutine readinput(array,bvec,cancel,outform)

    real(dp), allocatable, intent(out) :: array(:,:), bvec(:,:)
    integer :: dim, bdim
    logical, intent(out) :: cancel
    character(len=8), intent(out) :: outform

    read(21,*) dim
    if (dim <= 0) then
      cancel = .true.
      return
    end if

    allocate(array(dim,dim))
    
    read(21,*) array

    read(21,*) bdim

    write(*,*) bdim

    allocate(bvec(dim,bdim))
    read(21,*) bvec
    array = transpose(array)

    read(21,*) outform
    
  end subroutine readinput


  subroutine writetofile(Rarray, res, outform)
    real(dp), intent(in) :: Rarray(:,:), res(:,:)
    character(len=8), intent(in) :: outform
    integer :: ii, dim
    character(len=26) :: outputform = '(1000000F10.2)'

    dim = size(Rarray, 1)
    
    select case (outform)
    case ("simpfile")
      do ii = 1,dim
        write(11,"(100ES23.15)") res(ii,:)
      end do
    case ("compfile")
      write(11,*) "Upper triangle matrix:"
      do ii = 1,dim
        write(11,outputform) Rarray(ii,:)
      end do
      do ii =1,dim
        write(11,"(A,I0,A2,100ES23.15)") "x_", ii, "=", res(ii,:)
      end do
    end select
    
    
  end subroutine writetofile
  

  subroutine writetoscreen(Rarray, res, outform)

    real(dp), intent(in) :: Rarray(:,:), res(:,:)
    character(len=8), intent(in) :: outform
    integer :: ii, dim
    character(len=26) :: outputform = '(1000000F10.2)'

    dim = size(Rarray, 1)
    select case (outform)
    case ("simpscrn")
      do ii = 1,dim
        write(*,"(100ES23.15)") res(ii,:)
      end do
    case ("compscrn")
      write(*,*) "Upper triangle matrix:"
      do ii = 1,dim
        write(*,outputform) Rarray(ii,:)
      end do
      do ii =1,dim
        write(*,"(A,I0,A2,100ES23.15)") "x_", ii, "=", res(ii,:)
      end do
    end select
    

  end subroutine writetoscreen
  

end module io
