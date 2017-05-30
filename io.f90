!> Contains Input/Output related routines
module io
  use accuracy
  implicit none

contains

  !> Reads in a Matrix to be solved together with several solution vectors.
  subroutine readinput(array,bvec,cancel,outform)

    !> Input Matrix
    !!
    !! Must be read in from ID: 21.
    !! Needs first the Dimension of the matrix and then the matrix row by row.
    real(dp), allocatable, intent(out) :: array(:,:)
    !> Contains solution vectors
    !!
    !! Must be read in from ID: 21.
    !! Needs first the amount of solution vectors and then all solution vectors.
    real(dp), allocatable, intent(out) :: bvec(:,:)
    integer :: dim
    integer :: bdim
    !> Abort condition
    !!
    !! If amount of variables <= 0, cancel = .true.
    logical, intent(out) :: cancel
    !> Defines output style
    !!
    !! Must be one of the following characters
    !! "simpfile": Writes only solution vectors into file
    !! "compfile": Writes full output into file
    !! "simpscrn": Writes only solution vectors into file
    !! "compscrn": Writes full output on screen
    character, intent(out) :: outform

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

  !> Output in file

  subroutine writetofile(RLarray, res, outform)
    !> Reads in Matrix
    !!
    !! Must contain upper and under triangle Matrix of the input Matrix.
    real(dp), intent(in) :: RLarray(:,:), res(:,:)
    !> Defines output style
    !!
    !! Must be one of the two following characters
    !! "simpfile": writes only solution vectors into file
    !! "compfile": writes full output into file
    character(len=8), intent(in) :: outform
    integer :: ii, dim
    character(len=26) :: outputform = '(1000000F10.2)'
    real(dp), allocatable :: Rarray(:,:)

    dim = size(RLarray, 1)

    allocate(Rarray(dim,dim))
    do ii = 1,dim
      Rarray(ii,ii:dim) = RLarray(ii,ii:dim)
    end do
        
    select case (outform)
    case ("simpfile")
      do ii = 1,dim
        write(11,"(100ES23.15)") res(ii,:)
      end do

    case ("compfile")
      write(11,*) "upper triangle matrix:"
      do ii = 1,dim
        write(11,outputform) Rarray(ii,:)
      end do
      do ii =1,dim
        write(11,"(A,I0,A2,100ES23.15)") "x_", ii, "=", res(ii,:)
      end do
    end select
    
    
  end subroutine writetofile

  !> Output on Screen

  subroutine writetoscreen(RLarray, res, outform)

    !> Read in Matrix
    !!
    !! Must contain upper and under triangle Matrix of the input matrix.
    real(dp), intent(in) :: RLarray(:,:)
    !> Solution vector
    !!
    !! Must contain the solutions for all variables and all input vectors.
    real(dp), intent(in) :: res(:,:)
    !> Defines output style
    !!
    !! Must be one of the following two characters:
    !! "simpscrn": write only solution vectors on screen
    !! "compscrn": write full output on screen
    character(len=8), intent(in) :: outform
    integer :: ii, dim
    character(len=26) :: outputform = '(1000000F10.2)'
    real(dp), allocatable :: Rarray(:,:)

    dim = size(RLarray, 1)

    allocate(Rarray(dim,dim))

    do ii = 1,dim
      Rarray(ii,ii:dim) = RLarray(ii,ii:dim)
    end do
    
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
