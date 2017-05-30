module eqsolver
  use accuracy
  implicit none

contains
  subroutine ludecompose(array, RLarray, Parray)
    real(dp), intent(in) :: array(:,:)
    real(dp), allocatable, intent(out) :: RLarray(:,:), Parray(:,:)
    real(dp), allocatable :: tmp(:,:), onemat(:,:)
    integer :: nmax, ii, jj, dim

    dim = size(array,1)
    allocate(onemat(dim,dim))
    onemat = 0
    do ii=1,dim
      onemat(ii,ii)=1
    end do

    allocate(RLarray(dim,dim))
    allocate(Parray(dim,dim))
    RLarray = array
    Parray = onemat

    allocate(tmp(1,dim))
    gauss: do ii = 1, dim
      
      nmax = ii - 1 + maxloc(abs(RLarray(ii:dim,ii)),dim=1) !pivot

      tmp(1,:) = RLarray(ii,:)
      RLarray(ii,:) = RLarray(nmax,:)
      RLarray(nmax,:) = tmp(1,:)
      
      tmp(1,:) = Parray(ii,:)
      Parray(ii,:) = Parray(nmax,:)
      Parray(nmax,:) = tmp(1,:)

      tmp(1,1:ii-1) = RLarray(ii,1:ii-1)
      RLarray(ii,1:ii-1) = RLarray(nmax,1:ii-1)
      RLarray(nmax,1:ii-1) = tmp(1,1:ii-1)

      do jj = ii + 1, dim
        tmp(1,1)= RLarray(jj,ii)/RLarray(ii,ii)
        
        RLarray(jj,:) = RLarray(jj,:)-(RLarray(ii,:) * tmp(1,1))
        RLarray(jj,ii) = RLarray(jj,ii)+(tmp(1,1))       
        
      end do
      
      
      if (abs(RLarray(ii,ii)) < 1E-12_dp) then
        write(*,*) "underterminated System"
        return
      end if  
    end do gauss

  end subroutine ludecompose


  
  subroutine substituteback(RLarray, Parray, bvec, res)

    real(dp), intent(in) :: RLarray(:,:), Parray(:,:), bvec(:,:)
    real(dp), intent(out), allocatable :: res(:,:)
    integer :: dim, bdim
    real(dp), allocatable :: preres(:,:), bveccalc(:,:)
    integer :: ii

    dim = size(RLarray, 1)
    bdim = size(bvec, 2)
    allocate(res(dim,bdim))
    allocate(preres(dim,bdim))
    allocate(bveccalc(dim,bdim))
    
    bveccalc = matmul(Parray(1:dim,:), bvec)

    !Forward
    preres(1,:) = bveccalc(1,:)
    do ii = 2, dim
      preres(ii,:) = bveccalc(ii,:) - matmul(RLarray(ii,1:ii-1),preres(1:ii-1,:))
    end do
    
    !Backward
    res(dim,:) = preres(dim,:)/RLarray(dim,dim)
    do ii = dim-1, 1, -1
      res(ii,:) = (preres(ii,:) - matmul(RLarray(ii,ii+1:dim),&
          &res(ii+1:dim,:)))/RLarray(ii,ii)
    end do

  end subroutine substituteback
  

end module eqsolver
