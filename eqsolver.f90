module eqsolver
  use accuracy
  implicit none

contains
  subroutine ludecompose(array, Rarray, Larray, Parray)
    real(dp), intent(in) :: array(:,:)
    real(dp), intent(out), allocatable :: Rarray(:,:), Larray(:,:), Parray(:,:)
    real(dp), allocatable :: tmp(:,:), onemat(:,:)
    integer :: dim, nmax, ii, jj

    dim = size(array,dim=1)
    allocate(Rarray(dim,dim))
    Rarray = array

    allocate(Larray(dim,dim))
    allocate(Parray(dim,dim))
    allocate(onemat(dim,dim))
    onemat = 0
    do ii=1,dim
      onemat(ii,ii)=1
    end do
    

    Larray = onemat
    Parray = onemat

    allocate(tmp(1,dim))
    gauss: do ii = 1, dim
      
      nmax = ii - 1 + maxloc(abs(Rarray(ii:dim,ii)),dim=1) !pivot

      tmp(1,:) = Rarray(ii,:)
      Rarray(ii,:) = Rarray(nmax,:)
      Rarray(nmax,:) = tmp(1,:)
      
      tmp(1,:) = Parray(ii,:)
      Parray(ii,:) = Parray(nmax,:)
      Parray(nmax,:) = tmp(1,:)

      do jj = ii + 1, dim
        tmp(1,1)= Rarray(jj,ii)/Rarray(ii,ii)
        
        Rarray(jj,:) = Rarray(jj,:)-(Rarray(ii,:) * tmp(1,1))
        Larray(jj,:) = Larray(jj,:)-(Larray(ii,:) * tmp(1,1))       
        
      end do
      
      if (abs(Rarray(ii,ii)) < 1E-12_dp) then
        write(*,*) "underterminated System"
        return
      end if  
    end do gauss

    deallocate(onemat)
    deallocate(tmp)

  end subroutine ludecompose


  
  subroutine substituteback(Rarray, Larray, Parray, bvec, res)

    real(dp), intent(in) :: Rarray(:,:), Larray(:,:), Parray(:,:), bvec(:)
    real(dp), intent(out), allocatable :: res(:)
    real(dp), allocatable :: preres(:), bveccalc(:)
    integer :: ii, dim

    dim = size(Rarray,dim=1)
    
    allocate(res(dim))
    allocate(preres(dim))
    allocate(bveccalc(dim))
    bveccalc = matmul(Parray, bvec)
    preres(1) = bveccalc(1)
    do ii = 2, dim
      preres(ii) = bveccalc(ii) - dot_product(Larray(ii,1:ii-1),preres(1:ii-1))
    end do
    res(dim) = preres(dim)/Rarray(dim,dim)
    do ii = dim-1, 1, -1
      res(ii) = (preres(ii) - dot_product(Rarray(ii,ii+1:dim),preres(ii+1:dim)))/Rarray(ii,ii)
    end do

    deallocate(preres)
    deallocate(bveccalc)
  end subroutine substituteback
  

end module eqsolver
