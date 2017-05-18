module eqsolver
  use accuracy
  implicit none

contains
  subroutine ludecompose(dim, array)
    real(dp), intent(inout) :: array(:,:)
    integer, intent(in) :: dim
    real(dp), allocatable :: tmp(:,:), onemat(:,:)
    integer :: nmax, ii, jj

    allocate(onemat(dim,dim))
    onemat = 0
    do ii=1,dim
      onemat(ii,ii)=1
    end do

    !Rarray in array(dim+1:2*dim,:)
    array(dim+1:2*dim,:) = array(1:dim,:)
    !Larray in array(2*dim+1:3*dim,:)
    array(2*dim+1:3*dim,:) = onemat
    !Parray in array(3*dim+1:4*dim,:)
    array(3*dim+1:4*dim,:) = onemat

    allocate(tmp(1,dim))
    gauss: do ii = 1, dim
      
      nmax = ii - 1 + maxloc(abs(array(dim+ii:2*dim,ii)),dim=1) !pivot

      tmp(1,:) = array(dim+ii,:)
      array(dim+ii,:) = array(dim+nmax,:)
      array(dim+nmax,:) = tmp(1,:)
      
      tmp(1,:) = array(3*dim+ii,:)
      array(3*dim+ii,:) = array(3*dim+nmax,:)
      array(3*dim+nmax,:) = tmp(1,:)

      do jj = ii + 1, dim
        tmp(1,1)= array(dim+jj,ii)/array(dim+ii,ii)
        
        array(dim+jj,:) = array(dim+jj,:)-(array(dim+ii,:) * tmp(1,1))
        array(2*dim+jj,:) = array(2*dim+jj,:)-(array(2*dim+ii,:) * tmp(1,1))       
        
      end do
      
      if (abs(array(dim+ii,ii)) < 1E-12_dp) then
        write(*,*) "underterminated System"
        return
      end if  
    end do gauss

    deallocate(onemat)
    deallocate(tmp)

  end subroutine ludecompose


  
  subroutine substituteback(dim, array, bvec, res)

    real(dp), intent(inout) :: array(:,:)
    real(dp), allocatable, intent(inout) :: bvec(:)
    real(dp), intent(out), allocatable :: res(:)
    integer, intent(in) :: dim
    real(dp), allocatable :: preres(:), bveccalc(:)
    integer :: ii
    
    allocate(res(dim))
    allocate(preres(dim))
    allocate(bveccalc(dim))
    
    bveccalc = matmul(array(3*dim+1:4*dim,:), bvec)

    !Forward
    preres(1) = bveccalc(1)
    do ii = 2, dim
      preres(ii) = bveccalc(ii) - dot_product(array(2*dim+ii,1:ii-1),preres(1:ii-1))
    end do

    !Backward
    res(dim) = preres(dim)/array(2*dim,dim)
    do ii = dim-1, 1, -1
      res(ii) = (preres(ii) - dot_product(array(dim+ii,ii+1:dim),&
          &res(ii+1:dim)))/array(dim+ii,ii)
    end do

    deallocate(bvec)
    deallocate(preres)
    deallocate(bveccalc)
  end subroutine substituteback
  

end module eqsolver
