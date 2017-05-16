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


end module eqsolver
