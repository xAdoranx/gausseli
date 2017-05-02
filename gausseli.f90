program gausseli
  Implicit none

  integer :: ii, jj, dim, nmax
  integer, allocatable :: num(:)
  integer, parameter :: dp = selected_real_kind(12,99)
  real(dp), allocatable :: array(:,:), arrayb(:,:), res(:), bvec(:), tmp(:,:)
  character(len=26) :: outputform = '(1000000F10.2)'

  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
   ! write(*,*) "To stop make Matrix Dimension negative or 0"
   ! write(*,*) "Please enter the Matrix Dimension of A in A*x=b"
    read(21,*) dim  
    if (dim <= 0) then
      exit mainloop
    end if
    
    allocate(array(dim,dim))
    allocate(bvec(dim))
    !write(*,*) "Please enter the Matrix values row by row:"
    read(21,*) array
    !write(*,*) "Please enter the values of the solution vector b:"
    read(21,*) bvec
    array = transpose(array)

    write(*,*) "Input Matrix:"
    do ii = 1, dim
      write(*,outputform) array(ii,:)
    end do

    allocate(arrayb(dim,dim+1))
    arrayb(:,1:dim)=array
    arrayb(:,dim+1)=bvec

    allocate(tmp(1,dim+1))
    gauss: do ii = 1, dim
      write(12,*) ii
      tmp(1,:) = arrayb(ii,:)
      nmax = ii - 1 + maxloc(abs(arrayb(ii:dim,ii)),dim=1) !pivot
      
      arrayb(ii,:) = arrayb(nmax,:)
      arrayb(nmax,:) = tmp(1,:)

      do jj = ii + 1, dim
        arrayb(jj,:) = arrayb(jj,:)-(arrayb(ii,:) * (arrayb(jj,ii)/arrayb(ii,ii)))
      end do
      
      if (abs(arrayb(ii,ii)) < 1E-12_dp) then
        write(*,*) "underterminated System"
        exit mainloop
      end if  
    end do gauss

    write(11,*) "Upper triangle matrix:"
    do ii = 1,dim
      write(11,outputform) arrayb(ii,:)
    end do
   
    allocate(res(dim))
    res(dim) = arrayb(dim,dim+1) / arrayb(dim,dim)
    fin: do ii = dim-1, 1, -1
      res(ii) = (arrayb(ii, dim + 1) - dot_product(arrayb(ii,ii+1:dim),res(ii:dim))) / arrayb(ii,ii)
    end do fin

    output: do ii = 1,dim
      write(11,"(A,I0,A2,F10.6)") "x_", ii, "=", res(ii)
    end do output

    deallocate(res)
    deallocate(tmp)
    deallocate(array)
    deallocate(arrayb)
    deallocate(bvec)
        
  end do mainloop
  close(11)
  close(21)

end program gausseli
