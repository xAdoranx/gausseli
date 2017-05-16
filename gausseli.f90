program gausseli
  use accuracy
  Implicit none

  integer :: ii, jj, dim, nmax
  integer, allocatable :: num(:)
  real(dp), allocatable :: onemat(:,:), Rarray(:,:), Larray(:,:), Parray(:,:), preres(:), res(:)&
      &, bvec(:), tmp(:,:)
  
  character(len=26) :: outputform = '(1000000F10.2)'

  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    read(21,*) dim  
    if (dim <= 0) then
      exit mainloop
    end if
    
    allocate(Rarray(dim,dim))
    allocate(bvec(dim))
    
    read(21,*) Rarray
    
    read(21,*) bvec
    Rarray = transpose(Rarray)

    write(*,*) "Input Matrix:"
    do ii = 1, dim
      write(*,outputform) Rarray(ii,:)
    end do

    allocate(onemat(dim,dim))
    allocate(Larray(dim,dim))
    allocate(Parray(dim,dim))
    onemat = 0
    do ii=1, dim
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
        cycle mainloop
      end if  
    end do gauss

    write(11,*) "Upper triangle matrix:"
    do ii = 1,dim
      write(11,outputform) Rarray(ii,:)
    end do
    write(11,*) "Lower triangle matrix:"
    do ii = 1,dim
      write(11,outputform) Larray(ii,:)
    end do
    
   
    allocate(res(dim))
    allocate(preres(dim))
    bvec = matmul(Parray, bvec)
    preres(1) = bvec(1)
    do ii = 2, dim
      preres(ii) = bvec(ii) - dot_product(Larray(ii,1:ii-1),preres(1:ii-1))
    end do
    res(dim) = preres(dim)/Rarray(dim,dim)
    do ii = dim-1, 1, -1
      res(ii) = (preres(ii) - dot_product(Rarray(ii,ii+1:dim),preres(ii+1:dim)))/Rarray(ii,ii)
      write(*,*) res(ii)
    end do
    

    output: do ii = 1,dim
      write(11,"(A,I0,A2,F10.6)") "x_", ii, "=", res(ii)
    end do output

    deallocate(res)
    deallocate(preres)
    deallocate(onemat)
    deallocate(tmp)
    deallocate(Rarray)
    deallocate(Larray)
    deallocate(Parray)
    deallocate(bvec)
        
  end do mainloop
  close(11)
  close(21)

end program gausseli
