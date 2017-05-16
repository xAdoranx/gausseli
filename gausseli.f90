program gausseli
  use accuracy
  use eqsolver
  Implicit none

  integer :: ii, jj, dim
  integer, allocatable :: num(:)
  real(dp), allocatable :: Rarray(:,:), Larray(:,:), Parray(:,:), preres(:), res(:), bvec(:),&
      & array(:,:)
  
  
  character(len=26) :: outputform = '(1000000F10.2)'

  open(21, file="gauss.inp", status="old", form="formatted", action="read")
  open(11, file="output.dat", status="replace", form="formatted", action="write")
  
  mainloop: do
    
    read(21,*) dim  
    if (dim <= 0) then
      exit mainloop
    end if

    allocate(array(dim,dim))
    allocate(bvec(dim))
    
    read(21,*) array
    
    read(21,*) bvec
    array = transpose(array)

    write(*,*) "Input Matrix:"
    do ii = 1, dim
      write(*,outputform) array(ii,:)
    end do

    call ludecompose(array, Rarray, Larray, Parray)

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
    deallocate(array)
    deallocate(Rarray)
    deallocate(Larray)
    deallocate(Parray)
    deallocate(bvec)
        
  end do mainloop
  close(11)
  close(21)

end program gausseli
