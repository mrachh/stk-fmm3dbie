      implicit real *8 (a-h,o-z)
      real *8, allocatable :: amat(:,:),bmat(:,:),rhs(:),soln(:)
      real *8, allocatable :: errs(:),soltrue(:),wtmp(:)

      call prini(6,13)


      did = 1.1d0
      
      n = 2000

      dd = 0.1d0
      nsys = 9*n
      n3 = 3*n

      allocate(amat(nsys,n3),bmat(n3,nsys),rhs(nsys),soln(nsys))
      allocate(soltrue(nsys),wtmp(n3))

      do i=1,nsys
        do j=1,n3
          bmat(j,i) = (hkrand(0)-0.5d0)/sqrt(n3+0.0d0)
        enddo
        rhs(i) = hkrand(0)
        soltrue(i) = hkrand(0)
        soln(i) = 0
      enddo

      do i=1,n3
        do j=1,nsys
          amat(j,i) = dd*(hkrand(0)-0.5d0)/sqrt(nsys+0.0d0)
        enddo
      enddo

      call dmatvec(n3,nsys,bmat,soltrue,wtmp)
      call dmatvec(nsys,n3,amat,wtmp,rhs)
      do i=1,nsys
        rhs(i) = rhs(i) + did*soltrue(i)
      enddo

      numit = 100
      allocate(errs(numit+1))
      rres = 0
      eps_gmres = 1.0d-10
      call dgmres_stok_blas(n,amat,bmat,did,numit,rhs,eps_gmres,
     1  niter,errs,rres,soln)
      
      erra = 0
      ra = 0
      do i=1,nsys
        ra = ra + abs(soltrue(i))**2
        erra = erra + abs(soltrue(i) - soln(i))**2
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in solution=*',erra,1)


      



      

      

      

      stop
      end
