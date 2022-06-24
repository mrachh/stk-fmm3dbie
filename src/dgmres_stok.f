      subroutine dgmres_stok_blas(n,amat,bmat,did,numit,rhs,eps_gmres,
     1   niter,errs,rres,soln)
c
c
c  This sburoutine solves the (did *I + A*B)soln = rhs using 
c  GMRES which avoids low threshold stagnation when did \neq 0
c  and using blas routines to apply A and B. 
c
c  The matrices A, and B are assumed to be of size 9*n \times 3*n,
c  and 3*n \times 9*n respectively. 
c
c  The routine would use whatever version of blas the source
c  code is compiled with
c 
c  Input arguments:
c    * n: integer
c        number of points
c    * amat: real *8 (9*n,3*n) 
c        the operator A in the representation above
c    * bmat: real *8 (3*n,9*n)
c        the operator B in the representation above
c    * did: real *8
c        strength of the identity term
c    * numit: integer
c        max number of iterations
c    * rhs: real *8 (9*n)
c        the boundary data
c    * eps_gmres: real *8
c        gmres tolerance
c   
c  Output arugments:
c    * niter: integer
c        number of gmres iterations required for relative residual
c          
c    * errs: real *8(numit+1)
c        errs(1:niter) is the estimate of relative residual as 
c        a function of iteration number
c 
c    * rres: real *8
c        relative residual for computed solution
c              
c    * soln: real *8(9*n)
c        density which solves the linear system 
c
      implicit none
      integer, intent(in) :: n,numit
      real *8, intent(in) :: amat(9*n,3*n),bmat(3*n,9*n), did
      real *8, intent(in) :: rhs(9*n),eps_gmres
      integer, intent(out) :: niter
      real *8, intent(out) :: errs(numit+1),rres,soln(9*n)
c
c
c       gmres variables
c
      real *8 dtmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l,i,j,n3,nsys
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:),wtmp0(:)

      complex *16 ztmp

      nsys = 9*n
      n3 = 3*n


      allocate(vmat(nsys,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(nsys),svec(numit+1),yvec(numit+1))
      allocate(wtmp0(n3))


c
c     NOTE: matrix equation should be of the form (did*I + A*B)x = y
c       the identity scaling (z) is defined via did below,
c       and K represents the action of the principal value 
c       part of the matvec
c

      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
        errs(i) = 0
      enddo
      errs(numit+1) = 0


c
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,nsys
        rb = rb + abs(rhs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nsys
        vmat(i,1) = rhs(i)/rb
      enddo
C$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

        call dmatvec(n3,nsys,bmat,vmat(1,it),wtmp0)
        call dmatvec(nsys,n3,amat,wtmp0,wtmp)
        

        do k=1,it
          dtmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)          
          do j=1,nsys
            dtmp = dtmp + wtmp(j)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = dtmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,nsys
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,nsys
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,nsys
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        dtmp = wnrm2

        call rotmat_gmres(hmat(it,it),dtmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j=1,nsys
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,nsys
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        evaluation routine  
c
        call dmatvec(n3,nsys,bmat,soln,wtmp0)
        call dmatvec(nsys,n3,amat,wtmp0,wtmp)

            
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,nsys
            rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2
          enddo
C$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
c
      return
      end
