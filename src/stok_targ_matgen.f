



      subroutine getnearquad_stok_s_grad(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2     ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     3     iquad,rfac0,nquad,wnear,wnear_g)
c     
c       this subroutine generates the near field quadrature
c       for the velocity of the representation u = (alpha S + beta D)
c       where S is the Stokes single layer potential and D
c       the double layer potential.
c       the near field is specified by the user 
c       in row sparse compressed format.
c
c
c
c
c       The quadrature is computed by the following strategy
c        targets within a sphere of radius rfac0*rs
c        of a chunk centroid is handled using adaptive integration
c        where rs is the radius of the bounding sphere
c        for the patch
c  
c       All other targets in the near field are handled via
c        oversampled quadrature
c
c       The recommended parameter for rfac0 is 1.25d0
c
c       Note: the code currently only works for stokeslet,
c       Stresslet implementation of gradient is remaining
c
c
c              
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders - integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            starting location of data on patch i
c  
c         iptype - integer(npatches)
c           type of patch
c           iptype = 1 -> triangular patch discretized with RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c          srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c             srcvals(10:12,i) - normals info
c 
c         ndtarg - integer
c            leading dimension of target array
c        
c         ntarg - integer
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c         ipatch_id - integer(ntarg)
c            id of patch of target i, id = -1, if target is off-surface
c
c         uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if on surface, otherwise
c            set to 0 by default
c            
c          eps - real *8
c             precision requested
c
c           iquadtype - integer
c              quadrature type
c              iquadtype = 1, use ggq for self + adaptive integration
c                 for rest
c 
c
c           nnz - integer
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer
c               radius parameter for near field
c
c           nquad - integer
c               number of near field entries corresponding to
c               each source-target pair 
c
c        output
c            wnear - real *8(6,nquad)
c               the desired near field quadrature
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(6,nquad),wnear_g(10,nquad)

      real *8, allocatable :: wnear1(:)

      integer ipars(2), ijloc(2,6),ijkpow(3,10)
      integer ndd,ndz,ndi
      real *8 dpars
      complex *16 zpars
      real *8 alpha, beta
      integer i,j,ii,l
      integer ipv

      procedure (), pointer :: fker
      external st3d_slp, st3d_slp_grad_comp

c
c
c        initialize the appropriate kernel function
c



      fker => st3d_slp
      
      allocate(wnear1(nquad))

      ndd = 0
      ndi = 2
      ndz = 0

      ijloc(1,1) = 1
      ijloc(2,1) = 1
      ijloc(1,2) = 1
      ijloc(2,2) = 2
      ijloc(1,3) = 1
      ijloc(2,3) = 3
      ijloc(1,4) = 2
      ijloc(2,4) = 2
      ijloc(1,5) = 2
      ijloc(2,5) = 3
      ijloc(1,6) = 3
      ijloc(2,6) = 3
      
      do l = 1,6
         i = ijloc(1,l)
         j = ijloc(2,l)
         ipars(1) = i
         ipars(2) = j
         fker => st3d_slp
         ipv = 0
         
         call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1        ipatch_id,uvs_targ,
     1        eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,
     1        row_ptr,col_ind,iquad,rfac0,nquad,wnear1)
         print *, "done with kernel l=",l
         
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)       
         do ii=1,nquad
            wnear(l,ii) = wnear1(ii)
         enddo
C$OMP END PARALLEL DO        
      enddo

      ijkpow(1:3,1) = (/ 3,0,0 /)
      ijkpow(1:3,2) = (/ 2,1,0 /)
      ijkpow(1:3,3) = (/ 2,0,1 /)
      ijkpow(1:3,4) = (/ 1,2,0 /)
      ijkpow(1:3,5) = (/ 1,1,1 /)
      ijkpow(1:3,6) = (/ 1,0,2 /)
      ijkpow(1:3,7) = (/ 0,3,0 /)
      ijkpow(1:3,8) = (/ 0,2,1 /)
      ijkpow(1:3,9) = (/ 0,1,2 /)
      ijkpow(1:3,10) = (/ 0,0,3 /)

      ndi = 3
      ndd = 0
      ndz = 0


      
      do l = 1,10
         fker => st3d_slp_grad_comp
         ipv = 1
         
         call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1        ipatch_id,uvs_targ,
     1        eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ijkpow(1,l),nnz,
     1        row_ptr,col_ind,iquad,rfac0,nquad,wnear1)
         print *, "done with kernel l=",l
         
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)       
         do ii=1,nquad
            wnear_g(l,ii) = wnear1(ii)
         enddo
C$OMP END PARALLEL DO        
      enddo

      return
      end
c
c
c
c
c        
      subroutine stok_cg_matgen(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2    ipatch_id,uvs_targ,eps,xmat,gmat)
c
cf2py  intent(in) npatches,norders,ixyzs,iptype
cf2py  intent(in) npts,srccoefs,srcvals,eps,ndtarg,ntarg,targs
cf2py  intent(in) ipatch_id,uvs_targ
cf2py  intent(out) xmat,gmat
c
c     this subroutine returns the discretization matrix
c     for the on-surface discretization of the stokes combined
c     field layer potential. The unknowns are ordered as:
c
c     \sigma_{1}(x_{1}), \sigma_{2}(x_{1}), \sigma_{3}(x_{1})...
c      \sigma_{1}(x_{2}), \sigma_{2}(x_{2})l \sigma_{3}(x_{2})..
c                           .
c                           .
c                           .
c                           .
c      \sigma_{1}(x_{n}), \sigma_{2}(x_{n}), \sigma_{3}(x_{n})
c
c     And the same ordering for the velocity components as well
c
c     For the gradients, the ordering is \partial_{1} u_{1} (x_{1}),
c      \partial_{2} u_{1} (x_{1}), partial_{3} u_{1} x(x_{1}),
c      \partial_{1} u_{2} (x_{1}), .. and so on
c
c     For the gradient info and targets on surface, this code
c     only returns the principal value part of the itnegral
c
c
c     Representation:
c        u = S \sigma 
c     
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders- integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c 
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c         srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c 
c          eps - real *8
c             precision requested for computing quadrature and fmm
c             tolerance
c          ndtarg - integer
c            leading dimension of target array
c          ntarg - integer
c            number of targets
c          targs - real *8 (ndtarg,ntarg)
c            target information
c          ipatch_id - integer(ntarg)
c            patch id for on surface targets, unused otherwise
c          uvs_targ - real *8 (2,ntarg)
c            local uv coordinates of surface targets, unused otherwise
c
c         output
c           xmat - real *8(3*ntarg,3*npts)
c              discretization matrix for the stokes boundary value
c              problem
c 
c           gmat - real *8(9*ntarg,3*npts)
c              discretization matrix for the stokes boundary value
c              problem
c
c
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),ndtarg,ntarg,ipatch_id(ntarg)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      real *8 targs(ndtarg,ntarg),uvs_targ(2,ntarg)
      real *8 dpars(2)
      real *8 xmat(3*ntarg,3*npts)
      real *8 gmat(9*ntarg,3*npts)

      real *8 uint


      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear_s(:,:), wts(:),wnear_g(:,:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi,rsurf
      real *8 rfac,rfac0,alpha,beta
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

      real *8 did,ra
      integer jj,l,nmat
      real *8 w11,w12,w13,w21,w22,w23,w31,w32,w33,ww(10)
      

      
      nmat = 3*npts
      done = 1
      pi = atan(done)*4


c
c
c        this might need fixing
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)

      nnz = ntarg*npatches
      allocate(row_ptr(ntarg+1),col_ind(nnz))

      do i=1,ntarg+1
        row_ptr(i) = (i-1)*npatches + 1
      enddo


      do i=1,ntarg
        do j=1,npatches
          col_ind((i-1)*npatches + j) = j
        enddo
      enddo

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = 0
      
c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear_s(6,nquad))
      allocate(wnear_g(10,nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
      do i=1,nquad
         do j = 1,6
            wnear_s(j,i) = 0
         enddo

         do j=1,10
           wnear_g(j,i) = 0
         enddo
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_stok_s_grad(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1     ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1     iquad,rfac0,nquad,wnear_s,wnear_g)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      allocate(wts(npts))


      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)



      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,npols,jj,w11,w12,w13)
C$OMP$PRIVATE(w21,w22,w23,w31,w32,w33,ww)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
             jj = jstart + l-1
             w11 = wnear_s(1,jquadstart+l-1)
             w12 = wnear_s(2,jquadstart+l-1)
             w13 = wnear_s(3,jquadstart+l-1)
             w21 = w12
             w22 = wnear_s(4,jquadstart+l-1)
             w23 = wnear_s(5,jquadstart+l-1)
             w31 = w13
             w32 = w23
             w33 = wnear_s(6,jquadstart+l-1)
             xmat(3*(i-1)+1,3*(jj-1)+1) = w11
             xmat(3*(i-1)+1,3*(jj-1)+2) = w12 
             xmat(3*(i-1)+1,3*(jj-1)+3) = w13 

             xmat(3*(i-1)+2,3*(jj-1)+1) = w21 
             xmat(3*(i-1)+2,3*(jj-1)+2) = w22 
             xmat(3*(i-1)+2,3*(jj-1)+3) = w23 

             xmat(3*(i-1)+3,3*(jj-1)+1) = w31 
             xmat(3*(i-1)+3,3*(jj-1)+2) = w32 
             xmat(3*(i-1)+3,3*(jj-1)+3) = w33

             ww(1:10) = wnear_g(1:10, jquadstart+l-1)
c    handle x component of gradient             
             w11 = -2*ww(1) + ww(4) + ww(6)
             w12 = -2*ww(2) + ww(7) + ww(9)
             w13 = -2*ww(3) + ww(8) + ww(10)
             
             w21 = w12
             w22 = -ww(1) -4*ww(4) -ww(6)
             w23 = -3*ww(5)
             
             w31 = w13
             w32 = w23
             w33 = -ww(1) -ww(4) - 4*ww(6)

             gmat(9*(i-1)+1,3*(jj-1)+1) = w11
             gmat(9*(i-1)+1,3*(jj-1)+2) = w12
             gmat(9*(i-1)+1,3*(jj-1)+3) = w13
             
             gmat(9*(i-1)+4,3*(jj-1)+1) = w21
             gmat(9*(i-1)+4,3*(jj-1)+2) = w22
             gmat(9*(i-1)+4,3*(jj-1)+3) = w33
             
             gmat(9*(i-1)+7,3*(jj-1)+1) = w31
             gmat(9*(i-1)+7,3*(jj-1)+2) = w32
             gmat(9*(i-1)+7,3*(jj-1)+3) = w33

c
c   Now handle the y component of the gradient
c
             w11 = -4*ww(2) - ww(7) - ww(9)
             w12 = ww(1) - 2*ww(4) + ww(6)
             w13 = -3*ww(5)

             w21 = w12
             w22 = ww(2) -2*ww(7) + ww(9)
             w23 = ww(3) -2*ww(8) + ww(10)
             
             w31 = w13
             w32 = w23
             w33 = -ww(2) - ww(7) - 4*ww(9)

             gmat(9*(i-1)+2,3*(jj-1)+1) = w11
             gmat(9*(i-1)+2,3*(jj-1)+2) = w12
             gmat(9*(i-1)+2,3*(jj-1)+3) = w13
             
             gmat(9*(i-1)+5,3*(jj-1)+1) = w21
             gmat(9*(i-1)+5,3*(jj-1)+2) = w22
             gmat(9*(i-1)+5,3*(jj-1)+3) = w33
             
             gmat(9*(i-1)+8,3*(jj-1)+1) = w31
             gmat(9*(i-1)+8,3*(jj-1)+2) = w32
             gmat(9*(i-1)+8,3*(jj-1)+3) = w33

c
c  Now handle the z component
c
             w11 = -4*ww(3)-ww(8)-ww(10)
             w12 = -3*ww(5)
             w13 = ww(1) + ww(4)-2*ww(6)
             
             w21 = w12
             w22 = -ww(3) -4*ww(8)-ww(10)
             w23 = ww(2) + ww(7) - 2*ww(9)

             w31 = w13
             w32 = w23
             w33 = ww(3) + ww(8) -2*ww(10)
             
             gmat(9*(i-1)+3,3*(jj-1)+1) = w11
             gmat(9*(i-1)+3,3*(jj-1)+2) = w12
             gmat(9*(i-1)+3,3*(jj-1)+3) = w13
             
             gmat(9*(i-1)+6,3*(jj-1)+1) = w21
             gmat(9*(i-1)+6,3*(jj-1)+2) = w22
             gmat(9*(i-1)+6,3*(jj-1)+3) = w33
             
             gmat(9*(i-1)+9,3*(jj-1)+1) = w31
             gmat(9*(i-1)+9,3*(jj-1)+2) = w32
             gmat(9*(i-1)+9,3*(jj-1)+3) = w33
             

          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
     

     
c     
      return
      end
c     
c
c
