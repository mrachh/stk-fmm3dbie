



!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:6) = tanget vector 1
!     *info(7:9) = tangent vector 2
!     *info(10:12) = normal vector
!
!

!      
!     We take the following conventions for the Stokes kernels
!
!     For a source y and target x, let r_i = x_i-y_i
!     and let r = sqrt(r_1^2 + r_2^2 + r_3^2)
!
!     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
!     (without the 1/4pi scaling) are
!
!     G_{ij}(x,y) = (r_i r_j)/(2r^3) + delta_{ij}/(2r)
!     P_j(x,y) = -r_j/r^3
!
!     The (Type I) stresslet, T_{ijk}, and its associated pressure
!     tensor, PI_{jk}, (without the 1/4pi scaling) are
!     
!     T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5
!     PI_{jk} = -2 delta_{jk} + 6 r_j r_k/r^5      
!
!


subroutine st3d_slp_grad_comp(src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)

  complex *16 :: zk
  real *8 :: val, dr(3),over4pi
  data over4pi/0.07957747154594767d0/
!f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val

  !
  ! returns one entry of the gradient of the 
  ! Stokes single layer potential
  ! kernel
  !
  !
  ! G_ij = (0.5*(targ_i-src_i)(targ_j-src_j)/|src-targ|^3 +
  !               0.5*delta_ij/|src-targ|)/4/pi
  !
  ! Returns val(i)= x^{\ell} y^{m} z^{n}/|src-targ|^5/8/pi
  !    with \ell + m + n = 3, and coeffs are ordered as
  !    i, \ell, m, n
  !    1, 3, 0, 0
  !    2, 2, 1, 0
  !    3, 2, 0, 1
  !    4, 1, 2, 0
  !    5, 1, 1, 1
  !    6, 1, 0, 2
  !    7, 0, 3, 0
  !    8, 0, 2, 1
  !    9, 0, 1, 2
  !    10, 0, 0, 3
  

  i = ipars(1)
  j = ipars(2)
  k = ipars(3)
  
  dr(1)=targ(1)-src(1)
  dr(2)=targ(2)-src(2)
  dr(3)=targ(3)-src(3)

  
  r=sqrt(dx**2+dy**2+dz**2)
  rinv = 1.0d0/r
  rinv3 = 0.5d0*rinv**5
  

  val = (dr(1)**i)*(dr(2)**j)*(dr(3)**k)*rinv5*over4pi

  return
end subroutine st3d_slp_grad_comp


