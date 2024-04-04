!########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use mond_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd,rcut
  real(dp):: Mpl, rmax, d_coma, Pi

  Pi = acos(-1.0d0)

  emass=2.*boxlen*0.5d0**nlevelmax
  xmass=boxlen/2.0
  ymass=boxlen/2.0
  zmass=boxlen/2.0
  !rcut =boxlen/2.0*0.75

  rmax=gravity_params(1)
  dmass=1.0/emass/(1.0+emass)**2

if (mond) then

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     rz=x(i,3)-zmass
     rr=sqrt(rx**2+ry**2+rz**2)
     d_coma=0.0795*rr*(1.841d-34+1.41d-36*rr + 2.41d-39*rr*rr+6.945d-42*rr**3)/((2.967d-25+3.014d-27*rr+3.89d-30*rr*rr)**2*(76176 + rr*rr)**2)
     if (rr .le. rmax) then
        d(i)=d_coma
     else 
        d(i) = 1.0d-5
     end if 
  end do

else

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     rz=x(i,3)-zmass
     rr=sqrt(rx**2+ry**2+rz**2)
     d_coma=(1.217d16 + 5.328d10*rr*rr)/(76176+rr*rr)**2
     if (rr .le. rmax) then
        d(i)=d_coma
     else 
        d(i) = 1.0d-5
     end if 
  end do
end if 
  

end subroutine rho_ana

