!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
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
  real(dp):: Mpl, rpl, dPl, Pi

  Pi = acos(-1.0d0)

  emass=2.*boxlen*0.5d0**nlevelmax
  xmass=boxlen/2.0
  ymass=boxlen/2.0
  zmass=boxlen/2.0
  !rcut =boxlen/2.0*0.75

  Mpl=gravity_params(1)
  rpl=gravity_params(2)
  dmass=1.0/emass/(1.0+emass)**2

  !write(*,*)"Rho analytical is being called", Mpl, rpl

  do i=1,ncell
     rx=x(i,1)
     ry=x(i,2)
     rz=x(i,3)
     rr=sqrt(rx**2+ry**2+rz**2)
     dPl=Mpl*3.0d0/(4*Pi*rpl**3)
     if (rr .le. rpl) then
	!write(*,*)"Density in Msun/Kpc^3 is", dpl
        d(i)=dPl 
     else 
        d(i) = 0.d0
     end if 
  end do

end subroutine rho_ana

