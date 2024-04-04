!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,force,dx,ncell)
  use amr_parameters
  use poisson_parameters  
  use poisson_commons
  use mond_commons
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::force ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp) :: Force_constant, r_inv, cos_theta !Added by IB.
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,constant

  real(dp),dimension(1:nvector)::r ! Cell center position.

  ! Self-gravity
  ! TODO: Check boundary conditions
  if (gravity_type==0) then
     
     ! Compute position relative to center of mass
     r = 0.0d0
     do idim=1,ndim
        do i=1,ncell
           x(i,idim) = x(i,idim) - multipole(idim+1)/multipole(1)
           r(i) = r(i) + x(i,idim)**2
        end do
     enddo
     do i=1,ncell
         r(i)=sqrt(r(i))
     end do

     if (connected_Mond) then
        if (Activate_g_ext) then !Added by IB based on Arxiv: 1509.08457.
            constant = -multipole(1)*nu_ext
            do i = 1, ncell
                r_inv = 1/r(i)
                cos_theta = sum(g_ext_dir(:)*x(i, :))*r_inv
                force(i, :) = -2.0d0*u_Lambda_eff*x(i, :) + &
                constant*r_inv*r_inv*(K_0_ext*cos_theta*g_ext_dir(:) + r_inv*x(i, :)*(1.0d0 + 0.5d0*K_0_ext*(1.0d0 - 3.0d0*cos_theta*cos_theta)))
            end do
        else
        ! Compute MONDian acceleration a_i = -sqrt(G*M*a0)*x_i/r^2  (a ~ 1/r)
        constant = -sqrt(multipole(1)*a0)
        !do idim=1,ndim
        !   do i=1,ncell
        !      force(i,idim) = (const/r(i)**2 - 2.0d0*u_Lambda_eff)*x(i,idim) !Added by IB.
        !   end do
        !enddo
        do i = 1, ncell !Above block simplified by IB.
            Force_constant = const/r(i)**2 - 2.0d0*u_Lambda_eff !Added by IB.
            force(i, :) = Force_constant*x(i, :)
        end do
end if
     else
        ! Compute Newtonian acceleration a_i = -G*M*x_i/r^3  (a ~ 1/r^2)
        constant = -multipole(1)
        !do idim=1,ndim
        !   do i=1,ncell
        !      force(i,idim) = constant * x(i,idim) / r(i)**3
        !   end do
        !enddo
        do i = 1, ncell !Above block simplified by IB.
            Force_constant = constant/r(i)**3
            force(i, :) = Force_constant*x(i, :)
        end do
     endif
  endif

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           force(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     gmass=gravity_params(1) ! GM
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        force(i,1)=-gmass*rx/rr**3
#if NDIM>1
        force(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        force(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(xx, rr, pp, ngrid)
  use amr_commons
  use poisson_commons
  use mond_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  real(dp),dimension(1:nvector, 1:ndim)::xx !Added by IT.
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp) :: M_inv, r_inv, cos_theta !Added by IB.
  real(dp):: fourpi
  real(dp)::constant

#if NDIM==1
  fourpi=4.D0*ACOS(-1.0D0)
  constant = multipole(1)*fourpi/2d0
  do i=1,ngrid
     pp(i)=constant*rr(i)
  end do
#endif
#if NDIM==2
  constant = multipole(1)*2d0
  do i=1,ngrid
     pp(i)=constant*log(rr(i))
  end do
#endif
#if NDIM==3
  if (connected_Mond) then
     ! Compute the MONDian potential, phi = sqrt(G M a0) log(r)
     ! This approximation requires:
     ! (a) being in the deep MOND-limit
     ! (b) r being large enough [this is also the case for Newtonian dynamics]

if (Activate_g_ext) then !Added by IB based on Arxiv: 1509.08457.
    constant = -multipole(1)*nu_ext
    M_inv = 1.0d0/multipole(1)
    do i = 1,ngrid
        r_inv = 1/rr(i)
        cos_theta = sum(g_ext_dir(:)*(xx(i, :) - multipole(2:ndim+1)*M_inv))*r_inv
        pp(i) = constant*r_inv*(1.0d0 + 0.5d0*K_0_ext*(1 - cos_theta*cos_theta)) + u_Lambda_eff*rr(i)*rr(i)
    end do
else
     constant = sqrt(multipole(1)*a0)
     do i=1,ngrid
        pp(i)=constant*log(rr(i)) + u_Lambda_eff*rr(i)*rr(i) !Added by IB.
     end do
end if
  else
     ! Compute the Newtonian potential, phi = -G M/r
     constant = -multipole(1)
     do i=1,ngrid
        pp(i)=constant/rr(i)
     end do
  endif
#endif
end subroutine phi_ana