!==================================================================================
!=== Hydro ic for Toomre/exponential radial density profile galaxies.           ===
!=== Sets up a galactic binary merger.                                          ===
!=== Author : D.Chapon                                                          ===
!=== Date : 2010/09/01                                                          ===
!==================================================================================
module merger_parameters!{{{
  use amr_commons
  !--------------------!
  ! Galactic merger IC !
  !--------------------!

  ! Gas disk masses, given in GMsun in namelist,
  ! then converted in user unit.
  ! !!!!!! The galaxy #1 must be the heaviest !!!!!
  real(dp)::Mgas_disk1 = 2.0D1
  real(dp)::Mgas_disk2 = 2.0D1 ! Set to 0.0 for isolated galaxy runs
  ! Galactic centers, 0-centered, given in user unit
  real(dp), dimension(3)::gal_center1 = 0.0D0
  ! Set gal_center2 to values larger than 2*Lbox for isolated galaxy runs
  real(dp), dimension(3)::gal_center2 = 0.0D0
  ! Galactic disks rotation axis
  real(dp), dimension(3)::gal_axis1= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::gal_axis2= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::gal_axsph1= (/ 0.0D0, 0.0D0, 0.D0/)
  real(dp), dimension(3)::gal_axsph2= (/ 0.0D0, 0.0D0, 0.D0/)
  ! Particle ic ascii file names for the galaxies.
  !~~~~~~~~~~~~~~~~~~~~~~~WARNING~~~~~~~~~~~~~~~~~~~~~~~!
  ! Assumed to be J=z-axis, 0-centered galaxy ic files. !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::ic_part_file_gal1='ic_part1'
  character(len=256)::ic_part_file_gal2='ic_part2'
  ! Rotation curve files for the galaxies
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Those circular velocity files must contain :                              !
  ! - Column #1 : radius (in kpc)                                              !
  ! - Column #2 : circular velocity (in km/s)                                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::Vcirc_dat_file1='Vcirc1.dat'
  character(len=256)::Vcirc_dat_file2='Vcirc2.dat'
  ! Galactic global velocities, given in km/s in namelist, 
  ! then converted in user unit.
  real(dp), dimension(3)::Vgal1 = 0.0D0
  real(dp), dimension(3)::Vgal2 = 0.0D0
  ! Gas disk typical/truncation radius/height, given in kpc in namelist,
  ! then converted in user unit.
  real(dp)::typ_radius1 = 3.0D0
  real(dp)::typ_radius2 = 3.0D0
  real(dp)::cut_radius1 = 10.0D0
  real(dp)::cut_radius2 = 10.0D0
  real(dp)::typ_height1 = 1.5D-1
  real(dp)::typ_height2 = 1.5D-1
  real(dp)::cut_height1 = 4.0D-1
  real(dp)::cut_height2 = 4.0D-1
  real(dp)::typ_aratio1 = 0.0D-1
  real(dp)::typ_aratio2 = 0.0D-1
  real(dp)::typ_power1 = 0.0D-1
  real(dp)::typ_power2 = 0.0D-1
  real(dp)::hskpc_min = 0.0D0
  real(dp)::kgasrad = 1.0D0
  real(dp)::kgashsc = 1.0D0

  real(dp)::typ_radius_inner1 = 2.15d0 !IB/IT added block.
  real(dp)::mass_frac_inner1 = 0.0d0
  real(dp)::mass_frac_outer1 = 1.0d0
  real(dp)::typ_radius_inner2 = 2.15d0 !IB added block.
  real(dp)::mass_frac_inner2 = 0.0d0
  real(dp)::mass_frac_outer2 = 1.0d0
    
  ! Gas density profile : radial =>'exponential' (default) or 'Toomre'
  !                       vertical => 'exponential' or 'gaussian'
  character(len=16)::rad_profile='exponential'
  character(len=16)::z_profile='sech_sq'
  ! Inter galactic gas density contrast factor
  real(dp)::IG_density_factor = 1.0D-5
  ! added by IT 2015-11-24--2015-12-07
  real(dp),dimension(2) :: scale_vcirc=1.d0
  real(dp),dimension(2) :: scmin_vcirc=0.5d0
  real(dp),dimension(2) :: scmax_vcirc=2.0d0
  real(dp),dimension(2) :: power_vcirc=0.d0
  real(dp),dimension(2) :: scale_hs=1.0d0
  real(dp),dimension(2) :: scale_hspwr=0.0d0
  real(dp) :: scale_vigal=1.d0,scale_a2=1.d0
  real(dp) :: T2_ISM=100000.d0
  ! added by IT 2017-06-01
  real(dp) :: devflat_dens=1.d0,devscal_dens=0.d0,&
       devflat_vel=1.d0,devscal_vel=0.d0,scale_objsize=1.d0
  integer :: flg_qqm3d=-1!<0:off, 0:dens, 1:vel, 2:dens+vel
  integer :: flg_dmin=0!<0:original, 0:default, >0:IG=dmin
  integer :: flg_axsph=0!
  integer :: flg_extrcfile=2!0:just r and vcirc, 1:also hscale, 2:also gradhscale
  integer :: flg_varhscale=1
  
  !~~~~~~~~~~~~~~~~~~~~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~!
  ! For isolated galaxy runs :                           !
  ! --------------------------                           !
  ! - set 'Mgas_disk2' to 0.0                            !
  ! - set 'gal_center2' to values larger than Lbox*2     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
end module merger_parameters!}}}

module merger_commons

  use merger_parameters
  integer :: Vcirc_ncols=2
  real(dp), dimension(:,:), allocatable::Vcirc_dat1, Vcirc_dat2

end module merger_commons

subroutine read_merger_params! {{{
  use merger_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: i,nradii=100,igal !added by IT 2018-12-18
  real(dp) :: rfrac,rdum1,rdum2,vdum1,vdum2,hsdum1,hsdum2,hgraddum1,hgraddum2
  logical::nml_ok=.true.
  character(LEN=80)::infile

  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  real(dp)::norm_u
  logical::vcirc_file1_exists, vcirc_file2_exists
  logical::ic_part_file1_exists, ic_part_file2_exists
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  !added by IT 2018-10-09
  real(dp)::arot1,brot1,crot1,arot2,brot2,crot2
  real(dp),dimension(3)::vdum0=(/0.d0, 0.d0, 1.d0/)
  real(dp) :: d2r=1.745329251994329577d-2
  

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/merger_params/ IG_density_factor &
       & ,gal_center1, gal_center2, Mgas_disk1, Mgas_disk2 &
       & ,typ_radius1, typ_radius2, cut_radius1, cut_radius2 &
       & ,typ_height1, typ_height2, cut_height1, cut_height2 &
       & ,typ_radius_inner1,mass_frac_inner1,typ_aratio1,typ_power1 &
       & ,typ_radius_inner2,mass_frac_inner2,typ_aratio2,typ_power2 &
       & ,hskpc_min,rad_profile, z_profile, Vcirc_dat_file1, Vcirc_dat_file2 &
       & ,kgasrad,kgashsc,scale_hs,scale_hspwr &
       & ,ic_part_file_gal1, ic_part_file_gal2 &
       & ,gal_axis1, gal_axis2, gal_axsph1, gal_axsph2, Vgal1, Vgal2 &
! added by IT 2015-11-24--2017-06-02--20180702
       & ,scale_vcirc,scmin_vcirc,scmax_vcirc,power_vcirc,scale_vigal,scale_a2,T2_ISM &
       & ,devflat_dens,devflat_vel,devscal_dens,devscal_vel &
       & ,scale_objsize,flg_qqm3d,flg_dmin,flg_axsph,flg_extrcfile,flg_varhscale


  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=merger_params,END=106)
  goto 107
106 write(*,*)' You need to set up namelist &MERGER_PARAMS in parameter file'
  call clean_stop
107 continue
  close(1)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !-------------------------------------------------
  ! This section deals with the galactic merger initial conditions
  !-------------------------------------------------
  ! Gaseous disks masses : GMsun => user unit
  Mgas_disk1 = Mgas_disk1 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  Mgas_disk2 = Mgas_disk2 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  ! Galaxy global velocities : km/s => user unit
  Vgal1 = Vgal1 * 1.0D5 / scale_v
  Vgal2 = Vgal2 * 1.0D5 / scale_v
  ! Gas disk typical radii (a) : kpc => user unit
  typ_radius1 = typ_radius1 * 3.085677581282D21 / scale_l
  typ_radius2 = typ_radius2 * 3.085677581282D21 / scale_l
  typ_radius_inner1 = typ_radius_inner1 * 3.085677581282D21 / scale_l !Added by IB/IT.
  typ_radius_inner2 = typ_radius_inner2 * 3.085677581282D21 / scale_l !Added by IB.
  ! Gas disk max radii : kpc => user unit
  cut_radius1 = cut_radius1 * 3.085677581282D21 / scale_l
  cut_radius2 = cut_radius2 * 3.085677581282D21 / scale_l
  ! Gas disk typical thicknesses (h) : kpc => user unit
  typ_height1 = typ_height1 * 3.085677581282D21 / scale_l
  typ_height2 = typ_height2 * 3.085677581282D21 / scale_l
  ! Gas disk max thicknesses(zmax) : kpc => user unit
  cut_height1 = cut_height1 * 3.085677581282D21 / scale_l
  cut_height2 = cut_height2 * 3.085677581282D21 / scale_l
  mass_frac_outer1 = 1.0d0 - mass_frac_inner1 !Safety catch added by IB/IT.
  mass_frac_outer2 = 1.0d0 - mass_frac_inner2 !Safety catch added by IB.

  !copy to dummy variables in amr_commons to avoid circular dependencies of .mod files upon compilation
  axis1_dummy = gal_axis1
  axis2_dummy = gal_axis2
  
  !Galactic rot axis definition
  if (flg_axsph>=1) then
     arot1=gal_axsph1(3)*d2r!phase1
     brot1=gal_axsph1(2)*d2r!latitude1
     crot1=gal_axsph1(1)*d2r!longitude1
     arot2=gal_axsph2(3)*d2r!phase2
     brot2=gal_axsph2(2)*d2r!latitude2
     crot2=gal_axsph2(1)*d2r!longitude2
     !overwrite gal_axis*
     call rotzxz(vdum0,gal_axis1,arot1,brot1,crot1)
     call rotzxz(vdum0,gal_axis2,arot2,brot2,crot2)
  endif

  select case (rad_profile)
     case ('Toomre')
        if(myid==1) write(*,*) "Chosen hydro radial density profile :'Toomre'"
     case ('exponential')
        if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
     case ('double_exp') !IB added the double_exp case.
        if(myid==1) write(*,*) "Chosen hydro radial density profile :'double_exp'"
     case default
        if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
        rad_profile='exponential'
  end select
  select case (z_profile)
     case ('gaussian')
        if(myid==1) write(*,*) "Chosen hydro vertical density profile :'gaussian'"
     case ('exponential')
        if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
     case ('sech_sq') !IB added the sech_sq case.
        if(myid==1) write(*,*) "Chosen hydro vertical density profile :'sech_sq'"
     case default
        if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
        z_profile='sech_sq'
  end select
  if (myid==1) then
     write(*,*) 'Radial Vcirc scaling params:'
     write(*,*) 'scale_vcirc(:) =',scale_vcirc(:)
     write(*,*) 'power_vcirc(:) =',power_vcirc(:)
     write(*,*) 'scmin_vcirc(:) =',scmin_vcirc(:)
     write(*,*) 'scmax_vcirc(:) =',scmax_vcirc(:)
  endif
  if(Mgas_disk2 .GT. Mgas_disk1)then
     if(myid==1)write(*,*)'Error: The galaxy #1 must be bigger than #2'
     nml_ok=.false.
  endif
  ! Check whether velocity radial profile files exist or not.
  inquire(file=trim(Vcirc_dat_file1), exist=vcirc_file1_exists)
  if(.NOT. vcirc_file1_exists) then
     if(myid==1)write(*,*)'Error: Vcirc_dat_file1 ''', trim(Vcirc_dat_file1), ''' doesn''t exist '
     nml_ok=.false.
  end if
  inquire(file=trim(Vcirc_dat_file2), exist=vcirc_file2_exists)
  if(.NOT. vcirc_file2_exists) then
     if(myid==1)write(*,*)'Error: Vcirc_dat_file2 ''', trim(Vcirc_dat_file2), ''' doesn''t exist '
     nml_ok=.false.
  end if
  if ((vcirc_file1_exists) .AND. (vcirc_file2_exists)) then
     call read_vcirc_files()
     !display radial profile params
     if (myid==1) then
        open(187,file='radprofile.log')
        do i=0,nradii
           rfrac=dble(i)/nradii
           rdum1=rfrac*cut_radius1
           rdum2=rfrac*cut_radius2
           call get_Vcirc(rdum1, 1, Vdum1, hsdum1, hgraddum1, scale_l)
           call get_Vcirc(rdum2, 2, Vdum2, hsdum2, hgraddum2, scale_l)
           !print '(a,6(1x,f12.6))','r1,r2, hs1,hs2, hsgrad1,hsgrad2 =',rdum1,rdum2,hsdum1,hsdum2,hgraddum1,hgraddum2
           write(187,'(a,6(1x,f12.6))') 'r1,r2, hs1,hs2, hsgrad1,hsgrad2 =',rdum1,rdum2,hsdum1,hsdum2,hgraddum1,hgraddum2
        enddo
        close(187)
     endif
  end if

  ! Check whether ic_part files exist or not.
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal1), exist=ic_part_file1_exists)
  if(.NOT. ic_part_file1_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file1 ''', trim(ic_part_file_gal1), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal2), exist=ic_part_file2_exists)
  if(.NOT. ic_part_file2_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file2 ''', trim(ic_part_file_gal2), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if

  ! galactic rotation axis
  norm_u = sqrt(gal_axis1(1)**2 + gal_axis1(2)**2 + gal_axis1(3)**2)
  if(norm_u .EQ. 0.0D0) then
     if(myid==1)write(*,*)'Error: Galactic axis(1) is zero '
     nml_ok=.false.
  else
    if(norm_u .NE. 1.0D0) then
       gal_axis1 = gal_axis1 / norm_u
    end if
  end if
  norm_u = sqrt(gal_axis2(1)**2 + gal_axis2(2)**2 + gal_axis2(3)**2)
  if(norm_u .EQ. 0.0D0) then
     if(myid==1)write(*,*)'Error: Galactic axis(2) is zero '
     nml_ok=.false.
  else
    if(norm_u .NE. 1.0D0) then
       gal_axis2 = gal_axis2 / norm_u
    end if
  end if

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

end subroutine read_merger_params
! }}}

include 'qqm4ramses.f90'

subroutine conjure_qqm3d(mode)
  use qqm3d_module
  use storage_module
  integer :: mode
  if (mode<0) then
     call close_afield
     call close_storage
     return
  endif
  call init_qqm3d(76,0)
  call init_storage(0)
  print *,'filled array 0 (density)'
  if (mode==0.or.mode==10.or.mode==20) return!density only
  call drive_qqm3d
  call init_storage(1)
  print *,'filled array 1 (vx)'
  call drive_qqm3d
  call init_storage(2)
  print *,'filled array 2 (vy)'
  call drive_qqm3d
  call init_storage(3)
  print *,'filled array 3 (vz)'
  return
end subroutine conjure_qqm3d

subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use merger_commons
!  use qqm3d_module, only : qqm_nx,qqm_ny,qqm_nz,afield
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i, ind_gal
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::v,M,rho,dzz,zint, HH, rdisk, dpdr,dmax
  real(dp)::r, rr, rr1, rr2, abs_z
  real(dp), dimension(3)::vgal, axe_rot, xx1, xx2, xx, xx_rad
  real(dp), dimension(3)::objcenter=(/0.d0,0.d0,0.d0/),devvel !IT2017
  real(dp)::objsize,devdens !IT2017
  real(dp)::rgal, sum,sum2,dmin,dmin1,dmin2,zmin,zmax,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::M_b,az,eps, uu, u_1, u_2 !uu, u_1 and u_2 added by IB.
  real(dp)::rmin,rmax,a2,aa,Vcirc,hscale,hscale1,hscale2, hgradscale,hgradscale1,hgradscale2,&
       sigmasq, sigmasqgrad, HH_max
  real(dp)::rho_0_1, rho_0_2, rho_0, weight, da1, Vrot, Sigma_0, Sigma_0_1, Sigma_0_2
  real(dp)::fac_vcirc
  logical, save:: init_nml=.false.
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined merger parameters from the namelist
  if (.not. init_nml) then
      call read_merger_params
      if (flg_qqm3d>=0) call conjure_qqm3d(flg_qqm3d)
      init_nml = .true.
  end if

  if(maxval(abs(gal_center1)) .GT. (boxlen/2.0D0))then
     write(*,*)'Error: galactic center (1) coordinates must be in the box [', &
     & (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif
  if((maxval(abs(gal_center2)) .GT. (boxlen/2.0D0)) .AND. (Mgas_disk2 .NE. 0.0D0))then
     write(*,*)'Error: galactic center (2) coordinates must be in the box [', &
     & (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif

  if (T2_ISM<=0.d0) then
     a2=T2_star / scale_T2 * scale_a2! sound speed squared times user-defined scaling factor (default=1.0)
  else
     a2=T2_ISM / scale_T2! sound speed squared based on ISM polytropic temp
  endif
  aa=sqrt(a2)
  !print *,'>>>> a2,aa =',a2,aa
  ! Galactic central gas densities
  pi=dacos(-1.0D0)
  select case (rad_profile)
      case ('exponential')
         rho_0_1 = 1.0D0 - exp(-cut_radius1 / typ_radius1) * (1.0D0 + cut_radius1 / typ_radius1)
         rho_0_2 = 1.0D0 - exp(-cut_radius2 / typ_radius2) * (1.0D0 + cut_radius2 / typ_radius2)
      case ('double_exp') !IB added the double_exp case.
         !rho_0_1 = 1.0D0 - exp(-cut_radius1 / typ_radius1) * (1.0D0 + cut_radius1 / typ_radius1)
         u_1 = cut_radius1 / typ_radius_inner1
         u_2 = cut_radius1 / typ_radius1
         rho_0_1 = 1.0D0 - mass_frac_inner1*exp(-u_1)*(1.0D0 + u_1) &
              - mass_frac_outer1*exp(-u_2)*(1.0D0 + u_2)
         u_1 = cut_radius2 / typ_radius_inner2
         u_2 = cut_radius2 / typ_radius2
         rho_0_2 = 1.0D0 - mass_frac_inner2*exp(-u_1)*(1.0D0 + u_1) &
              - mass_frac_outer2*exp(-u_2)*(1.0D0 + u_2)
      case ('Toomre')
         rho_0_1 = sqrt(1.0D0 + (cut_radius1/typ_radius1)**2) - 1.0D0
         rho_0_2 = sqrt(1.0D0 + (cut_radius2/typ_radius2)**2) - 1.0D0
  end select
  select case (z_profile)
      case ('exponential')
          !rho_0_1 = rho_0_1 * (1.0D0 - exp(-cut_height1 / typ_height1))
          !rho_0_2 = rho_0_2 * (1.0D0 - exp(-cut_height2 / typ_height2))
      case ('gaussian')
          !rho_0_1 = rho_0_1 * (dsqrt(pi/2.0D0) * erf(cut_height1 / (dsqrt(2.0D0)*typ_height1)))
          !rho_0_2 = rho_0_2 * (dsqrt(pi/2.0D0) * erf(cut_height2 / (dsqrt(2.0D0)*typ_height2)))
      case ('sech_sq') !IB added the sech_sq case.
          !rho_0_1 = rho_0_1 * tanh(cut_height1/typ_height1)
          !rho_0_2 = rho_0_2 * tanh(cut_height2/typ_height2)
  end select
  !rho_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1**2 * typ_height1 * rho_0_1)
  !rho_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2**2 * typ_height2 * rho_0_2)

  select case (rad_profile) !IT added this block to replace former rho_0_1 evaluation method.
      case ('exponential')!IT modified exponential case
         !rho_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1*typ_radius1 * typ_height1 * rho_0_1)
         Sigma_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1*typ_radius1 * rho_0_1)
         rho_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1*typ_radius1 * typ_height1 * rho_0_1)
      case ('Toomre')
         rho_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1*typ_radius1 * typ_height1 * rho_0_1)
      case ('double_exp')
         uu = mass_frac_inner1*typ_radius_inner1*typ_radius_inner1 + mass_frac_outer1*typ_radius1*typ_radius1
         Sigma_0_1 = Mgas_disk1 / (4.0D0 * pi * uu * rho_0_1)
         rho_0_1 = Mgas_disk1 / (4.0D0 * pi * uu * typ_height1 * rho_0_1)
  end select

  select case (rad_profile) !IB added this block to replace former rho_0_2 evaluation method.
      case ('exponential')
         !rho_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2*typ_radius2 * typ_height2 * rho_0_2)
         Sigma_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2*typ_radius2 * rho_0_2)
         rho_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2*typ_radius2 * typ_height2 * rho_0_2)
      case ('Toomre')
         rho_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2*typ_radius2* typ_height2 * rho_0_2)
      case ('double_exp')
         uu = mass_frac_inner2*typ_radius_inner2*typ_radius_inner2 + mass_frac_outer2*typ_radius2*typ_radius2
         Sigma_0_2 = Mgas_disk2 / (4.0D0 * pi * uu * rho_0_2)
         rho_0_2 = Mgas_disk2 / (4.0D0 * pi * uu * typ_height2 * rho_0_2)
  end select

  ! Intergalactic gas density
  select case (rad_profile)
      case ('exponential')
          !dmin1 = rho_0_1 * exp(-cut_radius1 / typ_radius1)
          !dmin2 = rho_0_2 * exp(-cut_radius2 / typ_radius2)
          dmin1 = Sigma_0_1 * exp(-cut_radius1 / typ_radius1)
          dmin2 = Sigma_0_2 * exp(-cut_radius2 / typ_radius2)
      case ('double_exp') !double_exp case added by IB.
          u_1 = cut_radius1 / typ_radius_inner1
          u_2 = cut_radius1 / typ_radius1
          dmin1 = Sigma_0_1*(mass_frac_inner1*exp(-u_1) + mass_frac_outer1*exp(-u_2))
          u_1 = cut_radius2 / typ_radius_inner2
          u_2 = cut_radius2 / typ_radius2
          dmin2 = Sigma_0_2*(mass_frac_inner2*exp(-u_1) + mass_frac_outer2*exp(-u_2))
      case ('Toomre')
          dmin1 = rho_0_1 / sqrt(1.0D0 + (cut_radius1/typ_radius1)**2)
          dmin2 = rho_0_2 / sqrt(1.0D0 + (cut_radius2/typ_radius2)**2)
  end select
       
  if (flg_varhscale>=1) then
     call get_Vcirc(cut_radius1, 1, Vcirc, hscale1, hgradscale1, scale_l)
     call get_Vcirc(cut_radius1, 2, Vcirc, hscale2, hgradscale2, scale_l)
  else
     hscale1=typ_height1
     hscale2=typ_height2
  endif
  
  select case (z_profile)
      case ('exponential')
         !dmin1 = dmin1 * exp(-cut_height1 / typ_height1)
         !dmin2 = dmin2 * exp(-cut_height2 / typ_height2)
         !call get_Vcirc(cut_radius1, 1, Vcirc, hscale, hgradscale1, scale_l)
         dmin1 = dmin1 * exp(-cut_height1 / hscale1)
         !call get_Vcirc(cut_radius2, 2, Vcirc, hscale, hgradscale2, scale_l)
         dmin2 = dmin2 * exp(-cut_height2 / hscale2)
      case ('gaussian')
         !dmin1 = dmin1 * exp(-0.5D0 * (cut_height1 / typ_height1)**2)
         !dmin2 = dmin2 * exp(-0.5D0 * (cut_height2 / typ_height2)**2)
         !call get_Vcirc(cut_radius1, 1, Vcirc, hscale, hgradscale, scale_l)
         dmin1 = dmin1 * exp(-0.5D0 * (cut_height1 / hscale1)**2)
         !call get_Vcirc(cut_radius2, 2, Vcirc, hscale, hgradscale, scale_l)
         dmin2 = dmin2 * exp(-0.5D0 * (cut_height2 / hscale2)**2)
      case ('sech_sq') !IB added the sech_sq case.
         !call get_Vcirc(cut_radius1, 1, Vcirc, hscale, hgradscale, scale_l)
         dmin1 = dmin1 / (hscale1*(cosh(cut_height1 / hscale1))**2)
         !call get_Vcirc(cut_radius2, 2, Vcirc, hscale, hgradscale, scale_l)
         dmin2 = dmin2 / (hscale2*(cosh(cut_height2 / hscale2))**2)
  end select

  if(Mgas_disk2 .NE. 0.0D0) then
      dmin = min(dmin1, dmin2)
  else
      dmin = dmin1
  end if
  dmin = IG_density_factor * dmin

  !print '(a,2(2x,e12.4))','-------- CONDINIT: dmin =',dmin
  
  ! Loop over cells
  do i=1,nn
     xx1=x(i,:)-(gal_center1(:)+boxlen/2.0D0)
     xx2=x(i,:)-(gal_center2(:)+boxlen/2.0D0)

    ! Compute angular velocity

    ! Distance between cell and both galactic centers
     rr1 = norm2(xx1)
     rr2 = norm2(xx2)

    ! Projected cell position over galactic centers axis
     da1 = dot_product(gal_center2 - gal_center1, xx1) / norm2(gal_center2 - gal_center1)**2

     if(da1 .LT. (Mgas_disk1 / (Mgas_disk1 + Mgas_disk2))) then ! cell belongs to galaxy #1
         ind_gal = 1
         rr = rr1
         xx = xx1
         axe_rot = gal_axis1
         vgal = Vgal1
         rgal = typ_radius1
         rdisk = cut_radius1
         HH = typ_height1
         HH_max = cut_height1
         rho_0 = rho_0_1
         Sigma_0 = Sigma_0_1
     else ! cell belongs to galaxy #2
         ind_gal = 2
         rr = rr2
         xx = xx2
         axe_rot = gal_axis2
         vgal = Vgal2
         rgal = typ_radius2
         rdisk = cut_radius2
         HH = typ_height2
         HH_max = cut_height2
         rho_0 = rho_0_2
         Sigma_0 = Sigma_0_2
     end if
     objsize=2.d0*rdisk*scale_objsize
     !Random perturbation section
     if (flg_qqm3d==0.or.flg_qqm3d==2.or.flg_qqm3d==10.or.flg_qqm3d==12 &
          .or.flg_qqm3d==20.or.flg_qqm3d==22) then
        call probe_qqm3d(1,0,xx(1),xx(2),xx(3),objsize,objcenter,devflat_dens,devscal_dens,devdens)
        devdens=max(devdens,1.d-30)
     endif
     if (flg_qqm3d==1.or.flg_qqm3d==2) then!absolute perturbation
        call probe_qqm3d(1,1,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(1))
        call probe_qqm3d(1,2,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(2))
        call probe_qqm3d(1,3,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(3))
        devvel(1:3)=devvel(1:3)*1.0D5/scale_v!km/s to code units
     else if (flg_qqm3d==11.or.flg_qqm3d==12) then!scaled perturbation
        call probe_qqm3d(1,1,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(1))
        call probe_qqm3d(1,2,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(2))
        call probe_qqm3d(1,3,xx(1),xx(2),xx(3),objsize,objcenter,0.d0,devscal_vel,devvel(3))
     else if (flg_qqm3d==21.or.flg_qqm3d==22) then!relative perturbation
        call probe_qqm3d(1,1,xx(1),xx(2),xx(3),objsize,objcenter,devflat_vel,devscal_vel,devvel(1))
        call probe_qqm3d(1,2,xx(1),xx(2),xx(3),objsize,objcenter,devflat_vel,devscal_vel,devvel(2))
        call probe_qqm3d(1,3,xx(1),xx(2),xx(3),objsize,objcenter,devflat_vel,devscal_vel,devvel(3))
     endif
     !----

     ! Cylindric radius : distance between the cell and the galactic rotation axis
     xx_rad = xx - dot_product(xx,axe_rot) * axe_rot
     r = norm2(xx_rad)

     ! vertical position absolute value
     abs_z = sqrt(rr*rr - r*r)

     !print *,'**** condition =',((r-dx/2.0D0).lt.rdisk) .and. ((abs_z-dx/2.0D0) .lt. HH_max)

     if(((r-dx/2.0D0).lt.rdisk) .and. ((abs_z-dx/2.0D0) .lt. HH_max)) then
        ! Cell in the disk : analytical density profile + rotation velocity
        weight = (min(r+dx/2.0D0,rdisk)-(r-dx/2.0D0))/dx
        if (weight .NE. 1.0D0) then
            r = r + (weight-1.0D0)*dx/2.0D0
        end if
        ! Circular velocity
        !Vcirc= find_Vcirc(r, ind_gal)


        
        call get_Vcirc(r, ind_gal, Vcirc, hscale, hgradscale, scale_l)
        
        if (flg_varhscale<=0) then
           hscale=HH
        endif
        
        !sigmasq=sigmasq*scale_a2; sigmasqgrad=sigmasqgrad*scale_a2
        !print *,'>>>> a2, sigmasq, Vcirc =',a2, sigmasq, Vcirc
        !print *,'**** Vcirc, sigmasq, sigmasqgrad =',Vcirc, sigmasq, sigmasqgrad
        !call clean_stop
        weight = weight*(min(abs_z+dx/2.0D0,HH_max)-(abs_z-dx/2.0D0))/dx
        ! Density
        select case (rad_profile)
            case ('exponential')
                !q(i,1)= rho_0 * exp(-r / rgal)
                q(i,1)= Sigma_0 * exp(-r / rgal)
            case ('double_exp') !IB added the double_exp case.
                if (ind_gal .eq. 1) then
                    u_1 = r/typ_radius_inner1
                    u_2 = r/typ_radius1
                    q(i, 1) = Sigma_0*(mass_frac_inner1*exp(-u_1) + mass_frac_outer1*exp(-u_2))
                else
                    u_1 = r/typ_radius_inner2
                    u_2 = r/typ_radius2
                    q(i, 1) = Sigma_0*(mass_frac_inner2*exp(-u_1) + mass_frac_outer2*exp(-u_2))
                end if
            case ('Toomre')
                q(i,1)= rho_0 / sqrt(1.0D0 + (r/rgal)**2)
        end select
        select case (z_profile)
            case ('exponential')
                !q(i,1)= q(i,1) * exp(-abs_z / HH)
                q(i,1)= q(i,1) * exp(-abs_z / hscale)
            case ('gaussian')
                !q(i,1)= q(i,1) * exp(-0.5D0 * (abs_z / HH)**2)
                q(i,1)= q(i,1) * exp(-0.5D0 * (abs_z / hscale)**2)
            case ('sech_sq') !IB added the sech_sq case.
                q(i,1)= q(i,1) /(hscale*(cosh(abs_z/hscale))**2)
        end select
        
        q(i,1) = max(weight * q(i,1), dmin)
        if (flg_qqm3d==0.or.flg_qqm3d==2.or.flg_qqm3d==10.or.flg_qqm3d==12) then
           q(i,1)=q(i,1)*devdens
        endif
        ! P=rho*cs^2
        q(i,ndim+2)=a2*q(i,1)
        !q(i,ndim+2)=sigmasq*q(i,1)
        ! V = Vrot * (u_rot^xx_rad)/r + Vx_gal        
        !  -> Vrot = sqrt(Vcirc*Vcirc - 3*Cs^2 + r/rho * grad(rho) * Cs)
        select case (rad_profile)
            case ('exponential')
               !Vrot = sqrt(max(Vcirc*Vcirc - 3.0D0*a2 - r/rgal * a2,0.0D0))
               !Vrot = sqrt(max(Vcirc*Vcirc - sigmasq*r/rgal + r*sigmasqgrad, 0.0D0))
               Vrot = sqrt(max(Vcirc*Vcirc - a2*(kgasrad*r/rgal + kgashsc*r/hscale*hgradscale), 0.0D0))
            case ('double_exp') !IB added the double_exp case and used the above formula.
               if (ind_gal .eq. 1) then
                  u_1 = r/typ_radius_inner1
                  u_2 = r/typ_radius1
                  uu = mass_frac_inner1*u_1*exp(-u_1) + mass_frac_outer1*u_2*exp(-u_2)
                  uu = uu/(mass_frac_inner1*exp(-u_1) + mass_frac_outer1*exp(-u_2))
                  !Vrot = sqrt(max(Vcirc*Vcirc - 3.0D0*a2 - uu * a2,0.0D0))
                  !Vrot = sqrt(max(Vcirc*Vcirc - sigmasq*uu + r*sigmasqgrad, 0.0D0))
                  Vrot = sqrt(max(Vcirc*Vcirc - a2*(kgasrad*uu + kgashsc*r/hscale*hgradscale), 0.0D0))
               else
                  u_1 = r/typ_radius_inner2
                  u_2 = r/typ_radius2
                  uu = mass_frac_inner2*u_1*exp(-u_1) + mass_frac_outer2*u_2*exp(-u_2)
                  uu = uu/(mass_frac_inner2*exp(-u_1) + mass_frac_outer2*exp(-u_2))
                  !Vrot = sqrt(max(Vcirc*Vcirc - 3.0D0*a2 - uu * a2,0.0D0))
                  !Vrot = sqrt(max(Vcirc*Vcirc - sigmasq*uu + r*sigmasqgrad, 0.0D0))
                  Vrot = sqrt(max(Vcirc*Vcirc - a2*(kgasrad*uu + kgashsc*r/hscale*hgradscale), 0.0D0))
               end if
            case ('Toomre')
               Vrot = sqrt(max(Vcirc*Vcirc - 3.0D0*a2 - r*r/(r*r + rgal*rgal)*a2,0.0D0))
        end select
        fac_vcirc = scale_vcirc(ind_gal)**(max(r/rdisk,1.d-12)**power_vcirc(ind_gal))
        fac_vcirc = max(min(fac_vcirc, scmax_vcirc(ind_gal)), scmin_vcirc(ind_gal))
        Vrot = weight * Vrot * fac_vcirc!scale_vcirc factor added by IT 2015-11-24
        !Vrot = Vrot * scale_vcirc(ind_gal)!scale_vcirc factor added by IT 2015-11-24
        q(i,ndim-1:ndim+1) = Vrot * vect_prod(axe_rot,xx_rad)/r + vgal
        if (flg_qqm3d==1.or.flg_qqm3d==2) then!abs. perturbation
           q(i,ndim-1:ndim+1)=q(i,ndim-1:ndim+1)+devvel(1:3)
        else if (flg_qqm3d==11.or.flg_qqm3d==12) then!scaled perturbation
           q(i,ndim-1:ndim+1)=q(i,ndim-1:ndim+1)+Vrot*devvel(1:3)
        else if (flg_qqm3d==21.or.flg_qqm3d==22) then!rel. perturbation *EXPERIMENTAL*
           q(i,ndim-1:ndim+1)=q(i,ndim-1:ndim+1)*devvel(1:3)
        endif
        
!        if(metal)q(i,6)=z_ave*0.02 ! z_ave is in solar units
        !if(metal)q(i,6)=z_ave*0.02*10.**(0.5-r/cut_radius1) ! z_ave is in solar units
        !IT/IB 2019-11-14: z_ave set to 1.0
        !Melallicity used as a tag for the galaxies; Andromeda assumed to have
        !twice the metallicity of the MW.
        if (metal) then
           if (ind_gal .eq. 1) then
              q(i,6)=z_ave*0.02*0.2d0
           else
              q(i,6)=z_ave*0.02*0.1d0
           endif
        endif
    else ! Cell out of the gaseous disk : density = peanut and velocity = V_gal
       if (flg_dmin>=1) then
          q(i,1)=dmin
       else if (flg_dmin==0) then
          q(i,1)=1d-9/scale_nH
       else
          q(i,1)=1d-6/scale_nH
       endif
!        q(i,1)=dmin
        q(i,ndim+2)=1d7/scale_T2*q(i,1)
        ! V = Vgal
        q(i,ndim-1:ndim+1)= vgal*scale_vigal
        if(metal)q(i,6)=0.0
     endif
  enddo
  
!  if (flg_qqm3d>=0) call conjure_qqm3d(-1)

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum : Omega = rho * V
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif

  ! kinetic energy
  ! Total system global velocity : 0
  u(1:nn,ndim+2)=0.0D0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  ! E = Ec + P / (gamma - 1)
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do


contains
!{{{
!function find_Vcirc(rayon, indice)
!implicit none
!real(dp), intent(in)		:: rayon
!integer, intent(in)		:: indice
!real(dp)			:: find_Vcirc
!real(dp)			:: vitesse, rayon_bin, vitesse_old, rayon_bin_old
!integer				:: k, indmax
!
!k=2
!if (indice .EQ. 1) then
!	indmax = size(Vcirc_dat1,1)
!	rayon_bin = Vcirc_dat1(k,1)
!	rayon_bin_old = Vcirc_dat1(k-1,1)
!	vitesse = Vcirc_dat1(k,2)
!	vitesse_old = Vcirc_dat1(k-1,2)
!else
!	indmax = size(Vcirc_dat2,1)
!	rayon_bin = Vcirc_dat2(k,1)
!	rayon_bin_old = Vcirc_dat2(k-1,1)
!	vitesse = Vcirc_dat2(k,2)
!	vitesse_old = Vcirc_dat2(k-1,2)
!end if
!do while (rayon .GT. rayon_bin)
!	if(k .GE. indmax) then
!		write(*,*) "Hydro IC error : Radius out of rotation curve !!!"
!		call clean_stop
!	end if
!	k = k + 1
!	if (indice .EQ. 1) then
!		vitesse_old = vitesse
!		vitesse = Vcirc_dat1(k,2)
!		rayon_bin_old = rayon_bin
!		rayon_bin = Vcirc_dat1(k,1)
!	else
!		vitesse_old = vitesse
!		vitesse = Vcirc_dat2(k,2)
!		rayon_bin_old = rayon_bin
!		rayon_bin = Vcirc_dat2(k,1)
!	end if
!end do
!
!find_Vcirc = vitesse_old + (rayon - rayon_bin_old) * (vitesse - vitesse_old) / (rayon_bin - rayon_bin_old)
!
!return
!
!end function find_Vcirc


function vect_prod(a,b)
implicit none
real(dp), dimension(3), intent(in)::a,b
real(dp), dimension(3)::vect_prod

vect_prod(1) = a(2) * b(3) - a(3) * b(2)
vect_prod(2) = a(3) * b(1) - a(1) * b(3)
vect_prod(3) = a(1) * b(2) - a(2) * b(1)

end function vect_prod



function norm2(x)
implicit none
real(dp), dimension(3), intent(in)::x
real(dp) :: norm2

norm2 = sqrt(dot_product(x,x))

end function norm2
!}}}

end subroutine condinit


subroutine get_Vcirc(radius, indice, vcirc, hscale, hgradscale, scale_l)
  use merger_commons!, only : Vcirc_dat1, Vcirc_dat2
  implicit none
  real(dp), intent(in)		:: radius
  integer, intent(in)		:: indice
  real(dp), intent(out)		:: vcirc, hscale, hgradscale, scale_l!, sigmasq, sigmasqgrad
  real(dp)			:: v_c, radius_bin, v_c_old, radius_bin_old, h, h_old, hgrad, hgrad_old
  real(dp)                      :: rscale,h_aratio,h_power,hdum
  !real(dp)			:: ssq, ssq_old, ssqgrad, ssqgrad_old
  integer			:: k, indmax
  
  k=2
  if (indice .EQ. 1) then
     indmax = size(Vcirc_dat1,1)
     radius_bin = Vcirc_dat1(k,1)
     radius_bin_old = Vcirc_dat1(k-1,1)
     v_c = Vcirc_dat1(k,2)
     v_c_old = Vcirc_dat1(k-1,2)
     h = Vcirc_dat1(k,3)
     h_old = Vcirc_dat1(k-1,3)
     hgrad = Vcirc_dat1(k,4)
     hgrad_old = Vcirc_dat1(k-1,4)
     h_aratio = typ_aratio1
     h_power = typ_power1
     rscale = typ_radius1 !neglect typ_radius_inner* here
     hdum = typ_height1
     !ssq = Vcirc_dat1(k,3)
     !ssq_old = Vcirc_dat1(k-1,3)
     !ssqgrad = Vcirc_dat1(k,4)
     !ssqgrad_old = Vcirc_dat1(k-1,4)
  else
     indmax = size(Vcirc_dat2,1)
     radius_bin = Vcirc_dat2(k,1)
     radius_bin_old = Vcirc_dat2(k-1,1)
     v_c = Vcirc_dat2(k,2)
     v_c_old = Vcirc_dat2(k-1,2)
     h = Vcirc_dat2(k,3)
     h_old = Vcirc_dat2(k-1,3)
     hgrad = Vcirc_dat2(k,4)
     hgrad_old = Vcirc_dat2(k-1,4)
     h_aratio = typ_aratio2
     h_power = typ_power2
     rscale = typ_radius2 !neglect typ_radius_inner* here
     hdum = typ_height2
     !ssq = Vcirc_dat2(k,3)
     !ssq_old = Vcirc_dat2(k-1,3)
     !ssqgrad = Vcirc_dat2(k,4)
     !ssqgrad_old = Vcirc_dat2(k-1,4)
  end if
  do while (radius .GT. radius_bin)
     if(k .GE. indmax) then
        write(*,*) "Hydro IC error : Radius out of rotation curve !!!"
        call clean_stop
     end if
     k = k + 1
     if (indice .EQ. 1) then
        radius_bin_old = radius_bin
        radius_bin = Vcirc_dat1(k,1)
        v_c_old = v_c
        v_c = Vcirc_dat1(k,2)
        h_old = h
        h = Vcirc_dat1(k,3)
        hgrad_old = hgrad
        hgrad = Vcirc_dat1(k,4)
        !ssq_old = ssq
        !ssq = Vcirc_dat1(k,3)
        !ssqgrad_old = ssqgrad
        !ssqgrad = Vcirc_dat1(k,4)
     else
        radius_bin_old = radius_bin
        radius_bin = Vcirc_dat2(k,1)
        v_c_old = v_c
        v_c = Vcirc_dat2(k,2)
        h_old = h
        h = Vcirc_dat2(k,3)
        hgrad_old = hgrad
        hgrad = Vcirc_dat2(k,4)
        !ssq_old = ssq
        !ssq = Vcirc_dat2(k,3)
        !ssqgrad_old = ssqgrad
        !ssqgrad = Vcirc_dat2(k,4)
     end if
  end do

  vcirc = v_c_old + (radius - radius_bin_old) * (v_c - v_c_old) / (radius_bin - radius_bin_old)
  if (h_aratio>0.d0) then !power-law flaring profile
     hscale = h_aratio*rscale * (radius/rscale)**h_power
     hgradscale = h_aratio * h_power * (radius/rscale)**(h_power-1.d0)
  else
     if (flg_extrcfile>=1) then! flaring profile from dice rotation curve file
        hscale = h_old + (radius - radius_bin_old) * (h-h_old) / (radius_bin-radius_bin_old)
        if (flg_extrcfile>=2) then
           hgradscale = hgrad_old + (radius-radius_bin_old) * (hgrad-hgrad_old) / (radius_bin-radius_bin_old)
        else
           hgradscale = 0.d0
        endif
     else! no flaring
        hscale = hdum
        hgradscale = 0.d0
     endif
     if (scale_hs(indice)>0.d0) then!optional flare-data correction
        hscale = hscale*scale_hs(indice)*(radius/rscale)**scale_hspwr(indice)
     endif
  endif

  if (h_aratio>0.d0 .or. flg_extrcfile>=1) then
     !hskpc_min does not make sense for constant, user-defined hscale
     if (hskpc_min<0.d0) then
        hscale = max(hscale, abs(hskpc_min)*3.085677581282D21/scale_l)
     else if (hskpc_min>0.d0) then
        hscale = sqrt(hscale*hscale + (hskpc_min*3.085677581282D21/scale_l)**2)
     endif
  endif
  !sigmasq = ssq_old + (radius - radius_bin_old) * (ssq - ssq_old) / (radius_bin - radius_bin_old)
  !sigmasqgrad = ssqgrad_old + (radius - radius_bin_old) * (ssqgrad - ssqgrad_old) / (radius_bin - radius_bin_old)

  return
end subroutine get_Vcirc

!------------------------------------------------------------------------------------- 
! Circular velocity files reading
! {{{

subroutine read_vcirc_files
  use merger_commons
  implicit none
  integer:: nvitesses, ierr, i
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Array management
  Vcirc_ncols = flg_extrcfile+2
  
  ! Galaxy #1
  nvitesses = 0
  open(unit=123, file=trim(Vcirc_dat_file1), iostat=ierr)
  do while (ierr==0)
	read(123,*,iostat=ierr)
	if(ierr==0) then
		nvitesses = nvitesses + 1  ! Number of samples
	end if
  end do
  allocate(Vcirc_dat1(nvitesses,4))
  Vcirc_dat1 = 0.0D0
  rewind(123)
  do i=1,nvitesses
	read(123,*) Vcirc_dat1(i,1:Vcirc_ncols)
  end do
  close(123)
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat1(:,1) = Vcirc_dat1(:,1) * 3.085677581282D21 / scale_l
  Vcirc_dat1(:,2) = Vcirc_dat1(:,2) * 1.0D5 / scale_v
  Vcirc_dat1(:,3) = Vcirc_dat1(:,3) * 3.085677581282D21 / scale_l
  !Vcirc_dat1(:,4) = Vcirc_dat1(:,4) ! no scaling required for dh/dr
  !Vcirc_dat1(:,3) = Vcirc_dat1(:,3) * 1.0D10 / (scale_v*scale_v)
  !Vcirc_dat1(:,4) = Vcirc_dat1(:,4) * 1.0D10 / (scale_v*scale_v) * scale_l/3.085677581282D21

  ! Galaxy #2
  nvitesses = 0
  open(unit=123, file=trim(Vcirc_dat_file2), iostat=ierr)
  do while (ierr==0)
	read(123,*,iostat=ierr)
	if(ierr==0) then
		nvitesses = nvitesses + 1 ! Number of samples
	end if
  end do
  allocate(Vcirc_dat2(nvitesses,4))
  Vcirc_dat2 = 0.0D0
  rewind(123)
  do i=1,nvitesses
	read(123,*) Vcirc_dat2(i,1:Vcirc_ncols)
  end do
  close(123)
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat2(:,1) = Vcirc_dat2(:,1) * 3.085677581282D21 / scale_l
  Vcirc_dat2(:,2) = Vcirc_dat2(:,2) * 1.0D5 / scale_v
  Vcirc_dat2(:,3) = Vcirc_dat2(:,3) * 3.085677581282D21 / scale_l
  !Vcirc_dat2(:,4) = Vcirc_dat2(:,4) ! no scaling required for dh/dr
  !Vcirc_dat2(:,3) = Vcirc_dat2(:,3) * 1.0D10 / (scale_v*scale_v)
  !Vcirc_dat2(:,4) = Vcirc_dat2(:,4) * 1.0D10 / (scale_v*scale_v) * scale_l/3.085677581282D21

end subroutine read_vcirc_files
! }}}
!--------------------------------------------------------------------------------------
!Vector rotation routines
!==== vec3rot
SUBROUTINE vxrot(vec0,vec1,phi)
!---- Drehung um x-Achse
  IMPLICIT NONE
  REAL(KIND=8) :: vec0(1:3),vec1(1:3),phi,c,s
  c = cos(phi)
  s = sin(phi)
  vec1(1) =   vec0(1)
  vec1(2) =             c*vec0(2) - s*vec0(3)
  vec1(3) =             s*vec0(2) + c*vec0(3)
  RETURN
END SUBROUTINE vxrot

SUBROUTINE vyrot(vec0,vec1,phi)
  !---- Drehung um y-Achse
  IMPLICIT NONE
  REAL(KIND=8) :: vec0(1:3),vec1(1:3),phi,c,s
  c = cos(phi)
  s = sin(phi)
  vec1(1) = c*vec0(1)           + s*vec0(3)
  vec1(2) =             vec0(2)
  vec1(3) =-s*vec0(1)           + c*vec0(3)
  RETURN
END SUBROUTINE vyrot

SUBROUTINE vzrot(vec0,vec1,phi)
!---- Drehung um z-Achse
  IMPLICIT NONE
  REAL(KIND=8) :: vec0(1:3),vec1(1:3),phi,c,s
  c = cos(phi)
  s = sin(phi)
  vec1(1) = c*vec0(1) - s*vec0(2)
  vec1(2) = s*vec0(1) + c*vec0(2)
  vec1(3) =                         vec0(3)
  RETURN
END SUBROUTINE vzrot

SUBROUTINE rotzxz(vec0,vec1,alpha,beta,gamma)
!---- Euler rotation extrinsic forward
!     (swap alpha and gamma for intrinsic rotation)
  IMPLICIT NONE
  REAL(KIND=8) :: vec0(1:3),vec1(1:3),alpha,beta,gamma
  REAL(KIND=8) :: veca(1:3),vecb(1:3)
  CALL vzrot(vec0,veca,alpha)
  CALL vxrot(veca,vecb,beta)
  CALL vzrot(vecb,vec1,gamma)
  RETURN
END SUBROUTINE rotzxz
