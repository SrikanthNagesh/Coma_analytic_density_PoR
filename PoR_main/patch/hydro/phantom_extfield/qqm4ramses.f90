MODULE rnd_array_module
  IMPLICIT NONE
  INTEGER :: rnd_ioseed=0,rnd_cmode=0,rnd_kseed=37,rnd_dimseed
  INTEGER :: rnd_callcount=0,rnd_initcount=0
  INTEGER :: rnd_verbose=0

CONTAINS
    
  SUBROUTINE wrap_get_rnd(nx,xxx)
! Version with seed hidden in common block
    IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
    INTEGER nx
    REAL(KIND=8) :: xxx(nx)
!      write(*,*) 'rnd_ioseed =',rnd_ioseed
    CALL get_rnd_array(rnd_ioseed,nx,xxx)
    RETURN
  END SUBROUTINE wrap_get_rnd

  SUBROUTINE wrap_one_rnd(x)
! Simple random number, no array
! Version with seed hidden in common block
    IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
    REAL(KIND=8) :: xxx(1),x
    CALL get_rnd_array(rnd_ioseed,1,xxx)
    x=xxx(1)
    RETURN
  END SUBROUTINE wrap_one_rnd

  SUBROUTINE wrap_ind_rnd(x,ix,irange)
! Simple integer random number, no array
! Version with seed hidden in common block
    IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
    INTEGER :: ix,irange,imax=2147483647
    REAL(KIND=8) :: xxx(1),x,dimax=4294967296.d0
    CALL get_rnd_array(rnd_ioseed,1,xxx)
    x=xxx(1)
    IF (irange==0) THEN
       ix=int(imax*x)
    ELSE IF (irange<0) THEN
       ix=floor((x-0.5d0)*dimax)
    ELSE
       ix=int(irange*x)
    ENDIF
    RETURN
  END SUBROUTINE wrap_ind_rnd

  SUBROUTINE get_rnd_array(iseed,nx,xxx)
    IMPLICIT NONE
    INTEGER :: iseed,nx
    REAL(KIND=8) :: xxx(nx),xdum(1)
    rnd_callcount=rnd_callcount+1
    IF (iseed.NE.-1) CALL init_random_seed(iseed)
    CALL RANDOM_NUMBER(xxx)
    IF (rnd_cmode>=1) THEN!experimental closed interval mode (req. one more rnd call)
       CALL RANDOM_NUMBER(xdum)
       IF (xdum(1)<0.5d0) THEN
          xxx(1:nx)=1.d0-xxx(1:nx)
       ENDIF
    ENDIF
    iseed=-1!seed is now initialized
    RETURN
  END SUBROUTINE get_rnd_array

  SUBROUTINE get_rnd_number(iseed,x)
! Simple random number, no array
    IMPLICIT NONE
    INTEGER :: iseed
    REAL(KIND=8) :: xxx(1),x
    CALL get_rnd_array(iseed,1,xxx)
    x=xxx(1)
    RETURN
  END SUBROUTINE get_rnd_number

  SUBROUTINE scale_rnd_array(nx,xsc,xlo,xhi)
! Version with seed hidden in common block
    IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
    INTEGER :: nx
    REAL(KIND=8) :: xxx(nx),xsc(nx),xlo,xhi
!      write(*,*) 'rnd_ioseed =',rnd_ioseed
    CALL get_rnd_array(rnd_ioseed,nx,xxx)
    xsc(:)=xlo+(xhi-xlo)*xxx(:)
    RETURN
  END SUBROUTINE scale_rnd_array

  SUBROUTINE init_random_seed(iseed)
    INTEGER :: i, n, clock, iseed
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    rnd_initcount=rnd_initcount+1
    CALL RANDOM_SEED(size = n)
    rnd_dimseed=n
    ALLOCATE(seed(n))
    if (rnd_verbose>=1) write(*,*) 'n =',n
    IF (iseed<=0) THEN
       CALL SYSTEM_CLOCK(COUNT=clock)
    ELSE
       clock = iseed
    ENDIF
!      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    seed = clock + rnd_kseed * (/ (i - 1, i = 1, n) /)
    if (rnd_verbose>=1) then
       do i=1,n
          print *,'i,seed(i) =',i,seed(i)
       enddo
    endif
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed

  SUBROUTINE iter_random_seed(iter,iseed)
    implicit none
    integer :: i,j,iter,iseed,oseed,n
    real(kind=8) :: dimax=4294967296.d0
    real(kind=8),dimension(:),allocatable :: xxx
    INTEGER, DIMENSION(:), ALLOCATABLE :: aseed
    call init_random_seed(iseed)
    print *,'iter,iseed,rnd_dimseed=',iter,iseed,rnd_dimseed
    iseed=-1
    allocate(aseed(1:rnd_dimseed))
    allocate(xxx(1:rnd_dimseed))
    do i=1,iter
       call RANDOM_NUMBER(xxx)
       do j=1,rnd_dimseed
          aseed(j)=floor((xxx(j)-0.5d0)*dimax)
       enddo
    enddo
    if (rnd_verbose>=1) then
       do j=1,rnd_dimseed
          print *,'i,aseed(j) =',j,aseed(j)
       enddo
    endif
    if (iter>0) call RANDOM_SEED(put=aseed)
    return
  END SUBROUTINE iter_random_seed

END MODULE rnd_array_module

MODULE deviates_module
  USE rnd_array_module
  IMPLICIT NONE
!... For Gaussian
  INTEGER :: iset=0
  REAL(KIND=8) :: gset,x1_gauss,x2_gauss
  LOGICAL :: force_iset=.false.
!... For Lorentzian
  REAL(KIND=8) :: xold=-1.d0

CONTAINS

  SUBROUTINE uniform_pm1(x)
    REAL(KIND=8) :: x
    CALL wrap_one_rnd(x)
    x=2.d0*(x-0.5d0)
    RETURN
  END SUBROUTINE uniform_pm1
  
  SUBROUTINE uniform_range(x,xmin,xmax)
    REAL(KIND=8) :: x,xmin,xmax
    CALL wrap_one_rnd(x)
    x=xmin+(xmax-xmin)*x
    RETURN
  END SUBROUTINE uniform_range

  SUBROUTINE uniform3d(length,vcen,v)
    REAL(KIND=8) :: length,vcen(3),v(3),x1,x2,x3
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    CALL wrap_one_rnd(x3)
    v(1)=length*(x1-0.5d0)+vcen(1)
    v(2)=length*(x2-0.5d0)+vcen(2)
    v(3)=length*(x3-0.5d0)+vcen(3)
    RETURN
  END SUBROUTINE uniform3d

  SUBROUTINE unisphere_deviate(radius,vcen,v)
    REAL(KIND=8) :: radius,vcen(3),v(3),x1,x2,x3,r,x,y,z
    REAl(KIND=8) :: tau=6.283185307179586476925d0! 2*pi
    REAl(KIND=8) :: third=1.d0/3.d0
    LOGICAL :: issrfc
    issrfc=radius<0
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    CALL wrap_one_rnd(x3)
    IF (issrfc) THEN
       r=1.d0
    ELSE
       r=x1**third
    ENDIF
    z=(1.d0-2.d0*x2)*r
    x=sqrt(r*r-z*z)*cos(tau*x3)
    y=sqrt(r*r-z*z)*sin(tau*x3)
    v(1)=abs(radius)*x+vcen(1)
    v(2)=abs(radius)*y+vcen(2)
    v(3)=abs(radius)*z+vcen(3)
    RETURN
  END SUBROUTINE unisphere_deviate

  SUBROUTINE gauss_deviate(xgauss)
!    USES gfortran intrinsic PRNG
!    integer :: aseed(12),i!debug
    REAL(KIND=8) :: xgauss
    REAL(KIND=8) :: fac,rsq,v1,v2
    REAL(KIND=8) :: x1,x2
!    iset=0
    if (iset.eq.0.or.force_iset) then
1      CALL wrap_one_rnd(x1)
       CALL wrap_one_rnd(x2)
!       call random_seed(get=aseed)
!       print *,'rnd_call+inicount=',rnd_callcount,rnd_initcount
!       do i=1,12
!          print *,'GAUSS_DEVIATE: i,aseed =',i,aseed(i)
!       enddo
! export
       x1_gauss=x1
       x2_gauss=x2
!       
       v1=2.d0*x1-1.d0
       v2=2.d0*x2-1.d0
       rsq=v1**2+v2**2
       if (rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
       fac=sqrt(-2.d0*log(rsq)/rsq)
       gset=v1*fac
       xgauss=v2*fac
       iset=1
    else
       xgauss=gset
       iset=0
!       print *,'x1,x2,xgauss =',x1,x2,xgauss
    endif
    return
  END SUBROUTINE gauss_deviate

  SUBROUTINE gauss3d(sigma,v)
!---- 3D Gaussian distribution of one random particle.
!     Becomes Maxwell-Boltzmann for sigma = sqrt(k*T/mu) = sqrt(R*T/M)
    REAL(KIND=8) :: xgauss(3),sigma,v(3)!,vsq(n)
    CALL gauss_deviate(xgauss(1))!; v(1)=xgauss*sigma
    CALL gauss_deviate(xgauss(2))!; v(2)=xgauss*sigma
    CALL gauss_deviate(xgauss(3))!; v(3)=xgauss*sigma
    v(:)=xgauss(:)*sigma
    RETURN
  END SUBROUTINE gauss3d

  SUBROUTINE gauss_truncated(maxsigma,xgauss)
    REAL(KIND=8) :: maxsigma,xgauss
    xgauss=2.d0*maxsigma
    IF (maxsigma<=0.d0) THEN!criterion turned off
       CALL gauss_deviate(xgauss)
       RETURN
    ENDIF
    DO WHILE (xgauss>maxsigma)
       CALL gauss_deviate(xgauss)
    END DO
    RETURN
  END SUBROUTINE gauss_truncated

  SUBROUTINE gauss3d_truncated(sigma,maxsigma,v)
!---- 3D Gaussian distribution of one random particle.
!     Becomes Maxwell-Boltzmann for sigma = sqrt(k*T/mu) = sqrt(R*T/M)
    REAL(KIND=8) :: xgauss(3),sigma,maxsigma,v(3)!,vsq(n)
    CALL gauss_truncated(maxsigma,xgauss(1))!; v(1)=xgauss*sigma
    CALL gauss_truncated(maxsigma,xgauss(2))!; v(2)=xgauss*sigma
    CALL gauss_truncated(maxsigma,xgauss(3))!; v(3)=xgauss*sigma
    v(:)=xgauss(:)*sigma
    RETURN
  END SUBROUTINE gauss3d_truncated

  SUBROUTINE lorentz_deviate(xlorentz)
    REAL(KIND=8) :: xlorentz,xgauss1,xgauss2
1   CALL gauss_deviate(xgauss1)
    IF (xold<=0.d0) THEN
       CALL gauss_deviate(xgauss2)
    ELSE
       xgauss2=xold
    ENDIF
    xold=xgauss1!recycle
    IF (xgauss2==0.d0) GOTO 1
! The ratio of two independent Gaussian deviates is Lorentz-distributed
    xlorentz=xgauss1/xgauss2
    RETURN
  END SUBROUTINE lorentz_deviate

  SUBROUTINE logistic_deviate(xlogistic)
    REAL(KIND=8) :: u,sigma,xlogistic
    REAL(KIND=8) :: x
    CALL wrap_one_rnd(x)
    xlogistic=log(x/(1.-x))
    RETURN
  END SUBROUTINE logistic_deviate

  SUBROUTINE sech_deviate(xsech)
    REAL(KIND=8) :: x,xsech
    REAL(KIND=8) :: pi=3.141592653589793238D0
    CALL wrap_one_rnd(x)
    xsech=2.d0/pi*log(tan(pi/2.d0*x))
    RETURN
  END SUBROUTINE sech_deviate

!**** Distributions with additional parameters ****

  SUBROUTINE voigt_deviate(sigma,gamma,u,v,xgauss,xlorentz,xvoigt)
! A Voigt profile is a convolution of Gaussian with Lorentzian profile
    REAL(KIND=8) :: xgauss,xlorentz,xvoigt,sigma,gamma,u,v
    CALL gauss_deviate(xgauss)
    CALL lorentz_deviate(xlorentz)
    xgauss=xgauss*sigma+u
    xlorentz=xlorentz*gamma+v
    xvoigt=xgauss+xlorentz
    RETURN
  END SUBROUTINE voigt_deviate

  SUBROUTINE weibull_distribution(lambda,k,xweibull)
    real(kind=8) :: xflat,xweibull,lambda,k
    call wrap_one_rnd(xflat)
    xweibull=1./lambda*(-log(1.d0-xflat))**(1.d0/k)
    return
  END SUBROUTINE weibull_distribution

END MODULE deviates_module

! S(q)uare-s(q)uare (M)iller 3D
! 3D random fluctuation generator based on the square-square alorithm
! by Gavin S. P. Miller: The definition and rendering of terrain maps
! ACM SIGGRAPH Computer Graphics, Band 20, Nr. 4, 1986, S. 39–48
! Extended to 3D by Ingo Thies in 2017.
module qqm3d_module
  use deviates_module
  implicit none
  integer :: nsize,fchar,qqm_nx,qqm_ny,qqm_nz,init_mode=0,deviate_mode=0,maxiter=100
  integer :: qqm_xmin,qqm_xmax,qqm_ymin,qqm_ymax,qqm_zmin,qqm_zmax,qqm_allhmode,qqm_balance
  real(kind=8),dimension(:,:,:),allocatable :: afield,bfield
  real(kind=8),dimension(2,2,2) :: aaux,baux
  real(kind=8),dimension(0:1) :: qqm_xrange,qqm_yrange,qqm_zrange!output image ranges
  real(kind=8),dimension(0:12) :: qqm_allscalh,qqm_powerspectrum
  real(kind=8) :: scalh_ini=1.d0,scalh,rfchar,hreduce=2.d0,hred_min,hred_max,hred_fac=1.d0,hred_pwr=1.d0
  real(kind=8) :: qqm_valmin,qqm_valmax,qqm_valavg,qqm_valmed,qqm_fracpos,qqm_valcorr
  real(kind=8),dimension(1:27) :: ax_ini
  real(kind=8) :: qqm_x0,qqm_xn,deviate_power=1.d0,x0_ini=0.d0,vlevel(3),vtol(3)
!  real(kind=8) :: qqm_k1=9.d0,qqm_k2=4.d0,qqm_k3=2.d0,qqm_k4=1.d0, qqm_k0=1.d0
  real(kind=8) :: qqm_k1=27.d0,qqm_k2=9.d0,qqm_k3=3.d0,qqm_k4=1.d0, qqm_k0=0.d0
  integer :: qqm_output_mode=0,selz,selslice
  integer :: qqm_verbose=0,qqm_debug=0,qqm_seed,flg_rndhred=0,maxit_hrnd,qqm_nb,qqm_border_mode
  logical :: didallocate=.false.
contains

  subroutine drive_qqm3d
    implicit none
    integer :: seed,seedmode,i,j,k,iter,irange,n1,n2,n3,hr,icount
    integer :: iter0
    real(kind=8) :: xdum
    call init_afield
    call init_corners
    iter0=fchar
    icount=0
    do iter=1,maxiter
       n1=2**(iter-1)+1
       n2=n1
       n3=n1
       if (qqm_allhmode>=1) then
          if (maxit_hrnd>0.and.iter>=maxit_hrnd.and.iter>iter0+1.and.qqm_allhmode>=9) then
             scalh=0.d0
          else
             scalh=qqm_allscalh(iter)
          endif
       else
          if (iter<=iter0) then
             scalh=0.d0
          else if (iter==iter0+1) then
             scalh=scalh_ini
          else if (maxit_hrnd>0.and.iter>=maxit_hrnd) then
             scalh=0.d0
          endif
       endif
       qqm_powerspectrum(iter)=scalh
       do i=1,n1
          do j=1,n2
             do k=1,n3
                icount=icount+1
                call qqmiller3d(i,j,k)
             enddo
          enddo
       enddo
       afield=bfield
       ! flg_rndreduce=1: randomize hreduce once per iteration
       if (flg_rndhred==1.or.flg_rndhred==3) then
          call uniform_range(hreduce,hred_min,hred_max)
       else if (flg_rndhred==-1.or.flg_rndhred==-3) then
          call wrap_one_rnd(xdum)!dummy call to keep seed sequence
       endif
       if (qqm_output_mode>=-1) print *,'iter,scalh=',iter,scalh
       if (iter>iter0.and.qqm_allhmode==0) scalh=scalh/hreduce
       if (2*n1>=qqm_nx) exit
    enddo
    if (qqm_balance>=1) then
       print *,'Balancing values...'
       xdum=0.d0
       do k=1,qqm_nz
          do j=1,qqm_ny
             do i=1,qqm_nx
                xdum=xdum+afield(i,j,k)
             enddo
          enddo
       enddo
       qqm_valcorr=-xdum/(qqm_nx*qqm_ny*qqm_nz)
       print *,'Correcting values by ',qqm_valcorr
       do k=1,qqm_nz
          do j=1,qqm_ny
             do i=1,qqm_nx
                afield(i,j,k)=afield(i,j,k)+qqm_valcorr
             enddo
          enddo
       enddo
    endif
    print *,'Array creation finished'
    return
  end subroutine drive_qqm3d

  subroutine write_qqm3d(iu1,iu2,iu3,iu4,iu5,iu6)
    integer,parameter :: nbmax=10000
    integer :: imod,bmod,nb,iybin(nbmax),ib,icount,flgmedian
    integer :: ix0,ix1,iy0,iy1,iz0,iz1
    integer :: iu1,iu2,iu3,iu4,iu5,iu6,i,j,slice,ii,jj,kslice
    integer :: imin,jmin,kmin,imax,jmax,kmax
    real(kind=8) :: x,x0,xn,ybin(3,nbmax),kw=1.d0
    real(kind=8) :: val,cumval,ss,ss1,scalhist,scalcum
    character(len=30) :: outfile
    if (qqm_border_mode==1) then
       qqm_xrange=(/0.d0,1.d0/)
       qqm_yrange=(/0.d0,1.d0/)
       qqm_zrange=(/0.d0,1.d0/)
    else
       qqm_xrange=(/0.5d0,dble(qqm_nx)+0.5d0/)
       qqm_yrange=(/0.5d0,dble(qqm_ny)+0.5d0/)
       qqm_zrange=(/0.5d0,dble(qqm_nz)+0.5d0/)
    endif
    if (qqm_output_mode>=1) then
       open(iu4,file='xqq3d-level1.out')
       open(iu5,file='xqq3d-level2.out')
       open(iu6,file='xqq3d-level3.out')
       ! dummy outputs
       write(iu4,*) '# vlevel,vtol=',vlevel(1),vtol(1)
       write(iu5,*) '# vlevel,vtol=',vlevel(2),vtol(2)
       write(iu6,*) '# vlevel,vtol=',vlevel(3),vtol(3)
    endif
    if (selslice>0.and.selslice<=qqm_nz) then
       iz0=selslice
       iz1=iz0
    else
       iz0=1
       iz1=qqm_nz
    endif
    ss=0.d0
    ss1=0.d0
    qqm_valmin=1.d99
    qqm_valmax=-1.d99
    icount=0
    imod=0
    bmod=0
!    print *,'selslice,iz0,iz1 =',selslice,iz0,iz1
    do slice=iz0,iz1
       if (qqm_output_mode>=2) write(outfile,'(a,i4.4,a)') 'xqq3d.',slice,'.out'
       if (qqm_output_mode>=2) print '(3a)','outfile = "',trim(outfile),'"'
       if (qqm_output_mode>=2) open(iu1,file=trim(outfile))
       do j=1,qqm_ny
          do i=1,qqm_nx
             icount=icount+1
             if (selz==1) then
                val=afield(slice,j,i)
             else if (selz==2) then
                val=afield(i,slice,j)
             else
                val=afield(i,j,slice)
             endif
             if (val<qqm_valmin) then
                qqm_valmin=min(qqm_valmin,val)
                imin=i
                jmin=j
                kmin=slice
             endif
             if (val>qqm_valmax) then
                qqm_valmax=max(qqm_valmax,val)
                imax=i
                jmax=j
                kmax=slice
             endif
             ss=ss+val
             if (val>0.d0) ss1=ss1+1.d0
             if (qqm_output_mode>=-1) call xbin(imod,bmod,val,qqm_x0,qqm_xn,qqm_nb,kw,ybin,iybin,ib)
             if (qqm_border_mode==1) then
                if (qqm_output_mode>=2) write(iu1,11) (dble(i)-0.5d0)/dble(qqm_nx),(dble(j)-0.5d0)/dble(qqm_ny),val
             else
                if (qqm_output_mode>=2) write(iu1,10) i,j,val
             endif
             if (qqm_output_mode>=1) then
                if (abs(val-vlevel(1))<=vtol(1)) then
                   if (qqm_border_mode==1) then
                      write(iu4,21) (dble(i)-0.5d0)/dble(qqm_nx),(dble(j)-0.5d0)/dble(qqm_ny),(dble(slice)-0.5d0)/dble(qqm_nz)
                   else
                      write(iu4,20) i,j,slice
                   endif
                else if (abs(val-vlevel(2))<=vtol(2)) then
                   if (qqm_border_mode==1) then
                      write(iu5,21) (dble(i)-0.5d0)/dble(qqm_nx),(dble(j)-0.5d0)/dble(qqm_ny),(dble(slice)-0.5d0)/dble(qqm_nz)
                   else
                      write(iu5,20) i,j,slice
                   endif
                else if (abs(val-vlevel(3))<=vtol(3)) then
                   if (qqm_border_mode==1) then
                      write(iu6,21) (dble(i)-0.5d0)/dble(qqm_nx),(dble(j)-0.5d0)/dble(qqm_ny),(dble(slice)-0.5d0)/dble(qqm_nz)
                   else
                      write(iu6,20) i,j,slice
                   endif
                endif
             endif
             imod=2
          enddo
          if (qqm_output_mode>=2.and.j<qqm_ny) write(iu1,*)
       enddo
       if (qqm_output_mode>=2) close(iu1)
    enddo
    qqm_valavg=ss/icount
    qqm_fracpos=ss1/icount
    if (qqm_output_mode>=0) open(iu3,file='xqq3d-hist.out')
    flgmedian=0
    scalcum=1.d0/icount
    scalhist=dble(qqm_nb)/(qqm_xn-qqm_x0)*scalcum
    do i=1,qqm_nb
       x=ybin(1,i)
       val=ybin(2,i)*scalhist
       cumval=ybin(3,i)*scalcum
       if (flgmedian==0.and.i>1.and.ybin(3,i)>=0.5d0*ybin(3,qqm_nb)) then
          if (ybin(3,i)==0.5d0*ybin(3,qqm_nb)) then
             qqm_valmed=ybin(1,i)
          else
             qqm_valmed=0.5d0*(ybin(1,i)+ybin(1,i-1))
          endif
          flgmedian=1
       endif
       if (qqm_output_mode>=0) write(iu3,40) x,val!,sqrt(ybin(2,i))*scalhist
    enddo
    if (qqm_output_mode>=0) close(iu3)
    if (qqm_output_mode>=-1) write(*,30) qqm_valmin,imin,jmin,kmin,&
         qqm_valmax,imax,jmax,kmax,qqm_valavg,qqm_valmed,&
         1.d2*qqm_fracpos,qqm_seed
    if (qqm_output_mode>=0) then
       open(iu2,file='xqq3d-stats.out')
       write(iu2,30) qqm_valmin,imin,jmin,kmin,qqm_valmax,imax,jmax,kmax,&
            qqm_valavg,qqm_valmed,1.d2*qqm_fracpos,qqm_seed
       close(iu2)
    endif
    if (qqm_output_mode>=1) then
       close(iu4)
       close(iu5)
       close(iu6)
    endif
10  format(i4.4,1x,i4.4,1x,e16.8)
11  format(f8.6,1x,f8.6,1x,e16.8)
20  format(i5.4,2(1x,i5.4))
21  format(f16.8,2(1x,f16.8))
30  format('Min. value =',F16.8,' at i,j,k =',3(1x,i4.4)/&
           'Max. value =',F16.8,' at i,j,k =',3(1x,i4.4)/&
           'Avg. value =',F16.8/&
           'Median     =',F16.8/&
           'Pos. fract.=',F6.2,'%'/&
           'Random seed=',I10)
40  format(f12.6,2(1x,1pe14.6e3))
    return
  end subroutine write_qqm3d

  subroutine qqmiller3d(ix,iy,iz)
    implicit none
    integer :: ix,iy,iz,i1,j1,k1,i2,j2,k2,hr,qr
    aaux(1,1,1)=afield(ix,iy,iz)
    aaux(2,1,1)=afield(ix+1,iy,iz)
    aaux(1,2,1)=afield(ix,iy+1,iz)
    aaux(2,2,1)=afield(ix+1,iy+1,iz)
    aaux(1,1,2)=afield(ix,iy,iz+1)
    aaux(2,1,2)=afield(ix+1,iy,iz+1)
    aaux(1,2,2)=afield(ix,iy+1,iz+1)
    aaux(2,2,2)=afield(ix+1,iy+1,iz+1)
    call cube1_step
    call cube2_step
    call cube3_step
    call cube4_step
    call cube5_step
    call cube6_step
    call cube7_step
    call cube8_step
    i1=2*(ix-1)+1
    j1=2*(iy-1)+1
    k1=2*(iz-1)+1
    i2=2*(ix-1)+2
    j2=2*(iy-1)+2
    k2=2*(iz-1)+2
!    print *,'i1,i2,j1,j2,k1,k2 =',i1,i2,j1,j2,k1,k2
    bfield(i1,j1,k1)=baux(1,1,1)
    bfield(i2,j1,k1)=baux(2,1,1)
    bfield(i1,j2,k1)=baux(1,2,1)
    bfield(i2,j2,k1)=baux(2,2,1)
    bfield(i1,j1,k2)=baux(1,1,2)
    bfield(i2,j1,k2)=baux(2,1,2)
    bfield(i1,j2,k2)=baux(1,2,2)
    bfield(i2,j2,k2)=baux(2,2,2)
    return
  end subroutine qqmiller3d
  
  subroutine init_corners
    real(kind=8),dimension(1:27) :: xx
    integer :: i
    ! set initial corner values
    do i=1,27
       if (init_mode==1) then
          call get_rnd_dh(xx(i))
          xx(i)=qqm_allscalh(0)*xx(i)+x0_ini
       else if (init_mode==2) then
          call get_rnd_dh(xx(i))
          xx(i)=qqm_allscalh(0)*xx(i)+ax_ini(i)+x0_ini
       else if (init_mode==11) then
          call get_rnd_dh(xx(i))!dummy call for seed consistency
          xx(i)=x0_ini
       else if (init_mode==12) then
          call get_rnd_dh(xx(i))!dummy call for seed consistency
          xx(i)=ax_ini(i)+x0_ini
       else if (init_mode==-1) then
          xx(i)=x0_ini
       else
          xx(i)=ax_ini(i)+x0_ini
       endif
    enddo
    !slice 1
    afield(1,1,1)=xx(1)
    afield(2,1,1)=xx(2)
    afield(3,1,1)=xx(3)
    afield(1,2,1)=xx(4)
    afield(2,2,1)=xx(5)
    afield(3,2,1)=xx(6)
    afield(1,3,1)=xx(7)
    afield(2,3,1)=xx(8)
    afield(3,3,1)=xx(9)
    !slice 2
    afield(1,1,2)=xx(10)
    afield(2,1,2)=xx(11)
    afield(3,1,2)=xx(12)
    afield(1,2,2)=xx(13)
    afield(2,2,2)=xx(14)
    afield(3,2,2)=xx(15)
    afield(1,3,2)=xx(16)
    afield(2,3,2)=xx(17)
    afield(3,3,2)=xx(18)
    !slice 3
    afield(1,1,3)=xx(19)
    afield(2,1,3)=xx(20)
    afield(3,1,3)=xx(21)
    afield(1,2,3)=xx(22)
    afield(2,2,3)=xx(23)
    afield(3,2,3)=xx(24)
    afield(1,3,3)=xx(25)
    afield(2,3,3)=xx(26)
    afield(3,3,3)=xx(27)
  end subroutine init_corners
  
  subroutine init_afield
    implicit none
    if (nsize>12.or.nsize<2) then
       print *,'nsize must be between 2 and 12'
       stop
    endif
    qqm_nx=2**nsize+2
    qqm_ny=qqm_nx
    qqm_nz=qqm_nx
    if (.not.didallocate) then
       allocate(afield(qqm_nx,qqm_ny,qqm_nz))
       allocate(bfield(qqm_nx,qqm_ny,qqm_nz))
       didallocate=.true.
    endif
    afield=0.d0
    bfield=0.d0
    print *,'qqm_nx =',qqm_nx
    return
  end subroutine init_afield

  subroutine close_afield
    if (didallocate) then
       deallocate(afield)
       deallocate(bfield)
       didallocate=.false.
    endif
    return
  end subroutine close_afield

  subroutine cube1_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k1
    fb=qqm_k2
    fc=qqm_k2
    fd=qqm_k3
    fe=qqm_k2
    ff=qqm_k3
    fg=qqm_k3
    fh=qqm_k4
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(1,1,1)=val
    return
  end subroutine cube1_step
  
  subroutine cube2_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k2
    fb=qqm_k1
    fc=qqm_k3
    fd=qqm_k2
    fe=qqm_k3
    ff=qqm_k2
    fg=qqm_k4
    fh=qqm_k3
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(2,1,1)=val
    return
  end subroutine cube2_step
  
  subroutine cube3_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k2
    fb=qqm_k3
    fc=qqm_k1
    fd=qqm_k2
    fe=qqm_k3
    ff=qqm_k4
    fg=qqm_k2
    fh=qqm_k3
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(1,2,1)=val
    return
  end subroutine cube3_step
  
  subroutine cube4_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k3
    fb=qqm_k2
    fc=qqm_k2
    fd=qqm_k1
    fe=qqm_k4
    ff=qqm_k3
    fg=qqm_k3
    fh=qqm_k2
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(2,2,1)=val
    return
  end subroutine cube4_step
  
  subroutine cube5_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k2
    fb=qqm_k3
    fc=qqm_k3
    fd=qqm_k4
    fe=qqm_k1
    ff=qqm_k2
    fg=qqm_k2
    fh=qqm_k3
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(1,1,2)=val
    return
  end subroutine cube5_step
  
  subroutine cube6_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k3
    fb=qqm_k2
    fc=qqm_k4
    fd=qqm_k3
    fe=qqm_k2
    ff=qqm_k1
    fg=qqm_k3
    fh=qqm_k2
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(2,1,2)=val
    return
  end subroutine cube6_step
  
  subroutine cube7_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k3
    fb=qqm_k4
    fc=qqm_k2
    fd=qqm_k3
    fe=qqm_k2
    ff=qqm_k3
    fg=qqm_k1
    fh=qqm_k2
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(1,2,2)=val
    return
  end subroutine cube7_step
  
  subroutine cube8_step
    implicit none
    integer :: irange,icount
    integer :: ax,ay,bx,by,cx,cy,dx,dy
    real(kind=8) ::  avg,val,x,div,fa,fb,fc,fd,fe,ff,fg,fh,xdum,xr
    fa=qqm_k4
    fb=qqm_k3
    fc=qqm_k3
    fd=qqm_k2
    fe=qqm_k3
    ff=qqm_k2
    fg=qqm_k2
    fh=qqm_k1
    div=fa+fb+fc+fd+fe+ff+fg+fh
    avg= (fa*aaux(1,1,1)+fb*aaux(2,1,1)+fc*aaux(1,2,1)+fd*aaux(2,2,1)&
         +fe*aaux(1,1,2)+ff*aaux(2,1,2)+fg*aaux(1,2,2)+fh*aaux(2,2,2))/div
    call get_rnd_dh(x)
    xr=1
    if (flg_rndhred==2) then
       call uniform_range(xr,hred_min/2.d0,hred_max/2.d0)
    else if (flg_rndhred==-2) then
       call wrap_one_rnd(xdum)!dummy call to keep seed sequence
    endif
    val=avg+xr*scalh*x
    if (qqm_debug==1) val=avg
    baux(2,2,2)=val
    return
  end subroutine cube8_step
  
  subroutine get_rnd_dh(xdh)
    implicit none
    real(kind=8) :: x,xdh,pwr
    if (deviate_mode==1) then
       call gauss_deviate(x)
    else if (deviate_mode==2) then
       call logistic_deviate(x)
    else if (deviate_mode==3) then
       call sech_deviate(x)
    else if (deviate_mode==-1) then
       x=0.d0!debug mode
    else
       call uniform_pm1(x)
    endif
    pwr=deviate_power
    call sgnpwr(x,xdh,pwr)
    return
  end subroutine get_rnd_dh

  subroutine sgnpwr(xin,xout,pwr)
    real(kind=8) :: xin,xout,pwr
    xout=sign(abs(xin)**pwr,xin)
    return
  end subroutine sgnpwr
end module qqm3d_module

module storage_module
  use qqm3d_module, only : qqm_nx,qqm_ny,qqm_nz,afield
  real(kind=8),dimension(:,:,:),allocatable :: afield0,afield1,afield2,afield3
  logical :: alloc_0=.false.,alloc_1=.false.,alloc_2=.false.,alloc_3=.false.
  logical :: store_0=.false.,store_1=.false.,store_2=.false.,store_3=.false.
contains

  subroutine init_storage(ia)
    integer :: ia
    if (ia==0.and..not.alloc_0) then
       allocate(afield0(qqm_nx,qqm_ny,qqm_nz))
       alloc_0=.true.
    endif
    if (ia==1.and..not.alloc_1) then
       allocate(afield1(qqm_nx,qqm_ny,qqm_nz))
       alloc_1=.true.
    endif
    if (ia==2.and..not.alloc_2) then
       allocate(afield2(qqm_nx,qqm_ny,qqm_nz))
       alloc_2=.true.
    endif
    if (ia==3.and..not.alloc_3) then
       allocate(afield3(qqm_nx,qqm_ny,qqm_nz))
       alloc_3=.true.
    endif
    if (ia==0.and..not.store_0) then
       afield0=afield
       store_0=.true.
    endif
    if (ia==1.and..not.store_1) then
       afield1=afield
       store_1=.true.
    endif
    if (ia==2.and..not.store_2) then
       afield2=afield
       store_2=.true.
    endif
    if (ia==3.and..not.store_3) then
       afield3=afield
       store_3=.true.
    endif
    return
  end subroutine init_storage

  subroutine close_storage
    if (alloc_0) then
       deallocate(afield0)
       alloc_0=.false.
       store_0=.false.
    endif 
    if (alloc_1) then
       deallocate(afield1)
       alloc_1=.false.
       store_1=.false.
    endif
    if (alloc_2) then
       deallocate(afield2)
       alloc_2=.false.
       store_2=.false.
    endif
    if (alloc_3) then
       deallocate(afield3)
       alloc_3=.false.
       store_3=.false.
    endif
  end subroutine close_storage
end module storage_module

SUBROUTINE xbin(imod,bmod,x,x0,xn,nb,kw,ybin,iybin,ib)
!******** BINNING ALGORITHM v 0.50 ****
!     XBIN Calculates the number of events within interval partition
!     [u_i,u_i+1[, where u_i = x0 + i*(xn-x0)/nb
!     for i=0...nb
!
!     Example for 4 bins and x0=-10, x1=10
!
!     bin#  | 1 | 2 | 3 | 4 |
!     x    -10 -5   0   5   10
!     u_i   0   1   2   3   4
!
!     Bins are labelled as follows:
!     Bin #j := [x_j-db/2, x_j+db/2[ where x_j = (u_{j-1}+u_j)/2
!     and db the bin width. u_j are not explicitely used in the routine.
!
!     OPTIONS
!     imod = setup option: 0 - setup bin positions ybin(1,*) (first each)
!                              (will then be set to 1 automatically)
!                         >0 - re-use from given ybin
!                          2 - same as 1 + cumulative bins
!                         <0 - just get bin index (and setup)
!     bmod = boundary option: 0 - strictly excludes right interval sup
!                             1 - counts event if equal to sup of nth bin
!     INPUT:
!     x    = event abscissa
!     x0,xn= overall interval (from left side of first bin to right side of last bin)
!     nb   = number of partitions (set negative for re-initialization)
!     kw   = weighting factor (normally equal 1)
!
!     OUTPUT:
!     ybin = output array
!            (1,i) = ith midpoint abscissa, i=1...nb
!            (2,i) = ith bin value
!            (3,i) = spare array (used for cumulative binning in imod=2)
!     iybin = unweighted integer binning array
!             bins are advanced by 1 if kw != 0, else unchanged
!     ib    = index of actual bin
!             Boundary: ib=0 for x<u_0, ib=nb+1 for x>(=) u_nb
!
!     AUXILIARY VARIABLES/SWITCHES:
!     Copy corresponding COMMON block into calling program
!     dokb  = ycnt weighting switch    (input)
!                      true: kb=kw*db
!                     false: kb=kw
!     ycnt  = (0) weighted (by kw AND bin width) total number < x0
!             (1)    "       "     "    within range
!             (2)    "       "     "    >(=) xn
!     iycnt = unweighted count of inside/outside hits (if kw!=0)
!     nkweq0= Number if calls with kw=0 (inside + outside bins)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: nb,imod,bmod,iybin(nb),ib,jb!,ncall
  REAL(KIND=8) :: x,u,x0,xn,db,kw,ybin(3,nb),kb
  LOGICAL add2y
!.... Auxiliary variables: Copy this block into calling program to use
  LOGICAL dokb!used for ycnt weighting only
  INTEGER :: iycnt(0:2),nkweq0
  REAL(KIND=8) :: ycnt(0:2)
  COMMON /xbaux/ ycnt,iycnt,nkweq0,dokb
!.... Debugging switch
  LOGICAL :: dotnk
  COMMON /xbntnk/ dotnk
  IF (x0.ge.xn) THEN
     WRITE(*,*) 'XBIN ERROR: x0 < xn violated! --> Abort'
     WRITE(*,*) 'x0 = ',x0,'          xn = ',xn
     STOP
  ENDIF
!---- Setup
  db=(xn-x0)/nb!calculate this always (e.g. for ybin's of different size)
  IF (imod.eq.0) THEN
     CALL setxb(x0,xn,nb,db, ybin,ycnt,iybin,iycnt, dotnk)
     nkweq0=0
     imod=1!unexpected reset may be fatal
  ENDIF
!---- Switches
  add2y=kw.ne.0.d0
  IF (add2y) nkweq0=nkweq0+1
  IF (dokb) THEN
     kb=kw*db
  ELSE
     kb=kw
  ENDIF
!---- Daily use
  u=x-x0
  ib=max(min(int(u/db+1.d0),nb+1),0)
  IF (bmod.eq.1 .and. x.eq.xn) ib=nb
!.... Return if index-only imod is set
  IF (imod.lt.0) RETURN
!.... Continue else 
  jb=max(min(ib,nb),1)
!      if(x0.eq.0.d0) then
  if (dotnk) then
     write(*,1000) 'x,x0,xn,yb,db,nb,ib =',x,x0,xn,ybin(1,jb),db,nb,ib
  endif
  IF (ib.le.0) THEN
     ycnt(0)=ycnt(0)+kb
     IF (add2y) iycnt(0)=iycnt(0)+1
  ELSE IF (ib.gt.nb) THEN
     ycnt(2)=ycnt(2)+kb
     IF (add2y) iycnt(2)=iycnt(2)+1
  ELSE
     ycnt(1)=ycnt(1)+kb
     IF (add2y) iycnt(1)=iycnt(1)+1
     ybin(2,ib)=ybin(2,ib)+kw
     IF (add2y) iybin(ib)=iybin(ib)+1
  ENDIF
!---- Create cumulative histogram (e.g. for final call of xbin)
  IF (imod.eq.2) THEN
     CALL cumbin(nb,ybin)
  ENDIF
1000 FORMAT(A,1X,5(F9.6,2X),2(2X,I6))
  RETURN
END SUBROUTINE xbin

SUBROUTINE setxb(x0,xn,nb,db, ybin,ycnt,iybin,iycnt, dotnk)
  IMPLICIT NONE
  INTEGER :: nb,i
  INTEGER :: iybin(nb),iycnt(0:2)
  REAL(KIND=8) :: x0,xn,db,ybin(3,nb),ycnt(0:2)
  LOGICAL :: dotnk
!      write(*,*) 'XBIN: Initializing...'
  DO i=1,nb
!--- Calculation of bin midpoints
!     ybin(1,i) abscissa equal to bin midpoint
!     ybin(2,*) use small non-zero value to avoid log crashing
     ybin(1,i)=x0+(dble(i)-0.5d0)*db
     ybin(2,i)=1.d-99
     ybin(3,i)=1.d-99
     iybin(i)=0
     if (dotnk) then
        write(*,*) 'i,xi =',i,ybin(1,i)
     endif
  ENDDO
  ycnt(0)=0.d0
  ycnt(1)=0.d0
  ycnt(2)=0.d0
  iycnt(0)=0
  iycnt(1)=0
  iycnt(2)=0
  RETURN
END SUBROUTINE setxb

SUBROUTINE cumbin(nb,ybin)
  IMPLICIT NONE
  INTEGER :: nb,i
  REAL(KIND=8) :: ybin(3,nb)
  ybin(3,1)=ybin(2,1)
  DO i=2,nb
     ybin(3,i)=ybin(3,i-1)+ybin(2,i)
  ENDDO
  RETURN
END SUBROUTINE cumbin

subroutine probe_intarray3d(nx,ny,nz,xrange,yrange,zrange,cdata,x,y,z,val)
  ! array with integer coordinates from 1--n will be probed as being an array
  ! with real coordinates xyzrange(0:1)
  implicit none
  integer :: nx,ny,nz,ix,iy,iz
  real(kind=8),dimension(nx,ny,nz) :: cdata
  real(kind=8),dimension(0:1) :: xrange,yrange,zrange
  real(kind=8) :: x,y,z,val
  real(kind=8) :: dxtmp,dytmp,dztmp
  real(kind=8) :: dx,dy,dz
  logical :: inxrange,inyrange,inzrange
  dx=(xrange(1)-xrange(0))/nx
  dy=(yrange(1)-yrange(0))/ny
  dz=(zrange(1)-zrange(0))/nz
  inxrange=x>=xrange(0).and.x<=xrange(1)
  inyrange=y>=yrange(0).and.y<=yrange(1)
  inzrange=z>=zrange(0).and.z<=zrange(1)
  !Example: xrange=(/0.,3./)
  !
  !  [...1...|...2...|...3...]
  !  0  0.5  1  1.5  2   2.5 3
  ix=nint((x-xrange(0))/dx+0.5d0)
  iy=nint((y-yrange(0))/dy+0.5d0)
  iz=nint((z-zrange(0))/dz+0.5d0)
!  print *,'ix,iy,iz raw  =',ix,iy,iz
  ix=max(min(ix,nx),1)
  iy=max(min(iy,ny),1)
  iz=max(min(iz,nz),1)
!  print *,'ix,iy,iz corr.=',ix,iy,iz
  if (inxrange.and.inyrange.and.inzrange) then
     val=cdata(ix,iy,iz)
  else
     val=0.d0
  endif
  return
end subroutine probe_intarray3d

subroutine linint_intarray3d(nx,ny,nz,xrange,yrange,zrange,cdata,x,y,z,val)
  ! array with integer coordinates from 1--n will be probed as being an array
  ! with real coordinates xyzrange(0:1)
  implicit none
  integer :: nx,ny,nz,ix,iy,iz,ix0,iy0,iz0,ix1,iy1,iz1,ixa,iya,iza
  real(kind=8),dimension(nx,ny,nz) :: cdata
  real(kind=8),dimension(0:1) :: xrange,yrange,zrange
  real(kind=8) :: x,y,z,val,s1,s2,s3,s4,c,c0
  real(kind=8) :: val000,val001,val010,val011
  real(kind=8) :: val100,val101,val110,val111
  real(kind=8) :: valx00,valx01,valx10,valx11
  real(kind=8) :: valxy0,valxy1
  real(kind=8) :: dxtmp,dytmp,dztmp
  real(kind=8) :: dx,dy,dz
  logical :: inxrange,inyrange,inzrange
  dx=(xrange(1)-xrange(0))/(nx-1)
  dy=(yrange(1)-yrange(0))/(ny-1)
  dz=(zrange(1)-zrange(0))/(nz-1)
  inxrange=x>=xrange(0).and.x<=xrange(1)
  inyrange=y>=yrange(0).and.y<=yrange(1)
  inzrange=z>=zrange(0).and.z<=zrange(1)
  !Example: xrange=(/0.,3./)
  !
  !  [...1...|...2...|...3...]
  !  0  0.5  1  1.5  2   2.5 3
  !  
  ix0=floor((x-xrange(0))/dx)+1
  iy0=floor((y-yrange(0))/dy)+1
  iz0=floor((z-zrange(0))/dz)+1
  ix1=floor((x-xrange(0))/dx)+2
  iy1=floor((y-yrange(0))/dy)+2
  iz1=floor((z-zrange(0))/dz)+2
!  print *,'i0x,iy0,iz0 raw  =',ix0,iy0,iz0
!  print *,'i1x,iy1,iz1 raw  =',ix1,iy1,iz1
  ix0=max(min(ix0,nx),1)
  iy0=max(min(iy0,ny),1)
  iz0=max(min(iz0,nz),1)
  ix1=max(min(ix1,nx),1)
  iy1=max(min(iy1,ny),1)
  iz1=max(min(iz1,nz),1)
  ixa=ix0; iya=iy0; iza=iz0
  if (ix1==nx) ixa=ix1-1
  if (iy1==ny) iya=iy1-1
  if (iz1==nz) iza=iz1-1
!  print *,'iax,iya,iza corr  =',ixa,iya,iza
!  print *,'i0x,iy0,iz0 corr  =',ix0,iy0,iz0
!  print *,'i1x,iy1,iz1 corr  =',ix1,iy1,iz1
  if (inxrange.and.inyrange.and.inzrange) then
     val000=cdata(ixa,iya,iza)
     val001=cdata(ixa,iya,iz1)
     val010=cdata(ixa,iy1,iza)
     val011=cdata(ixa,iy1,iz1)
     val100=cdata(ix1,iya,iza)
     val101=cdata(ix1,iya,iz1)
     val110=cdata(ix1,iy1,iza)
     val111=cdata(ix1,iy1,iz1)
     c0=xrange(0)+(ixa-1)*dx
     c=x-c0
     s1=(cdata(ix1,iya,iza)-cdata(ixa,iya,iza))/dx
     s2=(cdata(ix1,iya,iz1)-cdata(ixa,iya,iz1))/dx
     s3=(cdata(ix1,iy1,iza)-cdata(ixa,iy1,iza))/dx
     s4=(cdata(ix1,iy1,iz1)-cdata(ixa,iy1,iz1))/dx
     valx00=c*s1+cdata(ixa,iya,iza)
     valx01=c*s2+cdata(ixa,iya,iz1)
     valx10=c*s3+cdata(ixa,iy1,iza)
     valx11=c*s4+cdata(ixa,iy1,iz1)
     c0=yrange(0)+(iya-1)*dy
     c=y-c0
     s1=(valx10-valx00)/dy
     s2=(valx11-valx01)/dy
     valxy0=c*s1+valx00
     valxy1=c*s2+valx01
     c0=zrange(0)+(iza-1)*dz
     c=z-c0
     s1=(valxy1-valxy0)/dz
     val=c*s1+valxy0
!     print *,'c,z,c0 =',c,z,c0
  else
     val=0.d0
  endif
  return
end subroutine linint_intarray3d

subroutine scale_perturbation(mode,nx,ny,nz,xrange,yrange,zrange,cdata,&
     x,y,z,vflat,vscale,val)
  implicit none
  integer :: mode,nx,ny,nz
  real(kind=8),dimension(nx,ny,nz) :: cdata
  real(kind=8),dimension(0:1) :: xrange,yrange,zrange
  real(kind=8) :: vflat,vscale,x,y,z,val,val0
  if (mode>=1) then
     call linint_intarray3d(nx,ny,nz,xrange,yrange,zrange,cdata,x,y,z,val0)
  else
     call probe_intarray3d(nx,ny,nz,xrange,yrange,zrange,cdata,x,y,z,val0)
  endif
  val=val0*vscale+vflat
  return
end subroutine scale_perturbation

subroutine init_qqm3d(iunit,iomode)
  use qqm3d_module
  integer :: iunit,iomode,seed,seedmode,i,count0,hmode
  real(kind=8) :: rdum,kdum,ldum,rfrest
  logical :: found1
  open(iunit,file='qqm3d.par')
  read(iunit,*)
  read(iunit,*) nsize,rfchar,scalh_ini,qqm_balance
  read(iunit,*) init_mode,deviate_mode,deviate_power
  read(iunit,*) x0_ini
  read(iunit,*) hreduce,hred_min,hred_max,flg_rndhred
  read(iunit,*) qqm_allhmode,maxit_hrnd,hred_fac,hred_pwr
  read(iunit,*) seed,seedmode
  read(iunit,*) selz,selslice
  read(iunit,*) vlevel(1),vtol(1)
  read(iunit,*) vlevel(2),vtol(2)
  read(iunit,*) vlevel(3),vtol(3)
  read(iunit,*) qqm_output_mode,qqm_border_mode
  read(iunit,*) 
  read(iunit,*) qqm_x0,qqm_xn,qqm_nb
  count0=0
  found1=.false.
  kdum=1.d0
  ldum=1.d0
  fchar=int(rfchar)
  rfrest=rfchar-fchar
  do i=0,12
     if (i<=1) read(iunit,*)
     read(iunit,*) rdum
     if (rdum<=0.d0.and..not.found1) count0=count0+1
     if (rdum>0.d0.and..not.found1) then
        kdum=hreduce**(1.d0-min(rdum,2.d0))
     endif
     found1=found1.or.rdum>0.d0
     if (qqm_allhmode<=1) then
        qqm_allscalh(i)=scalh_ini*rdum
     else if (qqm_allhmode==2) then
        qqm_allscalh(i)=scalh_ini*rdum/hreduce**(i-1)
     else if (qqm_allhmode==3) then
        qqm_allscalh(i)=scalh_ini*rdum/hreduce**(i-1-count0)
     else if (qqm_allhmode==4) then
        qqm_allscalh(i)=scalh_ini*kdum*rdum/hreduce**(i-1-count0)
     else if (qqm_allhmode==9) then
        if (i<=fchar) then
           qqm_allscalh(i)=0.d0
        else if (i==fchar+1) then
           qqm_allscalh(i)=scalh_ini*(1.d0-rfrest)/hreduce**(i-1)
        else
           qqm_allscalh(i)=scalh_ini/hreduce**(i-1)
        endif
     else if (qqm_allhmode==10) then
        if (i<=fchar) then
           qqm_allscalh(i)=0.d0
        else if (i==fchar+1) then
           ldum=scalh_ini*hreduce**rfrest
           qqm_allscalh(i)=ldum*(1.d0-rfrest)
        else
           qqm_allscalh(i)=ldum/hreduce**(i-1-fchar)
        endif
     endif
!     print *,'i,qqm_allscalh(i) =',i,qqm_allscalh(i)
  enddo
  read(iunit,*)
  do i=1,27
     read(iunit,*) ax_ini(i)
  enddo
  close(iunit)
  if (seedmode==0) then
     seed=0
     print *,'Using randomized seed'
  else
     print *,'Using seed =',seed
  endif
  qqm_seed=seed
  call init_random_seed(seed)
  rnd_ioseed=-1
  if (iomode>=0) call drive_qqm3d
  if (iomode>=1) call write_qqm3d(20,25,30,40,41,42)
  return
end subroutine init_qqm3d

subroutine probe_qqm3d(mode,sel,x,y,z,boxsize,boxcenter,vflat,vscale,val)
  use qqm3d_module
  use storage_module
  integer :: mode,sel
  real(kind=8) :: x,y,z,boxsize,vflat,vscale,val
  real(kind=8),dimension(3) :: boxcenter
  real(kind=8),dimension(0:1) :: xrange,yrange,zrange
  xrange=(/boxcenter(1)-0.5d0*boxsize,boxcenter(1)+0.5d0*boxsize/)
  yrange=(/boxcenter(2)-0.5d0*boxsize,boxcenter(2)+0.5d0*boxsize/)
  zrange=(/boxcenter(3)-0.5d0*boxsize,boxcenter(3)+0.5d0*boxsize/)
  if (sel==0) then
     call init_storage(0)
     call scale_perturbation(mode,qqm_nx,qqm_ny,qqm_nz,xrange,yrange,zrange,afield0,&
          x,y,z,vflat,vscale,val)
  else if (sel==1) then
     call init_storage(0)
     call scale_perturbation(mode,qqm_nx,qqm_ny,qqm_nz,xrange,yrange,zrange,afield1,&
          x,y,z,vflat,vscale,val)
  else if (sel==2) then
     call init_storage(0)
     call scale_perturbation(mode,qqm_nx,qqm_ny,qqm_nz,xrange,yrange,zrange,afield2,&
          x,y,z,vflat,vscale,val)
  else if (sel==3) then
     call init_storage(0)
     call scale_perturbation(mode,qqm_nx,qqm_ny,qqm_nz,xrange,yrange,zrange,afield3,&
          x,y,z,vflat,vscale,val)
  else
     !main array
     call scale_perturbation(mode,qqm_nx,qqm_ny,qqm_nz,xrange,yrange,zrange,afield,&
          x,y,z,vflat,vscale,val)
  endif
  return
end subroutine probe_qqm3d
