subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  ! From Msun/kpc^3 to g/cm^3
  ! 1 Msun/kpc^3 = (1.98855e+33 g) / (3.08567758e+21 cm)^3 = 6.76838229e-32 g/cm^3
  scale_d = (1.98855D+33 / (3.08567758D+21)**3)
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  ! G = 6.674D-8 cm^3 g^-1 s^-2
  ! G = 4.302D-6 kpc Msun^-1 (km/s)^2
  scale_t = 1.0/sqrt(6.674D-8 * scale_d)   ! 1.85e+18
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  ! 1 kpc = 3.08567758 x 10^21 cm
  scale_l = 3.08567758D+21
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v convert velocity in user units into cm/s
  ! scale_v = 1.0D+5  ! 1km/s = 10^5 cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
