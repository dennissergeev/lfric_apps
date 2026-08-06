!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Evolves the internal_flux field as a Markov process driven by
!>        random spherical harmonics.
!>
!> @details Kernel that updates the prescribed internal flux at the surface
!!          following equations 4-6 of Showman, Tan & Zhang (2019, ApJ,
!!          arXiv:1807.08433). The instantaneous random spatial pattern
!!          (equation 6) is
!!
!!            F = f_amp * sum_{m=1}^{n_f} N_{n_f}^m(sin(phi)) * cos[m(lambda + psi_m)]
!!
!!          where phi and lambda are latitude and longitude, N_n^m are the
!!          normalised associated Legendre polynomials, f_amp is the forcing
!!          amplitude and psi_m is a randomly chosen phase, different for
!!          each zonal wavenumber m and redrawn every call (identically on
!!          every PE) at the algorithm layer, since the paper notes that new
!!          phases are drawn each time equation 4 is evaluated. This
!!          instantaneous pattern is then blended into internal_flux (which
!!          therefore plays the role of the paper's S_h) as a Markov/AR(1)
!!          process, equation 4:
!!
!!            internal_flux(t+dt) = r * internal_flux(t) + sqrt(1-r^2) * F
!!
!!          with the memory coefficient r given by equation 5:
!!
!!            r = 1 - dt / tau_for
!!
!!          where dt is the model timestep and tau_for is the forcing
!!          decorrelation timescale. Both dt and tau_for only enter via the
!!          scalar r, which is computed once at the algorithm layer (it does
!!          not vary from point to point) and passed in here. Persistence of
!!          internal_flux, and hence of this Markov process, across
!!          timesteps and restarts relies entirely on internal_flux being a
!!          checkpointed prognostic field; no separate state is kept.
!!
!!          The Legendre polynomials are evaluated locally at each point using
!!          the same stable recurrence relations as get_Pnm_star_kernel_mod,
!!          rather than via a precomputed field, since only a single row
!!          (n = n_f) of the (small) triangular matrix is required.
module random_harmonics_forcing_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_INTEGER, GH_REAL,       &
                              GH_READ, GH_READWRITE,     &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, pi
use kernel_mod,        only : kernel_type

implicit none

private

public :: random_harmonics_forcing_kernel_type
public :: random_harmonics_forcing_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer. Note that, as for
! spectral_2_cs_kernel_mod, the per-mode random_phase array is not
! representable in PSyclone kernel metadata (GH_ARRAY is not yet supported,
! see PSyclone ticket #1312) and is instead passed directly to
! random_harmonics_forcing_code by a hand-written "psykal-lite" invoker,
! invoke_random_harmonics_forcing_kernel_type in psykal_lite_phys_mod.
type, extends(kernel_type) :: random_harmonics_forcing_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                     &
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! internal_flux
    arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! latitude
    arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! longitude
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! forcing_wavenumber
    arg_type(GH_SCALAR, GH_REAL, GH_READ),                                & ! forcing_amplitude
    arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 & ! ar1_memory_coeff
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: random_harmonics_forcing_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers            Number of layers
!> @param[in,out] internal_flux      Internal flux at the surface (the paper's S_h)
!> @param[in]     latitude           Latitude of each point (radians)
!> @param[in]     longitude          Longitude of each point (radians)
!> @param[in]     forcing_wavenumber Characteristic total wavenumber, n_f
!> @param[in]     forcing_amplitude  Forcing amplitude, f_amp
!> @param[in]     ar1_memory_coeff   AR(1) memory coefficient, r = 1 - dt/tau_for
!> @param[in]     random_phase       Random phase psi_m for m = 1, forcing_wavenumber
!> @param[in]     ndf_1              No. of DOFs per cell for internal_flux space
!> @param[in]     undf_1             No. of unique DOFs for internal_flux space
!> @param[in]     map_1              Dofmap for cell at base of column for internal_flux space
!> @param[in]     ndf_2              No. of DOFs per cell for latitude space
!> @param[in]     undf_2             No. of unique DOFs for latitude space
!> @param[in]     map_2              Dofmap for cell at base of column for latitude space
!> @param[in]     ndf_3              No. of DOFs per cell for longitude space
!> @param[in]     undf_3             No. of unique DOFs for longitude space
!> @param[in]     map_3              Dofmap for cell at base of column for longitude space
subroutine random_harmonics_forcing_code(nlayers,                    &
                                         internal_flux,              &
                                         latitude,                   &
                                         longitude,                  &
                                         forcing_wavenumber,         &
                                         forcing_amplitude,          &
                                         ar1_memory_coeff,           &
                                         random_phase,               &
                                         ndf_1, undf_1, map_1,       &
                                         ndf_2, undf_2, map_2,       &
                                         ndf_3, undf_3, map_3)

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_1, undf_1
  integer(i_def), intent(in) :: map_1(ndf_1)
  integer(i_def), intent(in) :: ndf_2, undf_2
  integer(i_def), intent(in) :: map_2(ndf_2)
  integer(i_def), intent(in) :: ndf_3, undf_3
  integer(i_def), intent(in) :: map_3(ndf_3)

  real(r_def), intent(inout) :: internal_flux(undf_1)
  real(r_def), intent(in)    :: latitude(undf_2)
  real(r_def), intent(in)    :: longitude(undf_3)

  integer(i_def), intent(in) :: forcing_wavenumber
  real(r_def),    intent(in) :: forcing_amplitude
  real(r_def),    intent(in) :: ar1_memory_coeff
  real(r_def),    intent(in) :: random_phase(forcing_wavenumber)

  ! Local variables
  integer(i_def) :: n, m, n_row, n_f_row
  real(r_def)    :: sin_lat, cos_lat
  real(r_def)    :: fact_1st, fact_2nd, coeff
  real(r_def)    :: harmonic_sum

  ! Normalised associated Legendre polynomials (including the spherical
  ! harmonic amplitude factor), stored in the same triangular layout as
  ! get_Pnm_star_kernel_mod: entry (n, m) is at index n*(n+1)/2 + m.
  real(r_def) :: pnm(0:(forcing_wavenumber+1)*(forcing_wavenumber+2)/2 - 1)

  sin_lat = sin(latitude(map_2(1)))
  cos_lat = cos(latitude(map_2(1)))

  ! Pnm(0,0)
  pnm(0) = 1.0_r_def / sqrt(4.0_r_def * pi)

  ! Sectoral terms: Pnm(n+1, n+1) from Pnm(n, n)
  n_row = 0
  do n = 0, forcing_wavenumber - 1
    n_row = n_row + n
    coeff = -sqrt((2.0_r_def*n + 3.0_r_def) / (2.0_r_def*n + 2.0_r_def))
    pnm(n_row + n + 1 + n + 1) = coeff * cos_lat * pnm(n_row + n)
  end do

  ! Near-sectoral terms: Pnm(n+1, n) from Pnm(n, n)
  n_row = 0
  do n = 0, forcing_wavenumber - 1
    n_row = n_row + n
    coeff = sqrt(2.0_r_def*n + 3.0_r_def)
    pnm(n_row + n + 1 + n) = coeff * sin_lat * pnm(n_row + n)
  end do

  ! Remaining terms via the three-term recurrence in n
  n_row = 1
  do n = 2, forcing_wavenumber
    n_row = n_row + n
    do m = 0, n - 2
      fact_1st = sqrt((4.0_r_def*n*n - 1.0_r_def) / real(n*n - m*m, r_def))
      fact_2nd = -sqrt((2.0_r_def*n + 1.0_r_def) / (2.0_r_def*n - 3.0_r_def) * &
                       ((n - 1.0_r_def)**2 - real(m*m, r_def)) /              &
                       real(n*n - m*m, r_def))
      pnm(n_row + m) = fact_1st * sin_lat * pnm(n_row - n + m) + &
                       fact_2nd * pnm(n_row - n - (n - 1) + m)
    end do
  end do

  ! Spherical harmonic normalisation factor for m >= 1
  n_row = 0
  do n = 1, forcing_wavenumber
    n_row = n_row + n
    do m = 1, n
      pnm(n_row + m) = pnm(n_row + m) * sqrt(2.0_r_def)
    end do
  end do

  ! Equation 6: sum over zonal wavenumbers m = 1, n_f at fixed total
  ! wavenumber n = n_f
  n_f_row = forcing_wavenumber * (forcing_wavenumber + 1) / 2

  harmonic_sum = 0.0_r_def
  do m = 1, forcing_wavenumber
    harmonic_sum = harmonic_sum + pnm(n_f_row + m) *                       &
      cos( real(m, r_def) * (longitude(map_3(1)) + random_phase(m)) )
  end do

  ! Equation 4: blend the new random pattern into internal_flux (S_h) as an
  ! AR(1) process, using the memory coefficient r from equation 5. The
  ! sqrt argument is clamped at zero: r can go negative if tau_for is set
  ! shorter than the model timestep, which equation 5 does not itself rule
  ! out, but 1 - r^2 must stay non-negative for the process to be defined.
  internal_flux(map_1(1)) = ar1_memory_coeff * internal_flux(map_1(1)) +   &
    sqrt(max(0.0_r_def, 1.0_r_def - ar1_memory_coeff * ar1_memory_coeff)) * &
    forcing_amplitude * harmonic_sum

end subroutine random_harmonics_forcing_code

end module random_harmonics_forcing_kernel_mod
