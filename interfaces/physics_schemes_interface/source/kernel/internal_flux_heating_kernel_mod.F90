!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies the vertical structure of the random harmonics internal
!>        flux forcing as a potential temperature increment.
!>
!> @details Kernel that builds the 3D forcing S(lambda,phi,p,t) described
!!          just before equation 4 of Showman, Tan & Zhang (2019, ApJ,
!!          arXiv:1807.08433):
!!
!!            S(lambda,phi,p,t) = S_v(p) * S_h(lambda,phi,t)
!!
!!          where S_v(p) is the (nondimensional) vertical structure of the
!!          forcing, described as varying linearly in log-pressure from one
!!          at the base of the domain to zero one scale height above:
!!
!!            S_v(p) = clip[0,1]( 1 - ln(p_base/p) )
!!
!!          Using p = p_zero * exner**(1/kappa), this is evaluated exactly
!!          (no separate scale-height constant needed) as
!!
!!            S_v(k) = clip[0,1]( 1 - (1/kappa) * ln(exner_base/exner(k)) )
!!
!!          with exner_base the Exner pressure at the base of the domain
!!          (level 0) in the same column. S_h(lambda,phi,t) here is the
!!          random_harmonics internal_flux pattern (see
!!          random_harmonics_forcing_kernel_mod), rescaled from its native
!!          units (the surface flux internal_flux_value, W m-2) to the
!!          heating amplitude requested for this 3D pathway (K s-1) via the
!!          scalar heating_ratio = internal_flux_heating_amplitude /
!!          internal_flux_value, computed once at the algorithm layer.
!!          Since both quantities are driven by the identical AR(1) random
!!          realisation (same random phases and memory coefficient), this
!!          rescaling is exact, not an approximation.
module internal_flux_heating_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_READ, GH_WRITE, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

public :: internal_flux_heating_kernel_type
public :: internal_flux_heating_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
type, extends(kernel_type) :: internal_flux_heating_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                   & ! dtheta_internal_flux
    arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                   & ! exner_in_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),& ! internal_flux
    arg_type(GH_SCALAR, GH_REAL, GH_READ),                           & ! kappa
    arg_type(GH_SCALAR, GH_REAL, GH_READ)                            & ! heating_ratio
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: internal_flux_heating_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers              Number of layers
!> @param[in,out] dtheta_internal_flux Potential temperature increment (K s-1)
!> @param[in]     exner_in_wth         Exner pressure on Wtheta levels
!> @param[in]     internal_flux        Random harmonics pattern at the surface (W m-2)
!> @param[in]     kappa                Ratio of gas constant to specific heat, R/cp
!> @param[in]     heating_ratio        internal_flux_heating_amplitude / internal_flux_value
!> @param[in]     ndf_wth              No. of DOFs per cell for Wtheta space
!> @param[in]     undf_wth             No. of unique DOFs for Wtheta space
!> @param[in]     map_wth              Dofmap for cell at base of column for Wtheta space
!> @param[in]     ndf_2d               No. of DOFs per cell for internal_flux space
!> @param[in]     undf_2d              No. of unique DOFs for internal_flux space
!> @param[in]     map_2d               Dofmap for cell at base of column for internal_flux space
subroutine internal_flux_heating_code(nlayers,                     &
                                      dtheta_internal_flux,        &
                                      exner_in_wth,                &
                                      internal_flux,               &
                                      kappa,                       &
                                      heating_ratio,               &
                                      ndf_wth, undf_wth, map_wth,  &
                                      ndf_2d, undf_2d, map_2d)

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: map_wth(ndf_wth)
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)

  real(r_def), intent(inout) :: dtheta_internal_flux(undf_wth)
  real(r_def), intent(in)    :: exner_in_wth(undf_wth)
  real(r_def), intent(in)    :: internal_flux(undf_2d)
  real(r_def), intent(in)    :: kappa
  real(r_def), intent(in)    :: heating_ratio

  ! Local variables
  integer(i_def) :: k
  real(r_def)    :: exner_base, s_v

  ! Exner pressure at the base of the domain (level 0) in this column
  exner_base = exner_in_wth(map_wth(1))

  ! dtheta_internal_flux is assumed to have already been set to zero before
  ! this kernel runs, so levels above the exit point below are left
  ! correctly at zero.
  do k = 0, nlayers

    ! S_v(p) = 1 - ln(p_base/p), using p = p_zero * exner**(1/kappa), so
    ! ln(p_base/p) = (1/kappa) * ln(exner_base/exner(k)). This is exactly 1
    ! at the base of the domain (k = 0, exner = exner_base) and decreases
    ! monotonically with height (exner falls off monotonically with
    ! height), reaching 0 exactly one scale height above the base.
    s_v = 1.0_r_def - &
      (1.0_r_def/kappa) * log(exner_base/exner_in_wth(map_wth(1) + k))

    ! S_v is monotonically decreasing in k, so once it has reached zero it
    ! will stay zero (or negative, before clipping) for all levels above -
    ! nothing more to do.
    if (s_v <= 0.0_r_def) exit

    dtheta_internal_flux(map_wth(1) + k) = &
      heating_ratio * internal_flux(map_2d(1)) * min(1.0_r_def, s_v)

  end do

end subroutine internal_flux_heating_code

end module internal_flux_heating_kernel_mod
