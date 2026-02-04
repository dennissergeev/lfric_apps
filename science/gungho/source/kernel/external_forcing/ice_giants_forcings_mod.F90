!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains forcing terms for use in the Ice Giants kernels.
!>
!> @details Support functions for Ice Giants temperature forcing following
!> Guendelman & Kaspi (2025), DOI: 10.1051/0004-6361/202554015; hereafter: GK25

module ice_giants_forcings_mod

  use constants_mod,               only: pi, r_def
  use external_forcing_config_mod, only: theta_relax_time_scale

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Local parameters
  !-------------------------------------------------------------------------------
  ! Parameters for variations w.r.t. pressure (sigma)
  real(kind=r_def), parameter :: temp_0 = 180.0_r_def  ! K
  real(kind=r_def), parameter :: temp_1 = 39.0_r_def  ! K
  real(kind=r_def), parameter :: alpha = 1.2_r_def  ! 1.25 in GK25
  real(kind=r_def), parameter :: beta = 5.0_r_def
  real(kind=r_def), parameter :: zeta = -0.06_r_def  ! 0.06 in GK25 (typo)
  ! Parameters for variations w.r.t. latitude
  real(kind=r_def), parameter :: delta_temp = -4.0_r_def  ! K
  real(kind=r_def), parameter :: width_lat = 20.0_r_def  ! degrees
  real(kind=r_def), parameter :: lat_centre = 35.0_r_def  ! degrees
  real(kind=r_def), parameter :: a_coef = 1.2_r_def  ! 1.5 in GK25
  real(kind=r_def), parameter :: b_coef = -5.0_r_def
  real(kind=r_def), parameter :: c_coef = 0.05_r_def

  public :: ice_giants_newton_frequency
  public :: ice_giants_equilibrium_theta

contains

!> @brief Function to calculate equilibrium theta profile for Ice Giants temperature forcing.
!> @param[in] exner         Exner pressure
!> @param[in] lat           Latitude (in radians)
!> @param[in] sigma         exner/exner0**(1.0/kappa)
!> @param[in] kappa         Ratio of Rd and cp
!> @return    theta_eq      Equilibrium theta
function ice_giants_equilibrium_theta(exner, lat, sigma, kappa) result(theta_eq)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in)    :: exner
  real(kind=r_def), intent(in)    :: sigma
  real(kind=r_def), intent(in)    :: lat, kappa
  ! Intermediate variables
  real(kind=r_def) :: temp_sigma
  real(kind=r_def) :: temp_lat
  real(kind=r_def) :: damping_factor
  ! Output
  real(kind=r_def) :: theta_eq ! Equilibrium theta

  ! Eq. A.1 in GK25
  temp_sigma = ( (temp_0 * sigma ** (alpha * kappa)) ** beta &
                 + (temp_1 * sigma**zeta) ** beta ) ** (1.0_r_def / beta)
  ! Eq. A.2 in GK25
  damping_factor = c_coef + a_coef * exp(b_coef * sigma**0.5_r_def)
  temp_lat = exp(-(((lat * 180.0_r_def / pi - lat_centre) / width_lat) ** 2.0_r_def)) &
           + exp(-(((lat * 180.0_r_def / pi + lat_centre) / width_lat) ** 2.0_r_def))
  ! Eq. A.3 in GK25
  ! lat_centre = 0.0
  ! temp_lat = cos(lat - lat_centre * pi / 180.0_r_def) ** 2.0_r_def

  ! Convert real temperature to potential
  ! temperature using the Exner function (exner*theta=Temp)
  theta_eq = (temp_sigma + delta_temp * temp_lat * damping_factor) / exner

end function ice_giants_equilibrium_theta

!> @brief Function to calculate the Newton relaxation frequency for Ice Giants idealised test case.
function ice_giants_newton_frequency() result(ice_giants_frequency)

  implicit none

  real(kind=r_def) :: ice_giants_frequency

  ! Set as a constant
  ice_giants_frequency = 1.0_r_def / (theta_relax_time_scale * 86400.0_r_def)  ! Table 1 in GK25

end function ice_giants_newton_frequency

end module ice_giants_forcings_mod
