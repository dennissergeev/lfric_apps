!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains forcing terms for use in the CAMEMBERT kernel.
!>
!> @details Support functions for a kernel that adds the CAMEMBERT Case 1
!!          temperature forcing for the sub-Neptunes GJ 1214b and K2-18b, based on
!!          Christie et al. (2022), CAMEMBERT: A Mini-Neptunes General Circulation Model
!!          Intercomparison, Protocol Version 1.0. A CUISINES Model Intercomparison Project,
!!          The Planetary Science Journal, 3, 11, 2022, DOI: 10.3847/PSJ/ac9dfe.

module camembert_forcings_mod

  use constants_mod,               only: r_def, i_def, pi
  use planet_config_mod,           only: kappa, p_zero, cp

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Local parameters (see Table 2 in Christie et al., 2023)
  !---------------------------------------------------------------------------
  ! Longitude of the substellar point
  real(kind=r_def), parameter :: LON_SUBSTELLAR = 0.0_r_def
  ! Equator-to-pole temperature difference
  real(kind=r_def), parameter :: DT_EQ_POLE = 600.0_r_def
  ! Pressure thresholds
  real(kind=r_def), parameter :: P_LOW = 10.0_r_def, P_HIGH = 1.0E6_r_def
  ! Height thresholds (h_1 is higher)
  real(kind=r_def), parameter :: H_1 = 3481339.3_r_def, H_2 = 364935.5_r_def
  ! Universal gas constant (J K-1 mol-1)
  real(kind=r_def), parameter :: GAS_CONST = 8.31446261815324_r_def 

  public :: camembert_equilibrium_theta
  public :: camembert_newton_frequency

contains

!> @brief Calculates equilibrium theta profile for the CAMEMBERT Case 1.
!>
!> @details Calculates the equilibrium (reference) profile of air potential
!!          temperature using eqn. 4 Christie et al. (2022).
!> @param[in] theta_init    Initial potential temperature profile (prescribed)
!> @param[in] exner         Exner function
!> @param[in] layer_height  Height of the kth layer above the surface
!> @param[in] lon           Longitude
!> @param[in] lat           Latitude
!> @return    theta_eq      Equilibrium theta profile
function camembert_equilibrium_theta(theta_init, exner, layer_height, lon, lat) result(theta_eq)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: theta_init   ! Initial potential temperature
  real(kind=r_def), intent(in) :: exner        ! Exner function
  real(kind=r_def), intent(in) :: layer_height ! Height [m]
  real(kind=r_def), intent(in) :: lon, lat     ! Longitude and latitue [rad]
  real(kind=r_def)             :: temp_eq      ! Equilibrium temperature
  real(kind=r_def)             :: theta_eq     ! Equilibrium theta
  real(kind=r_def)             :: dt_eq_factor ! Pressure-dependent DT factor

  ! Local
  real(kind=r_def) :: pressure ! Actual pressure

  ! Calculate pressure-dependent factors
  pressure = p_zero * exner ** (1.0_r_def / kappa)
  ! dt_eq_factor = DT_EQ_POLE * &
  !   max(0.0_r_def, min(1.0_r_def, 1.0_r_def - log10(pressure/ P_LOW) / log10(P_HIGH / P_LOW)))
  dt_eq_factor = DT_EQ_POLE * &
    max(0.0_r_def, min(1.0_r_def, (layer_height - H_2) / (H_1 - h_2)))

  ! Night side
  temp_eq = -0.5_r_def * dt_eq_factor
  if ( abs(lon - LON_SUBSTELLAR) <= 0.5_r_def*pi ) then
    ! Day side
    temp_eq = temp_eq + dt_eq_factor * abs(cos(lon) * cos(lat))
  end if

  theta_eq = theta_init + temp_eq / exner

end function camembert_equilibrium_theta

!> @brief Calculate the inverse radiative time scale the CAMEMBERT Case 1.
!> @param[in] exner Exner pressure
!> @return    freq  Newton cooling relaxation frequency
function camembert_newton_frequency(exner) result(freq)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: exner
  ! Returns
  real(kind=r_def) :: freq
  ! Locals
  real(kind=r_def) :: pressure
  real(kind=r_def) :: tau

  pressure = p_zero * exner ** (1.0_r_def / kappa)

  ! Eq. 2 in Christie et al. (2022)
  tau = (10.0_r_def**2.5_r_def) * (pressure**0.75_r_def)
  tau = max(1.0E4_r_def, min(1.0E7_r_def, tau))

  ! Eq. 3 in Christie et al. (2022)
  ! Using the definition of the universal gas constant:
  ! rd * mu = GAS_CONST
  tau = tau * 4.0_r_def * cp / (7.0_r_def * GAS_CONST)

  freq = 1.0_r_def / tau
  
end function camembert_newton_frequency

end module camembert_forcings_mod
