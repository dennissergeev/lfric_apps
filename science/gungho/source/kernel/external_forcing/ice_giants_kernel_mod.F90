!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Ice Giants temperature forcing
!>
!> @details Kernel that adds Ice Giants temperature forcing kernel following
!!          Guendelman & Kaspi (2025), DOI: 10.1051/0004-6361/202554015; hereafter: GK25
module ice_giants_kernel_mod

  use argument_mod,             only: arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_READWRITE,     &
                                      GH_SCALAR,                 &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_READ, CELL_COLUMN
  use constants_mod,            only: r_def, i_def
  use sci_chi_transform_mod,    only: chi2llr
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  use fs_continuity_mod,        only: Wtheta, Wchi
  use ice_giants_forcings_mod,  only: ice_giants_newton_frequency, &
                                      ice_giants_equilibrium_theta
  use kernel_mod,               only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: ice_giants_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                          &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),                    &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,      Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                                 &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                                  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ice_giants_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ice_giants_code

contains

!> @brief Adds the Ice Giants Temperature Forcing based on GK25
!> @param[in]     nlayers      The number of layers
!> @param[in,out] dtheta       Potential temperature increment data
!> @param[in]     theta        Potential temperature data
!> @param[in]     exner_in_wth The Exner pressure in Wtheta
!> @param[in]     chi_1        First component of the chi coordinate field
!> @param[in]     chi_2        Second component of the chi coordinate field
!> @param[in]     chi_3        Third component of the chi coordinate field
!> @param[in]     panel_id     A field giving the ID for mesh panels
!> @param[in]     kappa        Ratio of Rd and cp
!> @param[in]     dt           The model timestep length
!> @param[in]     ndf_wth      The number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wth     The number of unique degrees of freedom for Wtheta
!> @param[in]     map_wth      Dofmap for the cell at the base of the column for Wtheta
!> @param[in]     ndf_chi      The number of degrees of freedom per cell for Wchi
!> @param[in]     undf_chi     The number of unique degrees of freedom for Wchi
!> @param[in]     map_chi      Dofmap for the cell at the base of the column for Wchi
!> @param[in]     ndf_pid      Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid     Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid      Dofmap for the cell at the base of the column for panel_id
subroutine ice_giants_code( nlayers,                    &
                            dtheta, theta, exner_in_wth,&
                            chi_1, chi_2, chi_3,        &
                            panel_id, kappa, dt,        &
                            ndf_wth, undf_wth, map_wth, &
                            ndf_chi, undf_chi, map_chi, &
                            ndf_pid, undf_pid, map_pid  &
                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid

  real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta
  real(kind=r_def), dimension(undf_wth), intent(in)    :: theta
  real(kind=r_def), dimension(undf_wth), intent(in)    :: exner_in_wth
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def),                      intent(in)    :: kappa
  real(kind=r_def),                      intent(in)    :: dt

  integer(kind=i_def), dimension(ndf_wth),  intent(in) :: map_wth
  integer(kind=i_def), dimension(ndf_chi),  intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in) :: map_pid

  ! Internal variables
  integer(kind=i_def) :: k, df, location, ipanel

  real(kind=r_def)    :: theta_eq, exner
  real(kind=r_def)    :: lat, lon, radius

  real(kind=r_def) :: exner0 ! Lowest-level Exner value
  real(kind=r_def) :: sigma  ! exner/exner0**(1.0/kappa)

  real(kind=r_def) :: coords(3)

  ! Local variables
  real(kind=r_def)    :: newton_frequency ! Relaxation timescale

  coords(:) = 0.0_r_def

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Calculate x, y and z at the centre of the lowest cell
  do df = 1, ndf_chi
    location = map_chi(df)
    coords(1) = coords(1) + chi_1( location )/ndf_chi
    coords(2) = coords(2) + chi_2( location )/ndf_chi
    coords(3) = coords(3) + chi_3( location )/ndf_chi
  end do

  call chi2llr(coords(1), coords(2), coords(3), ipanel, lon, lat, radius)

  exner0 = exner_in_wth(map_wth(1))

  ! Set the relaxation timescale
  newton_frequency = ice_giants_newton_frequency()

  do k = 0, nlayers

    exner = exner_in_wth(map_wth(1) + k)

    sigma = (exner / exner0) ** (1.0_r_def / kappa)

    theta_eq = ice_giants_equilibrium_theta(exner, lat, sigma, kappa)

    dtheta(map_wth(1) + k) = -newton_frequency &
        * (theta(map_wth(1) + k) - theta_eq) * dt

  end do

end subroutine ice_giants_code

end module ice_giants_kernel_mod
