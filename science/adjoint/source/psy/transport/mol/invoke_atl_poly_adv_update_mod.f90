!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Regular invocation of the line-by-line kernel has incorrect
!!        halo depths, so need PSyKAl lite solution. Addressed by
!!        PSyclone issue #2886.
MODULE invoke_atl_poly_adv_update_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: invoke_atl_poly_adv_update_kernel_type

CONTAINS

    SUBROUTINE invoke_atl_poly_adv_update_kernel_type( advective,               &
                                                       ls_reconstruction,       &
                                                       wind,                    &
                                                       wind_dir,                &
                                                       ls_recon_stencil_extent, &
                                                       wind_dir_stencil_extent )

      USE constants_mod, ONLY: r_def, i_def
      USE field_mod, ONLY: field_type, field_proxy_type
      USE atl_poly_adv_update_kernel_mod, ONLY: atl_poly_adv_update_code
      USE mesh_mod, ONLY: mesh_type
      USE stencil_2D_dofmap_mod, ONLY: stencil_2D_dofmap_type, STENCIL_2D_CROSS

      TYPE(field_type), intent(in) :: advective, ls_reconstruction, wind, wind_dir
      INTEGER(KIND=i_def), intent(in) :: ls_recon_stencil_extent, wind_dir_stencil_extent
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: wind_dir_data
      REAL(KIND=r_def), pointer, dimension(:) :: wind_data
      REAL(KIND=r_def), pointer, dimension(:) :: ls_reconstruction_data
      REAL(KIND=r_def), pointer, dimension(:) :: advective_data
      TYPE(field_proxy_type) advective_proxy, ls_reconstruction_proxy, wind_proxy, wind_dir_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_ls_reconstruction(:,:), map_w2(:,:), map_wtheta(:,:)
      INTEGER(KIND=i_def) ndf_wtheta, undf_wtheta, ndf_adspc1_ls_reconstruction, undf_adspc1_ls_reconstruction, ndf_w2, undf_w2
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh
      INTEGER(KIND=i_def) wind_dir_max_branch_length
      INTEGER(KIND=i_def), pointer :: wind_dir_stencil_size(:,:)
      INTEGER(KIND=i_def), pointer :: wind_dir_stencil_dofmap(:,:,:,:)
      TYPE(stencil_2D_dofmap_type), pointer :: wind_dir_stencil_map
      INTEGER(KIND=i_def) ls_reconstruction_max_branch_length
      INTEGER(KIND=i_def), pointer :: ls_reconstruction_stencil_size(:,:)
      INTEGER(KIND=i_def), pointer :: ls_reconstruction_stencil_dofmap(:,:,:,:)
      TYPE(stencil_2D_dofmap_type), pointer :: ls_reconstruction_stencil_map
      !
      ! Initialise field and/or operator proxies
      !
      advective_proxy = advective%get_proxy()
      advective_data => advective_proxy%data
      ls_reconstruction_proxy = ls_reconstruction%get_proxy()
      ls_reconstruction_data => ls_reconstruction_proxy%data
      wind_proxy = wind%get_proxy()
      wind_data => wind_proxy%data
      wind_dir_proxy = wind_dir%get_proxy()
      wind_dir_data => wind_dir_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = advective_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => advective_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Initialise stencil dofmaps
      !
      ls_reconstruction_stencil_map => ls_reconstruction_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,&
&ls_recon_stencil_extent)
      ls_reconstruction_max_branch_length = ls_recon_stencil_extent + 1_i_def
      ls_reconstruction_stencil_dofmap => ls_reconstruction_stencil_map%get_whole_dofmap()
      ls_reconstruction_stencil_size => ls_reconstruction_stencil_map%get_stencil_sizes()
      wind_dir_stencil_map => wind_dir_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,wind_dir_stencil_extent)
      wind_dir_max_branch_length = wind_dir_stencil_extent + 1_i_def
      wind_dir_stencil_dofmap => wind_dir_stencil_map%get_whole_dofmap()
      wind_dir_stencil_size => wind_dir_stencil_map%get_stencil_sizes()
      !
      ! Look-up dofmaps for each function space
      !
      map_wtheta => advective_proxy%vspace%get_whole_dofmap()
      map_adspc1_ls_reconstruction => ls_reconstruction_proxy%vspace%get_whole_dofmap()
      map_w2 => wind_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = advective_proxy%vspace%get_ndf()
      undf_wtheta = advective_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc1_ls_reconstruction
      !
      ndf_adspc1_ls_reconstruction = ls_reconstruction_proxy%vspace%get_ndf()
      undf_adspc1_ls_reconstruction = ls_reconstruction_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for w2
      !
      ndf_w2 = wind_proxy%vspace%get_ndf()
      undf_w2 = wind_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_halo_cell(1)
      !
      ! Call kernels and communication routines
      !
      IF (advective_proxy%is_dirty(depth=1)) THEN
        CALL advective_proxy%halo_exchange(depth=1)
      END IF
      !
      IF (ls_reconstruction_proxy%is_dirty(depth=ls_recon_stencil_extent)) THEN
        CALL ls_reconstruction_proxy%halo_exchange(depth=ls_recon_stencil_extent)
      END IF
      !
      IF (wind_dir_proxy%is_dirty(depth=wind_dir_stencil_extent)) THEN
        CALL wind_dir_proxy%halo_exchange(depth=wind_dir_stencil_extent)
      END IF
      !
      !$omp parallel default(shared), private(cell)
      !$omp do schedule(static)
      DO cell=loop0_start,loop0_stop
        !
        CALL atl_poly_adv_update_code( nlayers,                                      &
                                       advective_data,                               &
                                       ls_reconstruction_data,                       &
                                       ls_reconstruction_stencil_size(:,cell),       &
                                       ls_reconstruction_max_branch_length,          &
                                       ls_reconstruction_stencil_dofmap(:,:,:,cell), &
                                       wind_data,                                    &
                                       wind_dir_data,                                &
                                       wind_dir_stencil_size(:,cell),                &
                                       wind_dir_max_branch_length,                   &
                                       wind_dir_stencil_dofmap(:,:,:,cell),          &
                                       ndf_wtheta,                                   &
                                       undf_wtheta,                                  &
                                       map_wtheta(:,cell),                           &
                                       ndf_adspc1_ls_reconstruction,                 &
                                       undf_adspc1_ls_reconstruction,                &
                                       map_adspc1_ls_reconstruction(:,cell),         &
                                       ndf_w2,                                       &
                                       undf_w2,                                      &
                                       map_w2(:,cell) )
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL advective_proxy%set_dirty()
      CALL advective_proxy%set_clean(1)
      CALL wind_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_atl_poly_adv_update_kernel_type

END MODULE invoke_atl_poly_adv_update_mod
