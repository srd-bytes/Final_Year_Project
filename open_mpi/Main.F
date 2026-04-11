! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! 3D MPI-3 and Fortran-based Adaptive Mesh Refinement version of Cubed Sphere
! based MHD code from Dinshaw Balsara's group. Uses Silo for I/O. 
! An MPI-twinned Fortran compiler, along with Silo libraries, is all one needs to run.

! Sectors and Patches:-
! Code operates on six spherical sectors, each of which may be split into
! many (logically) cubical patches. Each patch has "ntp_zonesXntp_zones"
! in the theta-phi directions and "nr_zones" in the radial direction.
! Each patch has "NGhost" number of ghost zones.

! Each Sector:-
! Within each sector, the base-level mesh can have "ntp_patches X ntp_patches" in
! the theta-phi directions and "nr_patches" in the radial direction. "ilevel_max"
! number of levels can be present; so that each base level zone can be split
! into up to "iscale**3" number of sub-cells.

! Global Mesh:-
! The above information can be used to build "ntp_global", the number of zones
! in any one theta-phi direction on any of the 6 sectors. Similarly, we have
! "nr_global" and "NGhost_global". Each sector has a complete 2D lattice
! on the unit sphere at the finest level, so that the spherical geometry can
! be extracted at any desired level from that.

! On the base mesh, we also retain "ntp_global_base" and "nr_global_base" so
! that building sibling relationships is very easy for one-level mesh.

! Active Patches:-
! "Max_base_patches_per_sector" is the max patches per sector on the base level.

! "Max_Theoretical_patches_per_sector" is the max number of theoretical patches
! that we can have in one sector on a "ilevel_max" global hierarchy.
! We keep track of whether that patch location is filled or not. If it
! is filled, we provide its patch number in that sector.

! Reality is that we only have "Max_patches_per_proc" on each processor.

! Processor numbering starts from 1; so remember to subtract 1 for MPI's C format.

! MHD Variables:-
! Thre are "n_cc_components" number of zone-centered variables which are
! reconstructed (based on order) to have up to "n_cc_modes". Each face has number of
! vector field "n_fc_components" which can each be reconstructed to have
! "n_fc_modes". After div-preserving reconstruction, we will have "n_cc_bfld_modes"
! number of modes for each component of the vector field.

! ----------------------------------------------------------------------

#include "directives"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      MODULE RIEMANN_3D_CS_AMR_SUPPORT

! Written by : Dinshaw Balsara

      IMPLICIT NONE

      include 'mpif.h'

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! List of public routines for this module:

      PUBLIC :: BUILD_BASE_PATCHES,

! Subroutines for creating MPI Windows and freeing them up.
     1      WINCREATE_MESSAGING, WINFREE_MESSAGING,

! Load Local Buffers; first load them in x and y-directions. Then load in z-direction.
     1      LOAD_LOCAL_BUFFERS_UXY, LOAD_LOCAL_BUFFERS_UZ,

     1      LOAD_LOCAL_BUFFERS_Flux_XY, LOAD_LOCAL_BUFFERS_Flux_Z,

! Fetch buffers from neighbors and initialize suitably.
     1      FETCH_BUFFERS_FILL_IN_UX, FETCH_BUFFERS_FILL_IN_UY,
     1      FETCH_BUFFERS_FILL_IN_UZ, FETCH_BUFFERS_FILL_IN_UY_TOP,

     1      FETCH_BUFFERS_FILL_IN_Flux_X, FETCH_BUFFERS_FILL_IN_Flux_Y,
     1      FETCH_BUFFERS_FILL_IN_Flux_Z,
!    1      FETCH_BUFFERS_FILL_IN_Flux_Y_TOP,
     1      SYNCHRONIZE_EFZ_QUAD_CORNER,

     1      LOAD_AMR_FLUXES_AT_BDY,

! Testing Subroutine, Print amr variable values in a sector.
     1      PRINT_VARIABLES_ONE_PATCH_PER_SECTOR,
     1      TEST_XY_MESSAGING, TEST_MESSAGING_ERROR_CALCULATION,

#if( ACCURACY_ANALYSIS == 1 )
     1      MAKE_GLOBAL_ACCURACY_CS,
#endif

! Silo subroutines.
     1      WRITE_SILO_ONE_BLOCK

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! START LIST of parameters that control Cube Sphere AMR mesh and its processing.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      INTEGER ntp_zones, nr_zones, NGhost, ntp_patches, nr_patches,

     1        ilevel_max, iscale,

     1        ntp_global, nr_global, NGhost_global,

     1        ntp_global_base, nr_global_base,

     1        Max_base_patches_per_sector, 
     1        Max_Theoretical_patches_per_sector, Max_patches_per_proc,

     1        n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1        n_fc_components, n_fc_modes, n_cc_bfld_modes

      PARAMETER

! Variables that specify the number of zones "ntp_zones, nr_zones" and the
! number of ghost zones "NGhost" for a single patch.
! For 480 cores at medium resolution:
! ntp_zones = 18, ntp_patches = 4, nr_zones = 60, nr_patches = 5,
! For 960 cores at medium resolution:
! ntp_zones = 18, ntp_patches = 4, nr_zones = 30, nr_patches = 10,
! For 960 cores at higher resolution:
! ntp_zones = 20; ntp_patches = 4; nr_zones = 36; nr_patches = 10;
! 
     1          ( ntp_zones = 18, nr_zones =  60, NGhost = 4,

! Number of patches "ntp_patchesXntp_patchesXnr_patches" within a sector on
! the global base-level mesh.
     1            ntp_patches = 4, nr_patches = 10,

! Maximum level in AMR hierarchy and scale between base and finest meshes.
     1            ilevel_max = 1, iscale = (2**ilevel_max) / 2,

! Number of zones within a sector, assuming fully filled levels, at the finest AMR level.
     1            ntp_global = ntp_zones * ntp_patches * iscale,
     1            nr_global = nr_zones * nr_patches * iscale,
     1            NGhost_global = NGhost * iscale,

! Number of zones within a sector on the base mesh.
     1            ntp_global_base = ntp_zones * ntp_patches,
     1            nr_global_base = nr_zones * nr_patches,

! Max theoretical patches within a sector for "ilevel_max" hierarchy.
! Also, number of patches at base level within a sector and max patches
! actually allowed on a processor.

     1      Max_base_patches_per_sector = ntp_patches**2 * nr_patches,

     1      Max_Theoretical_patches_per_sector =  ! Sum geometric series with ratio 8.
     1      Max_base_patches_per_sector * ( 8**ilevel_max - 1) / 7,

     1      Max_patches_per_proc = Max_base_patches_per_sector,  ! See "BUILD_BASE_PATCHES".

! MHD parameters.
     1            n_cc_components = NFIELD + 1 + NFLUID,
     1            n_fc_components = 1, n_fc_modes = 1,

#if ( NUMORDER == 2)
     1            n_cc_modes = 4, n_cc_modes_one_higher = 10,
     1            n_cc_bfld_modes = 7
#elif ( NUMORDER == 3)
     1            n_cc_modes = 10, n_cc_modes_one_higher = 20,
     1            n_cc_bfld_modes = 16
#elif ( NUMORDER == 4)
     1            n_cc_modes = 20, n_cc_modes_one_higher = 35,
     1            n_cc_bfld_modes = 30
#endif

     1            )

#if( USE_B0B1_FORMULATION == 1 )
      CHARACTER*SILO_VARLEN silo_allvars ( NFIELD + 1 + NFLUID + 3 
     1                                   + SILO_EXTRA_DATA)
#else
      CHARACTER*SILO_VARLEN silo_allvars ( NFIELD + 1 + NFLUID
     1                                   + SILO_EXTRA_DATA)
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Parent-Child Sibling information is stored on each processor so that
! the full composite mesh is easily accessed.

      INTEGER

! (-1,-1,-1) if the patch is inactive. If active data is present,
! "iactive_FromSector" gives (sector#, patch#, processsor#).

     1   iactive_FromSector ( 6, Max_Theoretical_patches_per_sector, 3),

! Uses "iactive_FromSector" to build the sector of the neighboring rt,lf,up,dn,tp,bt
! siblings in "isibling_sector_t". Also gives the number of the neighbor
! in "isibling_number_t" and the processor on which this sibling's data
! resides "isibling_proc_t".

     1   isibling_sector_t ( 6, Max_Theoretical_patches_per_sector, 7),
     1   isibling_number_t ( 6, Max_Theoretical_patches_per_sector, 7),
     1   isibling_proc_t ( 6, Max_Theoretical_patches_per_sector, 7),

! The desired usage is that we have "Max_patches_per_proc" number of patches
! on each processor. For each such active patch, we want its sector number
! and its logical patch number in that sector.

     1   iactive_FromProc ( Max_patches_per_proc, 2)

! ----------------------------------------------------------------------

! Minimal amounts of coordinate information is stored from which we can
! regenerate the geometry.

! "coordinates_gr_global" contains the coordinates of the vertices of the
! finest mesh that covers the 2D surface of the unit sphere. As a result, it has
! information for all six sectors. This is built once at the start of execution.

! From this information, each patch at each level in the AMR hierarchy can
! extract its own "vert_coords_patch, imask_active_patch" which allows the
! geometry for that patch to be quickly rebuilt.

      REAL coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

! "rb_save_global, dr_save_global, rc_save_global" gives us radial coordinate
! information on the finest level. Actually, only "rb_save_global" is used to
! build "vert_coords_patch".

      REAL rb_save_global ( 0 - NGhost_global:
     1                      nr_global + NGhost_global)

      REAL, DIMENSION ( 1  - NGhost_global: nr_global + NGhost_global)::
     1     dr_save_global, rc_save_global

! ----------------------------------------------------------------------

! Minimal amounts of geometry information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor. This is so that all other
! geometry information can be quickly regenerated from it.

! Vertices of each of the patches.
! "vert_coords_amr" for all patches is matched by "vert_coords_patch" for single patches.

      REAL vert_coords_amr
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3, Max_patches_per_proc)

! Mask to decide whether that zone can be touched; some zones are tri-corners.
! "imask_active_amr" for all patches is matched by "imask_active_patch" for single patches.

      INTEGER imask_active_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost, Max_patches_per_proc)

! The stencil info. is only two-dimensional and worth storing.
! "istencil_amr" for all patches is matched by "istencil_patch" for single patches.

      INTEGER istencil_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost: ntp_zones + NGhost, 2, 12,
     1                       Max_patches_per_proc)

! For each patch, keep track whether its four corners are tri-corners or quad-corners.
! "itri_patch_amr" for all patches empties into
! "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU" for the patches.

      INTEGER itri_patch_amr ( 4, Max_patches_per_proc)

! Sector and level information for each patch.
! "isector_ilevel_amr" for all patches is matched by "isector" which is the sector
! number for this patch.

      INTEGER isector_ilevel_amr ( 2, Max_patches_per_proc)

! Which "ix_patch, iy_patch, iz_patch" is this patch?
! "ixyz_patch_amr" for all patches empties into "ix_patch, iy_patch, iz_patch"
! for each patch.

      INTEGER ixyz_patch_amr ( 3, Max_patches_per_proc)

! Is top or bottom patch. i.e., has border at rmin or rmax. Just needed for ADD
! and PAD boundary conditions.

      INTEGER is_bottom_or_top ( 2, Max_patches_per_proc)

! --------------------------------------------------
! The above variables represent the minimal amount of geometry information that
! we need to store in order to rebuild all the geometry variables.
! But we also strategically store a little more geometry
! information that is used just for the zone update.

! "xf_quadnorml_amr" for all patches is matched by "xf_quadnorml_patch" for single patches.

      REAL xf_quadnorml_amr ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! "xf_area_amr" for all patches is matched by "xf_area_patch" for single patches.

      REAL xf_area_amr ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! "yf_quadnorml_amr" for all patches is matched by "yf_quadnorml_patch" for single patches.

      REAL yf_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 
     1   3, Max_patches_per_proc)

! "yf_area_amr" for all patches is matched by "yf_area_patch" for single patches.

      REAL yf_area_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! "zf_quadnorml_amr" for all patches is matched by "zf_quadnorml_patch" for single patches.

      REAL zf_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER, Max_patches_per_proc)

! "zf_area_amr" for all patches is matched by "zf_area_patch" for single patches.

      REAL zf_area_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! ------------------------------
! "v_scale_length_amr" for all patches is matched by "v_scale_length" for single patches.

      REAL v_scale_length_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! ------------------------------
! "v_volume_ோamr" for all patches is matched by "v_volume_patch" for single patches.

      REAL v_volume_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! "v_centroid_amr" for all patches is matched by "v_centroid_patch" for single patches.

      REAL v_centroid_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ----------------------------------------------------------------------

! Minimal amounts of MHD information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor.

! Each and every "var_name_amr" variable for the fluid variables will have
! a counterpart "var_name" for the patches.

! Zone-centered MHD variables for all patches on this processor.

      REAL U_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, Max_patches_per_proc)

#if( USE_B0B1_FORMULATION == 1 )

! the silo read/write variable is sized as "n_cc_components + B0 ( 3) + entropy".
      REAL U_silo_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost, 
     1                  1 - NGhost:  nr_zones + NGhost, 
     1   n_cc_components + 3 + SILO_EXTRA_DATA, Max_patches_per_proc)

#else

! the silo read/write variable is sized as "n_cc_components + entropy".
      REAL U_silo_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost, 
     1                  1 - NGhost:  nr_zones + NGhost, 
     1   n_cc_components + SILO_EXTRA_DATA, Max_patches_per_proc)

#endif

      REAL B0field_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! Zone-centered entropy density variable for all patches on this processor.

      REAL entropygr_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! Extra variables for silo output. divergence and entropy, tempfloor flags.

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost, 
     1                  1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1   divgr_amr, flag_entprs_amr, flag_Tfloor_amr

! Face-centered magnetic field variables for all patches on this processor.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1   Bxf_amr, Byf_amr, Bzf_amr

! Zone-centered flattener for for all patches on this processor.

      REAL flatten_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! Zone-centered Source for for all patches on this processor.

      REAL Source_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, Max_patches_per_proc)

! ----------

! Next predicted timestep for each patch on this processor.

      REAL dtcurp1_amr ( Max_patches_per_proc)

! Reduced minimum of next predicted timestep from all patches on this processor.

      REAL dtcurp1_reduce

! Min and max of variables for each patch on this processor.

      REAL, DIMENSION ( n_cc_components + 1, Max_patches_per_proc) ::
     1     V_min_amr, V_max_amr

! Reduced min and max variables from all patches on this processor.

      REAL, DIMENSION ( n_cc_components + 1) ::
     1    V_min_reduce, V_max_reduce

! ----------
#if ( ACCURACY_ANALYSIS == 1)

      REAL, DIMENSION ( n_cc_components + 1, Max_patches_per_proc) ::    
     1     U_L1_amr, U_Linf_amr, U_zc_Linf_amr

      REAL, DIMENSION ( n_cc_components + 1) ::    
     1     U_L1,        U_Linf,        U_zc_Linf,
     1     U_L1_reduce, U_Linf_reduce, U_zc_Linf_reduce 

#endif

! ----------------------------------------------------------------------

! Minimal amounts of flux and electric field information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor.
! The defects in these fluxes and electric fields are also stored here.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components, Max_patches_per_proc) ::
     1      Flux_x_amr, Flux_y_amr, Flux_z_amr,
     1      Defect_Flux_x_amr, Defect_Flux_y_amr, Defect_Flux_z_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1      EntFlux_x_amr, EntFlux_y_amr, EntFlux_z_amr,
     1      Defect_EntFlux_x_amr, Defect_EntFlux_y_amr,
     1      Defect_EntFlux_z_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1      Efx_amr, Efy_amr, Efz_amr,
     1      Defect_Efx_amr, Defect_Efy_amr, Defect_Efz_amr

! ----------------------------------------------------------------------
! Buffer arrays that store "0:N_Ghost" layers to the Right and Left of current patch.
! Here we set "n_cc_components + 4" so that we can move entropy and field information
! also in one single message.

      REAL, DIMENSION ( 0: NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_r_amr, U_snd_l_amr

      REAL, DIMENSION ( 0: NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4) :: U_recev_r, U_recev_l

! ----------
! Layers to the Up and Down of current patch.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0: NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_u_amr, U_snd_d_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0: NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4) :: U_recev_u, U_recev_d

! ----------
! Layers to the Top and Bottom of current patch.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost, 0: NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_t_amr, U_snd_b_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost, 0: NGhost,
     1   n_cc_components + 4) :: U_recev_t, U_recev_b

! ----------------------------------------------------------------------
! Buffer arrays that store one layer of fluxes and electric fields to the
! Right and Left of current patch. Here we set "n_cc_components + 3" so that
! we can move entropy flux and "Ey and Ez" in one single message.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 3, Max_patches_per_proc) ::
     1   Flux_snd_r_amr, Flux_snd_l_amr,
     1   Flux_defect_r_amr, Flux_defect_l_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 3) ::
     1   Flux_snd_r_patch, Flux_snd_l_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost, 
     1   0 - NGhost: nr_zones + NGhost, n_cc_components + 3) ::
     1   Flux_recev_r, Flux_recev_l

! ----------
! Layers to the Up and Down of current patch. "Ez and Ex" saved in that order.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 3, Max_patches_per_proc) ::
     1   Flux_snd_u_amr, Flux_snd_d_amr,
     1   Flux_defect_u_amr, Flux_defect_d_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 3) ::
     1   Flux_snd_u_patch, Flux_snd_d_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, n_cc_components + 3) ::
     1   Flux_recev_u, Flux_recev_d

! ----------
! Layers to the Top and Bottom of current patch. "Ex and Ey" saved in that order.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost,
     1  n_cc_components + 3, Max_patches_per_proc) ::
     1   Flux_snd_t_amr, Flux_snd_b_amr,
     1   Flux_defect_t_amr, Flux_defect_b_amr

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost,
     1  n_cc_components + 3) ::
     1   Flux_snd_t_patch, Flux_snd_b_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost, n_cc_components + 3) ::
     1   Flux_recev_t, Flux_recev_b

! ----------------------------------------------------------------------

! A layer "NGhost" deep of conserved variables (U) and fields will have to
! be messaged in the x-right (XR) direction. Similar for the x-left (XL) direction.
! Similar for the y-up (YU) direction and y-down (YD) direction.
! Similar for the z-top (ZT) direction and z-bottom (ZB) direction.

! We also build the analogous MPI handles for messaging of fluxes.

      INTEGER(KIND=MPI_ADDRESS_KIND) :: 
     1        window_size_xy, window_size_z, window_offset

      INTEGER Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
     1        Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
     1        Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle,

     1        Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
     1        Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
     1        Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle,

     1        type_mpi_real  ! Holds the number of bytes for real datatype.

! ----------------------------------------------------------------------





! ----------------------------------------------------------------------

! ON-PROCESSOR, SINGLE-PATCH GEOMETRY VARIABLES:-

! Local data, regenerated from "coordinates_gr_global" and "rb_save_global", that
! contains the full geometry information that is being processed on this processor.

! First, this data is used to build "vert_coords_patch" and "imask_active_patch".

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: nr_zones + NGhost)

! --------------------------------------------------

! Next, from the above two, all other geometry and quadrature information is built.

      REAL xf_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER)

      REAL xf_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER)

      REAL xf_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL yf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER)

      REAL yf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER)

      REAL yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL zf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER)

      REAL zf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER)

      REAL zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER)

      REAL zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost),
     1     zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

! ------------------------------

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------

      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

      INTEGER, DIMENSION ( n_cc_modes_one_higher) ::
     1         ii_mode, jj_mode, kk_mode

      INTEGER, DIMENSION ( n_cc_bfld_modes) ::
     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes_one_higher)

! ------------------------------

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------

! ON-PROCESSOR, SINGLE-PATCH FLUID VARIABLES, INCLUDING _SKINNY VERSIONS:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL U_skinny ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

      REAL entropygr_skinny ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   Bxf_skinny, Byf_skinny, Bzf_skinny

      REAL B0field ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_bfld_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL Source ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,  
     1   n_cc_components)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr, workgr

      REAL, DIMENSION ( n_cc_components + 1) :: V_min_local, V_max_local

! Extra variables for silo output; entropy, tempfloor flags.

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost, 
     1                  1 - NGhost:  nr_zones + NGhost) ::
     1   flag_entprs, flag_Tfloor

! ------------------------------
! Variables used for div-free reconstruction of magnetic fields:-

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           NUMORDER * NUMORDER) ::
     1   Bxf_quadpts, Byf_quadpts, Bzf_quadpts

! ------------------------------
! VARIABLES NEEDED FOR 1 STAGE PREDICTOR TIMESTEPPING:-
! (Currently, this is only specialized for second order.)

      REAL dU_dt_predgr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL dEnt_dt_predgr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL U_TAvg_modes_gr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_bfld_modes)

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Fluid fluxes and resolved states from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components) ::
     1      Flux_x, Flux_y, Flux_z

! Entropy flux from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

! Electric fields at the edges of the mesh.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Efx, Efy, Efz

! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE BUILD_BASE_PATCHES ()

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! for base level patches.
! Every processor will build this data locally so that it knows how to access
! data that might be on other sectors.

! When setting the number of processors and "Max_patches_per_proc", pay
! special attention to the logical flow of this subroutine.
! The check should be :-
! "6 * Max_base_patches_per_sector == nimages * npatch_per_proc".

! For load balance, "nimages" must be a multiple of 6; the number of sectors.

! So, once you know the total number of MPI processors "nimages",
! make sure that "Max_patches_per_proc .GE. npatch_per_proc".

! Once the above information is built, we set up "iactive_FromProc".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me,
     1        ixx, iyy, izz, isector, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1,
     1        isector_rt, isector_lf, isector_up, isector_dn,
     1        ixxn, iyyn, iin, ipatch

      INTEGER neighbor_sector_list ( 6, 4)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

! ----------

! We will have a total of "6 * Max_base_patches_per_sector" patches on the base mesh.
! If that does not divide evenly into "nimages", we have a load imbalance
! on the base level mesh. This should be treated as unacceptable.

! Please note that this check also forces "nimages" to be a multiple of six!

      IF ( ( MOD ( 6 * Max_base_patches_per_sector, nimages) .NE. 0)
     1   .AND. ( me .EQ. 1) ) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES : patches unbalanced"
      WRITE ( 6, *)"6 * Max_base_patches_per_sector, nimages = ",
     1             6 * Max_base_patches_per_sector, nimages
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------

! Find number of patches we can put per processor. It is bounded by
! "Max_patches_per_proc" so check for that.

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

      IF ( ( npatch_per_proc .GT. Max_patches_per_proc) .AND.
     1     ( me .EQ. 1)) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES : not enough patches"
      WRITE ( 6, *)"npatch_per_proc, Max_patches_per_proc = ",
     1              npatch_per_proc, Max_patches_per_proc
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------------------------------------------------------------------
! Set up "neighbor_sector_list ( 6, 4)". First index is parent sector.
! Second index is 1 (right), 2 (left), 3 (up), 4 (down) sectors relative
! to parent sector. z-direction does not need this mapping.

       neighbor_sector_list ( 1, 1) = 3
       neighbor_sector_list ( 1, 2) = 6
       neighbor_sector_list ( 1, 3) = 2
       neighbor_sector_list ( 1, 4) = 5

       neighbor_sector_list ( 2, 1) = 1
       neighbor_sector_list ( 2, 2) = 4
       neighbor_sector_list ( 2, 3) = 3
       neighbor_sector_list ( 2, 4) = 6

       neighbor_sector_list ( 3, 1) = 2
       neighbor_sector_list ( 3, 2) = 5
       neighbor_sector_list ( 3, 3) = 1
       neighbor_sector_list ( 3, 4) = 4

       neighbor_sector_list ( 4, 1) = 5
       neighbor_sector_list ( 4, 2) = 2
       neighbor_sector_list ( 4, 3) = 6
       neighbor_sector_list ( 4, 4) = 3

       neighbor_sector_list ( 5, 1) = 6
       neighbor_sector_list ( 5, 2) = 3
       neighbor_sector_list ( 5, 3) = 4
       neighbor_sector_list ( 5, 4) = 1

       neighbor_sector_list ( 6, 1) = 4
       neighbor_sector_list ( 6, 2) = 1
       neighbor_sector_list ( 6, 3) = 5
       neighbor_sector_list ( 6, 4) = 2

! ----------------------------------------------------------------------

! The first "Max_base_patches_per_sector" on each sector will be active and their
! sector number (first index) will point to themselves.

      DO isector = 1, 6

      iactive_FromSector ( isector, 1: Max_base_patches_per_sector, 1) =
     1                     isector

      END DO

! ----------

! The first "Max_base_patches_per_sector" on each sector will have their patch number
! (second index) be their ordinal number within the sector.

! The processor on which the patch lands is given by the third index.
! Each sector has "Max_Base_patches_per_sector" number of patches. But the patches
! in a single sector can be distributed across multiple processors.
! The data will be contiguous in order to facilitate efficient messaging.

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.

          ii = ii + 1

          iactive_FromSector ( :, ii, 2) = ii

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.
        END DO  ! End of "DO iyy = 1, ntp_patches" loop.
      END DO  ! End of "DO izz = 1, nr_patches" loop.

! ----------

      IF ( ( ii .NE. Max_base_patches_per_sector) .AND.
     1     ( me .EQ. 1)) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES :"
      WRITE ( 6, *)"ii .NE. Max_base_patches_per_sector = ",
     1              ii, Max_base_patches_per_sector
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------

! The processor on which the patch lands is given by the third index.
! Each sector has "Max_Base_patches_per_sector" number of patches. But the patches
! in a single sector can be distributed across multiple processors.
! The data will be contiguous in order to facilitate efficient messaging.

      DO isector = 1, 6

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.

          ii = ii + 1

          iactive_FromSector ( isector, ii, 3) = 
     1    1 + (ii + (isector - 1) * jj - 1) / npatch_per_proc

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.
        END DO  ! End of "DO iyy = 1, ntp_patches" loop.
      END DO  ! End of "DO izz = 1, nr_patches" loop.

      END DO  ! End of "DO isector = 1, 6" loop.

! ----------------------------------------------------------------------
! Now start building "isibling_number_t, isibling_sector_t, isibling_proc_t".
! There are 7 neighbors because the first neighbor is to the self.

! First we fill in the self identification.

      isibling_sector_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 1)

      isibling_number_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 2)

      isibling_proc_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 3)

! ----------
! Loop over sectors and find their neighbors.

      DO isector = 1, 6

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
      izzp1 = izz + 1
      izzm1 = izz - 1

        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
        iyyp1 = iyy + 1
        iyym1 = iyy - 1

          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.
          ixxp1 = ixx + 1
          ixxm1 = ixx - 1

          ii = ii + 1

! ----------------------------------------
! Find right neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( ixxp1 .LE. ntp_patches) THEN
! Right neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 2) = isector

            isibling_number_t ( isector, ii, 2) = ii + 1

            isibling_proc_t ( isector, ii, 2) = 
     1      1 + (ii + (isector - 1) * jj) / npatch_per_proc

! ----------

          ELSE
! Right neighbor is out-of-sector.

! ----------

! This is the right neighboring sector.
            isector_rt = neighbor_sector_list ( isector, 1)

! On right neighboring sector, x and y get flipped; see figure in notes.
! Find "iin" as the right sibling's patch number on sector "isector_rt".
            ixxn = iyy
            iyyn = ntp_patches
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 2) = isector_rt

            isibling_number_t ( isector, ii, 2) =
     1      iactive_FromSector ( isector_rt, iin, 2)

            isibling_proc_t ( isector, ii, 2) =
     1      iactive_FromSector ( isector_rt, iin, 3)

! ----------

          END IF  ! End of "IF ( ixxp1 .LE. ntp_patches) THEN".

! ----------------------------------------
! Find left neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( ixxm1 .GE. 1) THEN
! Left neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 3) = isector

            isibling_number_t ( isector, ii, 3) = ii - 1

            isibling_proc_t ( isector, ii, 3) =
     1      1 + (ii + (isector - 1) * jj - 2) / npatch_per_proc

! ----------

          ELSE
! Left neighbor is out-of-sector.

! ----------

! This is the left neighboring sector.
            isector_lf = neighbor_sector_list ( isector, 2)

! On left neighboring sector, x and y are not flipped but y-sequence
! is reversed; see figure in notes.
! Find "iin" as the left sibling's patch number on sector "isector_lf".
            ixxn = 1
            iyyn = ntp_patches - iyy + 1
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 3) = isector_lf

            isibling_number_t ( isector, ii, 3) = 
     1      iactive_FromSector ( isector_lf, iin, 2)

            isibling_proc_t ( isector, ii, 3) =
     1      iactive_FromSector ( isector_lf, iin, 3)

! ----------

          END IF  ! End of "IF ( ixxm1 .GE. 1) THEN".

! ----------------------------------------
! Find upper neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( iyyp1 .LE. ntp_patches) THEN
! Upper neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 4) = isector

            isibling_number_t ( isector, ii, 4) = ii + ntp_patches

            isibling_proc_t ( isector, ii, 4) =
     1      1 + (ii + ntp_patches + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Upper neighbor is out-of-sector.

! ----------

! This is the upper neighboring sector.
            isector_up = neighbor_sector_list ( isector, 3)

! On the upper neighboring sector, x and y get flipped; see figure in notes.
! Find "iin" as the upper sibling's patch number on sector "isector_up".
            ixxn = ntp_patches
            iyyn = ixx
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 4) = isector_up

            isibling_number_t ( isector, ii, 4) = 
     1      iactive_FromSector ( isector_up, iin, 2)

            isibling_proc_t ( isector, ii, 4) =
     1      iactive_FromSector ( isector_up, iin, 3)

! ----------

          END IF  ! End of "IF ( iyyp1 .LE. ntp_patches) THEN".

! ----------------------------------------
! Find down neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( iyym1 .GE. 1) THEN
! Down neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 5) = isector

            isibling_number_t ( isector, ii, 5) = ii - ntp_patches

            isibling_proc_t ( isector, ii, 5) =
     1      1 + (ii - ntp_patches + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Down neighbor is out-of-sector.

! ----------

! This is the down neighboring sector.
            isector_dn = neighbor_sector_list ( isector, 4)

! On the down neighboring sector, x and y are not flipped but x-sequence
! is reversed; see figure in notes.
! Find "iin" as the down sibling's patch number on sector "isector_dn".
            ixxn = ntp_patches - ixx + 1
            iyyn = 1
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 5) = isector_dn

            isibling_number_t ( isector, ii, 5) = 
     1      iactive_FromSector ( isector_dn, iin, 2)

            isibling_proc_t ( isector, ii, 5) =
     1      iactive_FromSector ( isector_dn, iin, 3)

! ----------

          END IF  ! End of "IF ( iyym1 .GE. 1) THEN".

! ----------------------------------------
! Find top neighbor; could be in-sector (easy) or non-existent.

          IF ( izzp1 .LE. nr_patches) THEN
! Top neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 6) = isector

            isibling_number_t ( isector, ii, 6) = ii + ntp_patches**2

            isibling_proc_t ( isector, ii, 6) = 
     1      1 + (ii + ntp_patches**2 + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Top neighbor is non-existent (will be filled in with boundary condition routine).

! ----------

            isibling_sector_t ( isector, ii, 6) = -1

            isibling_number_t ( isector, ii, 6) = -1

            isibling_proc_t ( isector, ii, 6) = -1

! ----------

          END IF  ! End of "IF ( izzp1 .LE. nr_patches) THEN".

! ----------------------------------------
! Find bottom neighbor; could be in-sector (easy) or non-existent.

          IF ( izzm1 .GE. 1) THEN
! Bottom neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 7) = isector

            isibling_number_t ( isector, ii, 7) = ii - ntp_patches**2

            isibling_proc_t ( isector, ii, 7) =
     1      1 + (ii - ntp_patches**2 + (isector - 1) * jj - 1) 
     1          / npatch_per_proc

! ----------

          ELSE
! Bottom neighbor is non-existent (will be filled in with boundary condition routine).

! ----------

            isibling_sector_t ( isector, ii, 7) = -1

            isibling_number_t ( isector, ii, 7) = -1

            isibling_proc_t ( isector, ii, 7) = -1

! ----------

          END IF  ! End of "IF ( izzm1 .GE. 1) THEN".

! ----------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.

        END DO  ! End of "DO iyy = 1, ntp_patches" loop.

      END DO  ! End of "DO izz = 1, nr_patches" loop.

      END DO  ! End of "DO isector = 1, 6" loop.

! ----------------------------------------------------------------------
! Now we set up "iactive_FromProc". Here we are looping over the list of
! patches that reside on each processor. These could be on any sector and
! have any logical patch number. At least for the base mesh, this array
! gives us the (sector#, logical patch#).

      iactive_FromProc ( :, :) = -1
   
      DO ii = 1, npatch_per_proc

      jj = ii + (me - 1) * npatch_per_proc

      isector = 1 + ( ( jj - 1) / Max_base_patches_per_sector)
      ipatch = jj - (isector - 1) * Max_base_patches_per_sector

      iactive_FromProc ( ii, 1) = isector  ! Which sector is this patch in?
      iactive_FromProc ( ii, 2) = ipatch  ! logical patch number in "isector".

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE BUILD_BASE_PATCHES

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LOAD_LOCAL_BUFFERS_UXY ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after each timestep iteration, so that the
! final, updated data is in "U_amr, entropygr_amr, Bxf_amr, Byf_amr,
! Bzf_amr". The date is stored in the big container arrays with this
! ordering of variables. This data is loaded into
! "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr, U_snd_d_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      U_snd_r_amr ( :, :, :, :, :) = 0.0
      U_snd_l_amr ( :, :, :, :, :) = 0.0
      U_snd_u_amr ( :, :, :, :, :) = 0.0
      U_snd_d_amr ( :, :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "U_snd_r_amr, U_snd_l_amr".
! Then store up and down y-data in "U_snd_u_amr, U_snd_d_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "U_snd_r_amr, U_snd_l_amr". First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, 1: n_cc_components, ii) =
     1    U_amr       ( ixr, iyy, izz, 1: n_cc_components, ii)

          U_snd_l_amr ( ixs, iyy, izz, 1: n_cc_components, ii) =
     1    U_amr       ( ixl, iyy, izz, 1: n_cc_components, ii)


          U_snd_r_amr   ( ixs, iyy, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixr, iyy, izz, ii)

          U_snd_l_amr   ( ixs, iyy, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variable "Bxf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 2, ii) =
     1    Bxf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 2, ii) =
     1    Bxf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 3, ii) =
     1    Byf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 3, ii) =
     1    Byf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 4, ii) =
     1    Bzf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 4, ii) =
     1    Bzf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "U_snd_u_amr, U_snd_d_amr". First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyu, izz, 1: n_cc_components, ii)

          U_snd_d_amr ( ixx, iys, izz, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyd, izz, 1: n_cc_components, ii)


          U_snd_u_amr       ( ixx, iys, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr       ( ixx, iys, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bxf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice

! --------------------

          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc"

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_UXY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LOAD_LOCAL_BUFFERS_UZ ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after each timestep iteration, so that the
! final, updated data is in "U_amr, entropygr_amr, Bxf_amr, Byf_amr,
! Bzf_amr". The date is stored in the big container arrays with this
! ordering of variables. This data is loaded into "U_snd_t_amr, U_snd_b_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        izt, izb, izs

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      U_snd_t_amr ( :, :, :, :, :) = 0.0
      U_snd_b_amr ( :, :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! Store top and bottom z-data in "U_snd_t_amr, U_snd_b_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "U_snd_t_amr, U_snd_b_amr".  First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyy, izt, 1: n_cc_components, ii)

          U_snd_b_amr ( ixx, iyy, izs, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyy, izb, 1: n_cc_components, ii)


          U_snd_t_amr       ( ixx, iyy, izs, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr       ( ixx, iyy, izs, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bxf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, 0".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_UZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LOAD_AMR_FLUXES_AT_BDY ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called before the synchronization of fluxes and 
! electric fields at the coinciding faces and edges, respectively.
! For that, the amr fluxes and electric fields at the coinciding faces and edges
! are populated in this subroutine from the "Flux_snd_X_amr" arrays.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        izt, izb,

     1        ixs, iys, izs

! ----------------------------------------------------------------------
! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "Flux_snd_r_amr, Flux_snd_l_amr".
! Then store up and down y-data in "Flux_snd_u_amr, Flux_snd_d_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_amr, Flux_snd_l_amr". First for the face centered fluxes,
! "Flux_x_amr, EntFlux_x_amr" then for the edge centered variables,
! "Efy_amr, Efz_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_x_amr, EntFlux_x_amr" data in
! "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Flux_x_amr,  &
!$OMP&    EntFlux_x_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice
          ixr = ntp_zones      ! location of right slice

! ----------

          Flux_x_amr     ( ixr, iyy, izz, 1: n_cc_components, ii) =
     1    Flux_snd_r_amr (      iyy, izz, 1: n_cc_components, ii)

          Flux_x_amr     ( ixl, iyy, izz, 1: n_cc_components, ii) =
     1    Flux_snd_l_amr (      iyy, izz, 1: n_cc_components, ii)

! ----------

          EntFlux_x_amr  ( ixr, iyy, izz, ii) =
     1    Flux_snd_r_amr (      iyy, izz, n_cc_components + 1, ii)

          EntFlux_x_amr  ( ixl, iyy, izz, ii) =
     1    Flux_snd_l_amr (      iyy, izz, n_cc_components + 1, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Efy_amr        ( ixr, iyy, izz, ii) =
     1    Flux_snd_r_amr (      iyy, izz, n_cc_components + 2, ii)

          Efy_amr        ( ixl, iyy, izz, ii) =
     1    Flux_snd_l_amr (      iyy, izz, n_cc_components + 2, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Efz_amr        ( ixr, iyy, izz, ii) =
     1    Flux_snd_r_amr (      iyy, izz, n_cc_components + 3, ii)

          Efz_amr        ( ixl, iyy, izz, ii) =
     1    Flux_snd_l_amr (      iyy, izz, n_cc_components + 3, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "Flux_snd_u_amr, Flux_snd_d_amr". First for the face centered fluxes,
! "Flux_y_amr, EntFlux_y_amr" then for the edge centered variables,
! "Efz_amr, Efx_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_y_amr, EntFlux_y_amr" data in
! "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Flux_y_amr,  &
!$OMP&    EntFlux_y_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice
          iyu = ntp_zones      ! location of up   slice

! ----------

          Flux_y_amr     ( ixx, iyu, izz, 1: n_cc_components, ii) =
     1    Flux_snd_u_amr ( ixx,      izz, 1: n_cc_components, ii)

          Flux_y_amr     ( ixx, iyd, izz, 1: n_cc_components, ii) = 
     1    Flux_snd_d_amr ( ixx,      izz, 1: n_cc_components, ii)

! ----------

          EntFlux_y_amr  ( ixx, iyu, izz, ii) =
     1    Flux_snd_u_amr ( ixx,      izz, n_cc_components + 1, ii)

          EntFlux_y_amr  ( ixx, iyd, izz, ii) =
     1    Flux_snd_d_amr ( ixx,      izz, n_cc_components + 1, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Efz_amr        ( ixx, iyu, izz, ii) =
     1    Flux_snd_u_amr ( ixx,      izz, n_cc_components + 2, ii)

          Efz_amr        ( ixx, iyd, izz, ii) =
     1    Flux_snd_d_amr ( ixx,      izz, n_cc_components + 2, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Efx_amr        ( ixx, iyu, izz, ii) =
     1    Flux_snd_u_amr ( ixx,      izz, n_cc_components + 3, ii)

          Efx_amr        ( ixx, iyd, izz, ii) =
     1    Flux_snd_d_amr ( ixx,      izz, n_cc_components + 3, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! Store top and bottom z--data in "Flux_snd_t_amr, Flux_snd_b_amr".


      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_t_amr, Flux_snd_b_amr". First for the face centered fluxes,
! "Flux_z_amr, EntFlux_z_amr" then for the edge centered variables,
! "Efx_amr, Efy_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_z_amr, EntFlux_z_amr" data in
! "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Flux_z_amr,  &
!$OMP&    EntFlux_z_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

! ----------

          Flux_z_amr     ( ixx, iyy, izt, 1: n_cc_components, ii) =
     1    Flux_snd_t_amr ( ixx, iyy,      1: n_cc_components, ii)

          Flux_z_amr     ( ixx, iyy, izb, 1: n_cc_components, ii) =
     1    Flux_snd_b_amr ( ixx, iyy,      1: n_cc_components, ii)

! ----------

          EntFlux_z_amr  ( ixx, iyy, izt, ii) =
     1    Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 1, ii)

          EntFlux_z_amr  ( ixx, iyy, izb, ii) =
     1    Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 1, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Efx_amr        ( ixx, iyy, izt, ii) =
     1    Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 2, ii)

          Efx_amr        ( ixx, iyy, izb, ii) =
     1    Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 2, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Efy_amr        ( ixx, iyy, izt, ii) =
     1    Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 3, ii)

          Efy_amr        ( ixx, iyy, izb, ii) =
     1    Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 3, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Assign Defect fluxes and electric fields at the coinciding faces and edges.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz )                             &

!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components,                                              &
!$OMP&    Defect_Flux_x_amr, Flux_x_amr, Defect_Efx_amr, Efx_amr,       &
!$OMP&    Defect_Flux_y_amr, Flux_y_amr, Defect_Efy_amr, Efy_amr,       &
!$OMP&    Defect_Flux_z_amr, Flux_z_amr, Defect_Efz_amr, Efz_amr,       &
!$OMP&    Defect_EntFlux_x_amr, EntFlux_x_amr,                          &
!$OMP&    Defect_EntFlux_y_amr, EntFlux_y_amr,                          &
!$OMP&    Defect_EntFlux_z_amr, EntFlux_z_amr)                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0,  nr_zones
        DO iyy = 0, ntp_zones
        DO ixx = 0, ntp_zones

! --------------------
! x-face.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_x_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1       - Flux_x_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_x_amr ( ixx, iyy, izz, ii) =
     1       - EntFlux_x_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "x-face".

! --------------------
! y-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_y_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1       - Flux_y_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_y_amr ( ixx, iyy, izz, ii) =
     1       - EntFlux_y_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "y-face".

! --------------------
! z-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_z_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1       - Flux_z_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1       - EntFlux_z_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "z-face".

! --------------------
! x-edge.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efx_amr ( ixx, iyy, izz, ii) =
     1       - Efx_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "x-edge".

! --------------------
! y-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efy_amr ( ixx, iyy, izz, ii) =
     1       - Efy_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "y-edge".

! --------------------
! z-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efz_amr ( ixx, iyy, izz, ii) =
     1       - Efz_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "z-edge".

! --------------------

        END DO ! End of "DO ixx = 0, ntp_zones".
        END DO ! End of "DO iyy = 0, ntp_zones".
        END DO ! End of "DO izz = 0,  nr_zones".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_AMR_FLUXES_AT_BDY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_XY ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the build of fluxes and electric fields,
! so that the fluxes and electric fields are synchronized at the coinciding
! faces.  The date is stored in the big container arrays, "Flux_x_amr, 
! Flux_y_amr, EntFlux_x_amr, EntFlux_y_amr, Efx_amr, Efy_amr, Efz_amr". 
! This data is loaded into 
! "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr, Flux_snd_d_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

! ----------------------------------------------------------------------
! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_r_amr ( :, :, :, :) = 0.0
      Flux_snd_l_amr ( :, :, :, :) = 0.0
      Flux_snd_u_amr ( :, :, :, :) = 0.0
      Flux_snd_d_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "Flux_snd_r_amr, Flux_snd_l_amr".
! Then store up and down y-data in "Flux_snd_u_amr, Flux_snd_d_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_amr, Flux_snd_l_amr". First for the face centered fluxes,
! "Flux_x_amr, EntFlux_x_amr" then for the edge centered variables,
! "Efy_amr, Efz_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_x_amr, EntFlux_x_amr" data in
! "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Flux_x_amr,  &
!$OMP&    EntFlux_x_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice
          ixr = ntp_zones      ! location of right slice

! ----------

          Flux_snd_r_amr (      iyy, izz, 1: n_cc_components, ii) =
     1        Flux_x_amr ( ixr, iyy, izz, 1: n_cc_components, ii) 

          Flux_snd_l_amr (      iyy, izz, 1: n_cc_components, ii) =
     1        Flux_x_amr ( ixl, iyy, izz, 1: n_cc_components, ii) 

! ----------

          Flux_snd_r_amr (      iyy, izz, n_cc_components + 1, ii) =
     1     EntFlux_x_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, n_cc_components + 1, ii) =
     1     EntFlux_x_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Flux_snd_r_amr (      iyy, izz, n_cc_components + 2, ii) =
     1           Efy_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, n_cc_components + 2, ii) =
     1           Efy_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Flux_snd_r_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "Flux_snd_u_amr, Flux_snd_d_amr". First for the face centered fluxes,
! "Flux_y_amr, EntFlux_y_amr" then for the edge centered variables,
! "Efz_amr, Efx_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_y_amr, EntFlux_y_amr" data in
! "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Flux_y_amr,  &
!$OMP&    EntFlux_y_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice
          iyu = ntp_zones      ! location of up   slice

! ----------

          Flux_snd_u_amr ( ixx,      izz, 1: n_cc_components, ii) =
     1        Flux_y_amr ( ixx, iyu, izz, 1: n_cc_components, ii) 

          Flux_snd_d_amr ( ixx,      izz, 1: n_cc_components, ii) =
     1        Flux_y_amr ( ixx, iyd, izz, 1: n_cc_components, ii) 

! ----------

          Flux_snd_u_amr ( ixx,      izz, n_cc_components + 1, ii) =
     1     EntFlux_y_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, n_cc_components + 1, ii) =
     1     EntFlux_y_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Flux_snd_u_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Flux_snd_u_amr ( ixx,      izz, n_cc_components + 3, ii) =
     1           Efx_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, n_cc_components + 3, ii) =
     1           Efx_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_Z ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the build of fluxes and electric fields,
! so that the fluxes and electric fields are synchronized at the coinciding
! faces.  The date is stored in the big container arrays, "Flux_z_amr, 
! EntFlux_z_amr, Efx_amr, Efy_amr". This data is loaded into
! "Flux_snd_t_amr, Flux_snd_b_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        izt, izb, izs

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_t_amr ( :, :, :, :) = 0.0
      Flux_snd_b_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! Store top and bottom z--data in "Flux_snd_t_amr, Flux_snd_b_amr".


      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_t_amr, Flux_snd_b_amr". First for the face centered fluxes,
! "Flux_z_amr, EntFlux_z_amr" then for the edge centered variables,
! "Efx_amr, Efy_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_z_amr, EntFlux_z_amr" data in
! "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Flux_z_amr,  &
!$OMP&    EntFlux_z_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

! ----------

          Flux_snd_t_amr ( ixx, iyy,      1: n_cc_components, ii) =
     1        Flux_z_amr ( ixx, iyy, izt, 1: n_cc_components, ii)

          Flux_snd_b_amr ( ixx, iyy,      1: n_cc_components, ii) =
     1        Flux_z_amr ( ixx, iyy, izb, 1: n_cc_components, ii)

! ----------

          Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 1, ii) =
     1     EntFlux_z_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 1, ii) =
     1     EntFlux_z_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 2, ii) =
     1           Efx_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 2, ii) =
     1           Efx_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_amr ( ixx, iyy,      n_cc_components + 3, ii) =
     1           Efy_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      n_cc_components + 3, ii) =
     1           Efy_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_UX ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its right, it does one-sided get of "U_snd_l" from that patch.
! For the patch to its left,  it does one-sided get of "U_snd_r" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its right, it does one-sided get of "U_snd_u" from that patch.
! For the patch to its left,  it does one-sided get of "U_snd_l" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn, ixxp1, iyyp1, ixx0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_xface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "U_recv_, U_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the right sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! right sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_l ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_d"
! is unpacked into "U_recev_r" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "U_recev_r", we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! right sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access right sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for right sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 2)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 2)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 2)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_l ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_r".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_l ( :, :, :, :, yr_patch_in_proc)" to "U_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_r, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy, ! remote data on yr.
     1               MPI_REAL8, Data_U_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_d".
! Index flipping is done so that "U_recev_d" is unpacked into "U_recev_r"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_u ( :, :, :, :, yr_patch_in_proc)" to "U_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_d, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_d" to "U_recev_r".
! Watch out, because "U_recev_d" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

        ixc = 1 + NGhost        ! begin index + end index.

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.
        iyn = ixc - ixx         ! yr_sector upper-y  is  my_sector lower-x.

        U_recev_r ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_d ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 0, NGhost

         ixc = 0 + NGhost       ! begin index + end index.

         izn = izz
         ixn = iyy              ! yr_sector - x  is  my_sector - y.
         iyn = ixc - ixx        ! yr_sector upper-y  is  my_sector lower-x.

         U_recev_r ( ixx, iyy, izz, n_cc_components + 2) =
     1 - U_recev_d ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 0, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (+ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

        ixc = 1 + NGhost        ! begin index + end index.

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.
        iyn = ixc - ixx         ! yr_sector upper-y  is  my_sector lower-x.

        U_recev_r ( ixx, iyy, izz, n_cc_components + 3) =
     1  U_recev_d ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

        ixc = 1 + NGhost        ! begin index + end index.

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.
        iyn = ixc - ixx         ! yr_sector upper-y  is  my_sector lower-x.

        U_recev_r ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_d ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "U_recev_r" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 + ntp_zones, ntp_zones + NGhost

! ----------

        ixc = ixx - ntp_zones     ! location in "U_recev_".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_r ( ixc, iyy, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_r     ( ixc, iyy, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_r" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 0 + ntp_zones, ntp_zones + NGhost

        ixc = ixx - ntp_zones     ! location in "U_recev_".

! ----------
! Exclude any copying or averaging at x- faces abutting tri-corners.
        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( ixx .EQ. ntp_zones) .AND.   ! Active Coinciding face.     
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN 

          Bxf_amr   ( ixx, iyy, izz, ii) =
     1  ( Bxf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_r ( ixc, iyy, izz, n_cc_components + 2)) * 0.5

        ELSE

          Bxf_amr   ( ixx, iyy, izz, ii) =                               
     1    U_recev_r ( ixc, iyy, izz, n_cc_components + 2) 

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 0 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_r" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost
      DO ixx = 1 + ntp_zones, ntp_zones + NGhost

        ixc = ixx - ntp_zones     ! location in "U_recev_".

        Byf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_r ( ixc, iyy, izz, n_cc_components + 3)

      END DO ! End of "DO ixx = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_r" to "Bzf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 + ntp_zones, ntp_zones + NGhost

        ixc = ixx - ntp_zones     ! location in "U_recev_".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_r ( ixc, iyy, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the left sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! left sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_l ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_r"
! is unpacked into "U_recev_l" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "U_recev_l", we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! left sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector)" will allow us to access left sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for left sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 3)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 3)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 3)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_l".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_r ( :, :, :, :, yr_patch_in_proc)" to "U_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_l, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_l ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_r".
! Index flipping is done so that "U_recev_r" is unpacked into "U_recev_l"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_l( :, :, :, :, yr_patch_in_proc)" to "U_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_r, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_r" to "U_recev_l".
! Watch out, because "U_recev_r" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

        ixc = 1 + NGhost                       ! begin index + end index.
        iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_l ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_r ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 0, NGhost

          ixc = 0 + NGhost                       ! begin index + end index.
          iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

          izn = izz
          ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
          iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

          U_recev_l ( ixx, iyy, izz, n_cc_components + 2) =
     1  - U_recev_r ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 0, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

          ixc = 1 + NGhost                       ! begin index + end index.
          iyc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

          izn = izz
          ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
          iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

          U_recev_l ( ixx, iyy, izz, n_cc_components + 3) =
     1  - U_recev_r ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1, NGhost

        ixc = 1 + NGhost                       ! begin index + end index.
        iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_l ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_r ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1, NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "U_recev_l" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, 0

! ----------

        ixc = ixx + NGhost        ! location in "U_recev_l".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_l ( ixc, iyy, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_l     ( ixc, iyy, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, 0".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_l" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 0 - NGhost, 0

        ixc = ixx + NGhost             ! location in "U_recev_l".

! ----------
! Exclude any copying or averaging at x- faces abutting tri-corners.
        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( ixx .EQ. 0) .AND.          ! Active Coinciding face.     
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Bxf_amr   ( ixx, iyy, izz, ii) =                               
     1  ( Bxf_amr   ( ixx, iyy, izz, ii)                                 
     1  + U_recev_l ( ixc, iyy, izz, n_cc_components + 2)) * 0.5

        ELSE                             ! Ghost face

          Bxf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_l ( ixc, iyy, izz, n_cc_components + 2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, 0".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_l" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, 0

! ----------

        ixc = ixx + NGhost        ! location in "U_recev_l".

        Byf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_l ( ixc, iyy, izz, n_cc_components + 3)

      END DO ! End of "DO ixx = 1 - NGhost, 0".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_l" to "Bzf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, 0

        ixc = ixx + NGhost        ! location in "U_recev_l".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_l ( ixc, iyy, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, 0".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_UX

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_UY ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "U_snd_d" from that patch.
! For the patch to its down, it does one-sided get of "U_snd_u" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "U_snd_r" from that patch.
! For the patch to its down, it does one-sided get of "U_snd_d" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "U_recv_, U_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_l"
! is unpacked into "U_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "U_recev_u we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc
! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_d ( :, :, :, :, yr_patch_in_proc)" to "U_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_u, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_l".
! Index flipping is done so that "U_recev_l" is unpacked into "U_recev_u"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_r ( :, :, :, :, yr_patch_in_proc)" to "U_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_l, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_l" to "U_recev_u".
! Watch out, because "U_recev_l" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_l ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (+ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 2) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

          iyc = 0 + NGhost        ! begin index + end index.

          izn = izz
          iyn = ixx               ! yr_sector - y  is  my_sector - x.
          ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

          U_recev_u ( ixx, iyy, izz, n_cc_components + 3) =
     1  - U_recev_l ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "U_recev_u" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

! ----------

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_u ( ixx, iyc, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_u     ( ixx, iyc, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bxf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

! ----------
! Exclude any copying or averaging at y- faces abutting tri-corners.
        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( iyy .EQ. ntp_zones) .AND.    ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Byf_amr   ( ixx, iyy, izz, ii) =
     1  ( Byf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_u ( ixx, iyc, izz, n_cc_components + 3)) * 0.5

        ELSE                                ! Ghost face.

          Byf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_u ( ixx, iyc, izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the down sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! down sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_u"
! is unpacked into "U_recev_d" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "U_recev_d", we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! down sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector)" will allow us to access down  sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for down sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 5)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 5)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 5)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_u ( :, :, :, :, yr_patch_in_proc)" to "U_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_d, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_u".
! Index flipping is done so that "U_recev_u" is unpacked into "U_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_d ( :, :, :, :, yr_patch_in_proc)" to "U_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_u, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_u" to "U_recev_d".
! Watch out, because "U_recev_u" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
        iyc = 1 + NGhost                       ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_d ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_u ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

         ixc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.
         iyc = 1 + NGhost                       ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         U_recev_d ( ixx, iyy, izz, n_cc_components + 2) =
     1 - U_recev_u ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
         iyc = 0 + NGhost                       ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         U_recev_d ( ixx, iyy, izz, n_cc_components + 3) =
     1 - U_recev_u ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
        iyc = 1 + NGhost                       ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_d ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_u ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "U_recev_d" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------

        iyc = iyy + NGhost        ! location in "U_recev_d".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_d ( ixx, iyc, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_d     ( ixx, iyc, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_d" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

        Bxf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_d ( ixx, iyc, izz, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

! ----------
! Exclude any copying or averaging at y- faces abutting tri-corners.
        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( iyy .EQ. 0) .AND.         ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Byf_amr   ( ixx, iyy, izz, ii) =
     1  ( Byf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_d ( ixx, iyc, izz, n_cc_components + 3)) * 0.5

        ELSE                             ! Ghost face.

          Byf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_d ( ixx, iyc, izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_d ( ixx, iyc, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_UY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_UY_TOP ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "U_snd_d" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "U_snd_r" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z 

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "U_recv_, U_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_l"
! is unpacked into "U_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "U_recev_u we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc
! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_d ( :, :, :, :, yr_patch_in_proc)" to "U_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_u, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_l".
! Index flipping is done so that "U_recev_l" is unpacked into "U_recev_u"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_r ( :, :, :, :, yr_patch_in_proc)" to "U_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_l, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_l" to "U_recev_u".
! Watch out, because "U_recev_l" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_l ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (+ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 2) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

          iyc = 0 + NGhost        ! begin index + end index.

          izn = izz
          iyn = ixx               ! yr_sector - y  is  my_sector - x.
          ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

          U_recev_u ( ixx, iyy, izz, n_cc_components + 3) =
     1  - U_recev_l ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "U_recev_u" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

! ----------

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_u ( ixx, iyc, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_u     ( ixx, iyc, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bxf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

! ----------
! Exclude any copying or averaging at y- faces abutting tri-corners.
        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( iyy .EQ. ntp_zones) .AND.    ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Byf_amr   ( ixx, iyy, izz, ii) =
     1  ( Byf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_u ( ixx, iyc, izz, n_cc_components + 3)) * 0.5

        ELSE                                ! Ghost face.

          Byf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_u ( ixx, iyc, izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_UY_TOP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_UZ ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the z-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank,
     1        ix_patch, iy_patch, iz_patch

      INTEGER num_of_reals_xy, num_of_reals_z

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "U_recv_, U_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! top sibling patch.

! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of 
! "U_snd_b ( :, :, :, :, yr_patch_in_proc)", just as we would do 
! in a Cartesian mesh code. This is saved in "U_recev_t" on the 
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) Now that the data is in "U_recev_t" we can unpack it into 
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! top sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access top sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for top sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 6)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 6)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 6)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "U_snd_b ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_t".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_b ( :, :, :, :, yr_patch_in_proc)" to "U_recev_t".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_ZB_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_t, num_of_reals_z, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,   ! remote data on yr.
     1               MPI_REAL8, Data_U_ZB_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_ZB_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) Unpack data from "U_recev_t" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 + nr_zones, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

! ----------

        izc = izz - nr_zones     ! location in "U_recev_t".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_t ( ixx, iyy, izc, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_t     ( ixx, iyy, izc, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO izz = 1 + nr_zones, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_t" to "Bxf_amr".

      DO izz = 1 + nr_zones, nr_zones  + NGhost                          
      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 0 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz - nr_zones     ! location in "U_recev_t".              
                                                                         
        Bxf_amr   ( ixx, iyy, izz, ii) =                             
     1  U_recev_t ( ixx, iyy, izc, n_cc_components + 2)  
                                                                         
      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 1 + nr_zones, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_t" to "Byf_amr".

      DO izz = 1 + nr_zones, nr_zones  + NGhost                          
      DO iyy = 0 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz - nr_zones     ! location in "U_recev_t".              
                                                                         
        Byf_amr   ( ixx, iyy, izz, ii) =                             
     1  U_recev_t ( ixx, iyy, izc, n_cc_components + 3)  
                                                                         
      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 1 + nr_zones, nr_zones  + NGhost".

! --------------------                                                   
! Unpack data from "U_recev_t" to "Bzf_amr".                             

      DO izz = 0 + nr_zones, nr_zones  + NGhost                          
      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz - nr_zones     ! location in "U_recev_t". 

! ----------                                                             
                                                                         
        IF ( ( izz .EQ. nr_zones) .AND.   ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Bzf_amr   ( ixx, iyy, izz, ii) =
     1  ( Bzf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_t ( ixx, iyy, izc, n_cc_components + 4)) * 0.5

        ELSE                              ! Ghost face.

          Bzf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_t ( ixx, iyy, izc, n_cc_components + 4)

        END IF

! ----------                                                             
                                                                         
      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 0 + nr_zones, nr_zones  + NGhost".  

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! bottom sibling patch.

! STEP 2.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of 
! "U_snd_t ( :, :, :, :, yr_patch_in_proc)", just as we would do 
! in a Cartesian mesh code. This is saved in "U_recev_b" on the 
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) Now that the data is in "U_recev_b" we can unpack it into 
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! bottom sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access bottom sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for bottom sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 7)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 7)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 7)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 2.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "U_snd_t ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_b".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_t ( :, :, :, :, yr_patch_in_proc)" to "U_recev_b".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_ZT_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_b, num_of_reals_z, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,   ! remote data on yr.
     1               MPI_REAL8, Data_U_ZT_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_ZT_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) Unpack data from "U_recev_b" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, 0
      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

! ----------

        izc = izz + NGhost     ! location in "U_recev_b".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_b ( ixx, iyy, izc, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_b     ( ixx, iyy, izc, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, 0".

! --------------------
! Unpack data from "U_recev_b" to "Bxf_amr".

      DO izz = 1 - NGhost, 0 
      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 0 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz + NGhost     ! location in "U_recev_b".              
                                                                         
        Bxf_amr   ( ixx, iyy, izz, ii) =                             
     1  U_recev_b ( ixx, iyy, izc, n_cc_components + 2)  
                                                                         
      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 1 - NGhost, 0".

! --------------------
! Unpack data from "U_recev_b" to "Byf_amr".

      DO izz = 1 - NGhost, 0
      DO iyy = 0 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz + NGhost     ! location in "U_recev_b".              
                                                                         
        Byf_amr   ( ixx, iyy, izz, ii) =                             
     1  U_recev_b ( ixx, iyy, izc, n_cc_components + 3)  
                                                                         
      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 1 - NGhost, 0".

! --------------------                                                   
! Unpack data from "U_recev_b" to "Bzf_amr".                             

      DO izz = 0 - NGhost, 0
      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost                           
                                                                         
        izc = izz + NGhost     ! location in "U_recev_b". 

! ----------                                                             
                                                                         
        IF ( ( izz .EQ. 0) .AND.   ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Bzf_amr   ( ixx, iyy, izz, ii) =
     1  ( Bzf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_b ( ixx, iyy, izc, n_cc_components + 4)) * 0.5

        ELSE                       ! Ghost face

          Bzf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_b ( ixx, iyy, izc, n_cc_components + 4)

        END IF

! ----------                                                             
                                                                         
      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO izz = 0 - NGhost, 0".  

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_UZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_X ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its right, it does one-sided get of "Flux_snd_l" from that patch.
! For the patch to its left,  it does one-sided get of "Flux_snd_r" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its right, it does one-sided get of "Flux_snd_u" from that patch.
! For the patch to its left,  it does one-sided get of "Flux_snd_l" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, iyc,
     1        ixn, iyn, izn, ixxp1, iyyp1, ixx0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z 

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_xface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) 
     1                * ( n_cc_components + 3 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the right sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! right sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_d"
! is unpacked into "Flux_recev_r" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "Flux_recev_r", we can unpack it into "Flux_x_amr, EntFlux_x_amr, 
! Efy_amr, Efz_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! right sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access right sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for right sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 2)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 2)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 2)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_r".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_l ( :, :, :, yr_patch_in_proc)" to "Flux_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_r, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_d".
! Index flipping is done so that "Flux_recev_d" is unpacked into "Flux_recev_r"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_u ( :, :, :, yr_patch_in_proc)" to "Flux_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_d, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_d" to "Flux_recev_r".
! Watch out, because "Flux_recev_d" does not have data in normal ordering.

!  my_sector x-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         izn = izz
         ixn = iyy              ! yr_sector - x  is  my_sector - y.

         Flux_recev_r ( iyy, izz, 1: n_cc_components + 1) =
     1 - Flux_recev_d ( ixn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-edge-centered data "Efy_amr".
!  my_sector y-edge == (+ve) yr_sector x-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.

        Flux_recev_r ( iyy, izz, n_cc_components + 2) =    ! YZ: Y at +2
     1  Flux_recev_d ( ixn, izn, n_cc_components + 3)      ! ZX: X at +3

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.

        Flux_recev_r ( iyy, izz, n_cc_components + 3) =    ! YZ: Z at +3
     1  Flux_recev_d ( ixn, izn, n_cc_components + 2)      ! ZX: Z at +2

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "Flux_recev_r" to
! "Flux_x_amr, EntFlux_x_amr, Efy_amr, Efz_amr".

! --------------------
! Unpack data from "Flux_recev_r" to "Flux_x_amr, EntFlux_x_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding face. 
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_r (      iyy, izz, 1: n_cc_components)) * 0.5

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_x_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_r (      iyy, izz, n_cc_components + 1)) * 0.5

        ELSE

          Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_r (      iyy, izz, 1: n_cc_components)

          EntFlux_x_amr( ixx, iyy, izz, ii) = 
     1    Flux_recev_r (      iyy, izz, n_cc_components + 1)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_r" to "Efy_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efy_amr ( ixx, iyy, izz, ii) =
     1  ( Efy_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_r ( iyy, izz, n_cc_components + 2)) * 0.5

        ELSE

          Efy_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_r ( iyy, izz, n_cc_components + 2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_r" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          IF ( ( ( iyy .EQ. 0 ) .AND. ( itri_patch_RD .EQ. 1)) .OR.
     1         ( ( iyy .EQ. ntp_zones) .AND.
     1                                ( itri_patch_RU .EQ. 1)) ) THEN

            Efz_amr ( ixx, iyy, izz, ii) =   ! accumulate tri-corner Efz; * 1/3 later.
     1      Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_r ( iyy, izz, n_cc_components + 3)

          ELSE

            Efz_amr ( ixx, iyy, izz, ii) =
     1    ( Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_r ( iyy, izz, n_cc_components + 3)) * 0.5

          END IF

        ELSE

          Efz_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_r ( iyy, izz, n_cc_components + 3)

        END IF

! ----------

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the left sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! left sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_r"
! is unpacked into "Flux_recev_l" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "Flux_recev_l", we can unpack it into "Flux_x_amr, EntFlux_x_amr, 
! Efy_amr, Efz_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! left sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access left sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for left sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 3)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 3)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 3)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_l".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_l, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_r".
! Index flipping is done so that "Flux_recev_r" is unpacked into "Flux_recev_l"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_l ( :, :, :, yr_patch_in_proc)" to "Flux_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_r, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_r" to "Flux_recev_l".
! Watch out, because "Flux_recev_r" does not have data in normal ordering.

!  my_sector x-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         Flux_recev_l ( iyy, izz, 1: n_cc_components + 1) =
     1 - Flux_recev_r ( iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-edge-centered data "Efy_amr".
!  my_sector y-edge == (-ve) yr_sector y-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.
         izn = izz

         Flux_recev_l ( iyy, izz, n_cc_components + 2) =    ! YZ: Y at +2
     1 - Flux_recev_r ( iyn, izn, n_cc_components + 2)      ! YZ: Y at +2

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        iyc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

        izn = izz
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        Flux_recev_l ( iyy, izz, n_cc_components + 3) =    ! YZ: Z at +3
     1  Flux_recev_r ( iyn, izn, n_cc_components + 3)      ! YZ: Z at +3

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "Flux_recev_l" to
! "Flux_x_amr, EntFlux_x_amr, Efy_amr, Efz_amr".

! --------------------
! Unpack data from "Flux_recev_l" to "Flux_x_amr, EntFlux_x_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding face. 
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_l (      iyy, izz, 1: n_cc_components)) * 0.5

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_x_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_l (      iyy, izz, n_cc_components + 1)) * 0.5

        ELSE

          Flux_x_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_l (      iyy, izz, 1: n_cc_components)

          EntFlux_x_amr( ixx, iyy, izz, ii) = 
     1    Flux_recev_l (      iyy, izz, n_cc_components + 1)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_l" to "Efy_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efy_amr ( ixx, iyy, izz, ii) =
     1  ( Efy_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_l ( iyy, izz, n_cc_components + 2)) * 0.5

        ELSE

          Efy_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_l ( iyy, izz, n_cc_components + 2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_l" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones) .AND. ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          IF ( ( ( iyy .EQ. 0 ) .AND. ( itri_patch_LD .EQ. 1)) .OR.
     1         ( ( iyy .EQ. ntp_zones) .AND. 
     1                                ( itri_patch_LU .EQ. 1)) ) THEN

            Efz_amr ( ixx, iyy, izz, ii) = ! accumulate tri-corner Efz; * 1/3 later.
     1      Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_l ( iyy, izz, n_cc_components + 3)

          ELSE

            Efz_amr ( ixx, iyy, izz, ii) =
     1    ( Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_l ( iyy, izz, n_cc_components + 3)) * 0.5

          END IF

        ELSE

          Efz_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_l ( iyy, izz, n_cc_components + 3)

        END IF

! ----------

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_d" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_u" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_r" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_d" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, 
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) 
     1                * ( n_cc_components + 3 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_l"
! is unpacked into "Flux_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "Flux_recev_u", we can unpack it into "Flux_y_amr, EntFlux_y_amr, 
! Efz_amr, Efx_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_l".
! Index flipping is done so that "Flux_recev_l" is unpacked into "Flux_recev_u"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_l, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_l" to "Flux_recev_u".
! Watch out, because "Flux_recev_l" does not have data in normal ordering.

!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         izn = izz
         iyn = ixx               ! yr_sector - y  is  my_sector - x.

         Flux_recev_u ( ixx, izz, 1: n_cc_components + 1) =
     1 - Flux_recev_l ( iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, n_cc_components + 2) =    ! ZX: Z at +2
     1  Flux_recev_l ( iyn, izn, n_cc_components + 3)      ! YZ: Z at +3

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-edge-centered data "Efx_amr".
!  my_sector x-edge == (+ve) yr_sector y-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, n_cc_components + 3) =    ! ZX: X at +3
     1  Flux_recev_l ( iyn, izn, n_cc_components + 2)      ! YZ: Y at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "Flux_recev_u" to
! "Flux_y_amr, EntFlux_y_amr, Efz_amr, Efx_amr".

! --------------------
! Unpack data from "Flux_recev_u" to "Flux_y_amr, EntFlux_y_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_u ( ixx,      izz, 1: n_cc_components)) * 0.5

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_y_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 1)) * 0.5

        ELSE

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_u ( ixx,      izz, 1: n_cc_components)

          EntFlux_y_amr( ixx, iyy, izz, ii) = 
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 1)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

! -----
! accumulate tri-corner Efz; * 1/3 later.

          IF ( ( ( ixx .EQ. 0) .AND. ( itri_patch_LU .EQ. 1)) .OR.
     1         ( ( ixx .EQ. ntp_zones) .AND.
     1                               ( itri_patch_RU .EQ. 1))) THEN

            Efz_amr      ( ixx, iyy, izz, ii) = 
     1      Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_u ( ixx,      izz, n_cc_components + 2)

! -----
! Exclude quad corner on "Efz" averaging for y-flux; this is handled separately.

          ELSE IF ( ( ixx .NE. 0) .AND. ( ixx .NE. ntp_zones)) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1    ( Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_u ( ixx,      izz, n_cc_components + 2)) * 0.5

! -----

          END IF

! ----------

        ELSE

! ----------

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 2)

! ----------

        END IF  ! IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) ...

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efx_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 3)) * 0.5

        ELSE

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the down sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! down sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_u"
! is unpacked into "Flux_recev_d" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "Flux_recev_d", we can unpack it into "Flux_y_amr, EntFlux_y_amr, 
! Efz_amr, Efx_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! down sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access down sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for down sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 5)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 5)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 5)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_u ( :, :, :, yr_patch_in_proc)" to "Flux_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_d, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_u".
! Index flipping is done so that "U_recev_u" is unpacked into "Flux_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_u" to "Flux_recev_d".
! Watch out, because "Flux_recev_u" does not have data in normal ordering.

!  my_sector y-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, 1: n_cc_components + 1) =
     1 - Flux_recev_u ( ixn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

         ixc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, n_cc_components + 2) =    ! ZX: Z at +2
     1   Flux_recev_u ( ixn, izn, n_cc_components + 2)      ! ZX: Z at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-edge-centered data "Efx_amr".
!  my_sector x-edge == (-ve) yr_sector x-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, n_cc_components + 3) =    ! ZX: X at +3
     1 - Flux_recev_u ( ixn, izn, n_cc_components + 3)      ! ZX: X at +3

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "Flux_recev_d" to
! "Flux_y_amr, EntFlux_y_amr, Efz_amr, Efx_amr".

! --------------------
! Unpack data from "Flux_recev_d" to "Flux_y_amr, EntFlux_y_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.   ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_d ( ixx,      izz, 1: n_cc_components)) * 0.5

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_y_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, n_cc_components + 1)) * 0.5

        ELSE

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_d ( ixx,      izz, 1: n_cc_components)

          EntFlux_y_amr( ixx, iyy, izz, ii) = 
     1    Flux_recev_d ( ixx,      izz, n_cc_components + 1)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_d" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

! -----
! accumulate tri-corner Efz; * 1/3 later.

          IF ( ( ( ixx .EQ. 0) .AND. ( itri_patch_LD .EQ. 1)) .OR. 
     1         ( ( ixx .EQ. ntp_zones) .AND.
     1                               ( itri_patch_RD .EQ. 1))) THEN

            Efz_amr      ( ixx, iyy, izz, ii) = 
     1      Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_d ( ixx,      izz, n_cc_components + 2)

! -----
! Exclude quad corner on "Efz" averaging for y-flux; this is handled separately.

          ELSE IF ( ( ixx .NE. 0) .AND. ( ixx .NE. ntp_zones)) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1    ( Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_d ( ixx,      izz, n_cc_components + 2)) * 0.5

! -----

          END IF

! ----------

        ELSE

! ----------

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_d ( ixx,      izz, n_cc_components + 2)

! ----------

        END IF  ! IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) ...

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_d" to "Efx_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, n_cc_components + 3)) * 0.5

        ELSE

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_d ( ixx,      izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! The tri-corner "Efz" are accumulated. Now correctly average the "Efz" 
! by dividing them by 3.

      DO ii = 1, npatch_per_proc

! ----------------------------------------

        itri_patch_LD = itri_patch_amr ( 1, ii)
        itri_patch_RD = itri_patch_amr ( 2, ii)
        itri_patch_RU = itri_patch_amr ( 3, ii)
        itri_patch_LU = itri_patch_amr ( 4, ii)

! ----------------------------------------

        DO izz = 1 - NGhost, nr_zones  + NGhost

! --------------------

          IF ( itri_patch_LD .EQ. 1)
     1    Efz_amr ( 0, 0, izz, ii) =
     1    Efz_amr ( 0, 0, izz, ii) / 3.0

          IF ( itri_patch_LU .EQ. 1)
     1    Efz_amr ( 0, ntp_zones, izz, ii) =
     1    Efz_amr ( 0, ntp_zones, izz, ii) / 3.0
        
          IF ( itri_patch_RD .EQ. 1)
     1    Efz_amr ( ntp_zones, 0, izz, ii) =
     1    Efz_amr ( ntp_zones, 0, izz, ii) / 3.0

          IF ( itri_patch_RU .EQ. 1)
     1    Efz_amr ( ntp_zones, ntp_zones, izz, ii) =
     1    Efz_amr ( ntp_zones, ntp_zones, izz, ii) / 3.0

! --------------------

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost"

! ----------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE SYNCHRONIZE_EFZ_QUAD_CORNER ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the synchronization of fluxes and electric
! fields at the coinciding faces. The "Efz" at the quad corner needs a
! further synchrnization so as to correctly average the "Efz" edge at the
! quad-corner.
!
! Consider the following quad-corner. The edge "Efz" is perpendicular to the screen.
!
!                            3  |  4
!                               |
!                        -------|-------
!                               |
!                            1  |  2
!                               
!
! The previous "FETCH..Flux_X ( )" subroutine synchronize at the x-face. 
! And it assigned the following new electric fields at the quad-corner:
! 
!            Efz_new_1 = ( Efz_1 + Efz_2 ) / 2
!            Efz_new_2 = ( Efz_1 + Efz_2 ) / 2
!
!            Efz_new_3 = ( Efz_3 + Efz_4 ) / 2
!            Efz_new_4 = ( Efz_3 + Efz_4 ) / 2
!
! Now the quad-corner electric fields are,
!
!        ( Efz_3 + Efz _4 ) / 2 |  ( Efz_3 + Efz _4 ) / 2
!                               |
!                        -------|-------
!                               |
!        ( Efz_1 + Efz _2 ) / 2 |  ( Efz_1 + Efz _2 ) / 2
!
!
! Note that the quad-corner Efz is not changed in the "FETCH..Flux_Y ( )" subroutine.
! This subroutine is called next, to handle the same.
!
! In this subroutine, the quad-corner Efz is synchronized in two distinct steps: 
!     i) LOAD  - Efz_new values.
!    ii) FETCH - synchronize along the y-face.
! 
! After the synchronization, the quad-corner Efz in all four shared patches,
! (1), (2), (3), (4) will have the value of:
!
!     "( Efz_1 + Efz _2 + Efz_3 + Efz_4 ) / 4"
!
! Thus completing the synchronization at the quad-corner Efz.
!
! --------------------
! The above mentioned technique works well in all the circumstances,
! except the for the "y-top" face when it interfaces with the next sector.
! This because of how the 'xy-axes' are defined in each sector.
!
! Hence, the Efz-quad corners which lie in the interfacing "y-top" face,
! are assigned separately; in step 5.
!
! ----------------------------------------------------------------------
!
! STEP 1 & 2: sync. "y-top"    quad corner Efz (y-top touching same sector patch)
! STEP 3    : sync. "y-bottom" quad corner Efz
! STEP 4 & 5: sync. "y-top"    quad corner Efz (y-top touching difft. sector patch)
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

      INTEGER ixc, iyc, ixn, iyn, izn, ixxp1, iyyp1, ixx0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_xface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------
! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_u_amr ( :, :, :, :) = 0.0
      Flux_snd_d_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! STEP 1) Loop over all the patches in this processor "me".
! Store up and down y-data in "Flux_snd_u_amr, Flux_snd_d_amr", only the "Efz"
! at the quad-corner.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - NGhost, nr_zones  + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

! ----------

          ixx = 0

          Flux_snd_u_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------

          ixx = ntp_zones

          Flux_snd_u_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, n_cc_components + 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO
! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Send complete; now synchronize and the recieve the Efz.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

! ----------------------------------------------------------------------
! STEP 2. Receive 'y-top' - "Flux_recev_u" if it is in the same sector.
! when "Flux_recev_u" lies in different sector, it has to be done
! separately at last; and it is done in step 5.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) Unpack data from "Flux_recev_u" to "Efz_amr" at quad corner.

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

        iyy = ntp_zones

! ----------
! Average "Efz" at the Left-Upper Quad corner.

        IF ( itri_patch_LU .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 2)) * 0.5

        ENDIF

! ----------
! Average "Efz" at the Right-Upper Quad corner.

        IF ( itri_patch_RU .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 2)) * 0.5

        ENDIF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) Receive 'y-bottom' from "Flux_recev_d". Unlike the 'y-top',
! here the values are received from the neighbouring sectors as well.
! Averaging of Efz from the neighbouring sector works for the 'y-bottom'
! Efz; due to how the xy-axes of neighbouring sectors defined.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 3.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! down sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access down sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for down sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 5)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 5)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 5)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------------------------------------
! STEP 3.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u_amr ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_u_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_d, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 3.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d_amr ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_u".
! Index flipping is done so that "Flux_recev_u" is unpacked into "Flux_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_u" to "Flux_recev_d".
! Watch out, because "Flux_recev_u" does not have data in normal ordering.

! Only the "Efz" is unpacked here.

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

         ixc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, n_cc_components + 2) =    ! ZX: Z at +2
     1   Flux_recev_u ( ixn, izn, n_cc_components + 2)      ! ZX: Z at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 3.d) Unpack data from "Flux_recev_d" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

      iyy = 0

! ----------
! Average "Efz" at the Left-Down Quad corner.

        IF ( itri_patch_LD .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, n_cc_components + 2)) * 0.5

        END IF

! ----------
! Average "Efz" at the Right-Down Quad corner.

        IF ( itri_patch_RD .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, n_cc_components + 2)) * 0.5

        END IF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4) Till here, the quad-corner "Efz" at the bottom-face is complete; 
! and the top-face is also completed for same-sector faces. 
! Now synchronize the top-face Efz; when it is shared with a different sector.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_r_amr ( :, :, :, :) = 0.0
      Flux_snd_l_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "Flux_snd_r_amr, Flux_snd_l_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_amr, Flux_snd_l_amr".
! Storing edge centered variable "Efz_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - NGhost, nr_zones  + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

! ----------

          iyy = 0

          Flux_snd_r_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------

          iyy = ntp_zones

          Flux_snd_r_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, n_cc_components + 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Send complete; now synchronize and the recieve the Efz.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! STEP 5. Receive 'y-top' - "Flux_recev_u" if it is coming from different
! sector. Note that the same sector synchronization for y-top is completed 
! in step 2.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 5.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------------------------------------
! STEP 5.b) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r_amr ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_l".
! Index flipping is done so that "Flux_recev_l" is unpacked into "Flux_recev_u"
! with the flipping.

      IF ( my_sector .NE. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_l, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_l" to "Flux_recev_u".
! Only the "Efz" data.
! Watch out, because "Flux_recev_l" does not have data in normal ordering.

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, n_cc_components + 2) =    ! ZX: Z at +2
     1  Flux_recev_l ( iyn, izn, n_cc_components + 3)      ! YZ: Z at +3

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------
! STEP 5.c) Unpack data from "Flux_recev_u" to "Efz_amr" at quad corner.

! NOTE: It is important to note here that the "Flux_recev_u" has 
! correctly averged "Efz" from all 4-contributions at the quadcorner.
! Hence, the values are --assigned-- here. Not averaged.

! ------------------------------------------------------------

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

        iyy = ntp_zones

! ----------
! Assign the average "Efz" at the Left-Upper Quad corner.

        IF ( itri_patch_LU .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 2)

        ENDIF

! ----------
! Assign the average "Efz" at the Right-Upper Quad corner.

        IF ( itri_patch_RU .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 2)

        ENDIF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .NE. yr_sector)".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE SYNCHRONIZE_EFZ_QUAD_CORNER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( 0 == 1)
! This subroutine is switched off. 
! For "U" we need such a subroutine to correctly get the ghost cell information.
! For "Flux, E-field" we need to average only on the active faces; and for that
! it is sufficient to call the above two "fetch" subroutines and then the
! efz-quad corner synchronization subroutine.
! Calling the following may cause float-point operation erros, when the source 
! terms are many magnitudes stronger than the actual evolving fields.


      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y_TOP ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_d" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_u" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_r" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_d" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, 
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) 
     1                * ( n_cc_components + 3 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_l"
! is unpacked into "Flux_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "Flux_recev_u", we can unpack it into "Flux_y_amr, EntFlux_y_amr, 
! Efz_amr, Efx_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_l".
! Index flipping is done so that "Flux_recev_l" is unpacked into "Flux_recev_u"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_l, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_l" to "Flux_recev_u".
! Watch out, because "Flux_recev_l" does not have data in normal ordering.

!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         izn = izz
         iyn = ixx               ! yr_sector - y  is  my_sector - x.

         Flux_recev_u ( ixx, izz, 1: n_cc_components + 1) =
     1 - Flux_recev_l ( iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, n_cc_components + 2) =    ! ZX: Z at +2
     1  Flux_recev_l ( iyn, izn, n_cc_components + 3)      ! YZ: Z at +3

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-edge-centered data "Efx_amr".
!  my_sector x-edge == (+ve) yr_sector y-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, n_cc_components + 3) =    ! ZX: X at +3
     1  Flux_recev_l ( iyn, izn, n_cc_components + 2)      ! YZ: Y at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "Flux_recev_u" to
! "Flux_y_amr, EntFlux_y_amr, Efz_amr, Efx_amr".

! --------------------
! Unpack data from "Flux_recev_u" to "Flux_y_amr, EntFlux_y_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_u ( ixx,      izz, 1: n_cc_components)) * 0.5

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_y_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 1)) * 0.5

        ELSE

          Flux_y_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_u ( ixx,      izz, 1: n_cc_components)

          EntFlux_y_amr( ixx, iyy, izz, ii) = 
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 1)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 2)) * 0.5

        ELSE

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 2)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efx_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, n_cc_components + 3)) * 0.5

        ELSE

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y_TOP

#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Z ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the z-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank,
     1        ix_patch, iy_patch, iz_patch

      INTEGER num_of_reals_xy, num_of_reals_z

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the top sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! top sibling patch.

! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of 
! "Flux_snd_b ( :, :, :, yr_patch_in_proc)", just as we would do 
! in a Cartesian mesh code. This is saved in "Flux_recev_t" on the 
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) Now that the data is in "Flux_recev_t" we can unpack it into 
! "Flux_z_amr, EntFlux_z_amr, Efx_amr, Efy_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! top sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access top sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for top sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 6)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 6)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 6)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "Flux_snd_b ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_t".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_b ( :, :, :, yr_patch_in_proc)" to "Flux_recev_t".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_ZB_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_t, num_of_reals_z, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_ZB_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_ZB_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) Unpack data from "Flux_recev_t" to
! "Flux_z_amr, EntFlux_z_amr, Efx_amr, Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1  ( Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_t ( ixx, iyy,      1: n_cc_components)) * 0.5

          EntFlux_z_amr ( ixx, iyy, izz, ii) = 
     1  ( EntFlux_z_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_t  ( ixx, iyy,      n_cc_components + 1)) * 0.5

        ELSE                              ! Ghost face.

          Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_t ( ixx, iyy,      1: n_cc_components)

          EntFlux_z_amr ( ixx, iyy, izz, ii) = 
     1    Flux_recev_t  ( ixx, iyy,      n_cc_components + 1)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_t" to "Efx_amr".

      DO iyy = 0 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost 

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones)) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_t ( ixx, iyy, n_cc_components + 2)) * 0.5

        ELSE                              ! Ghost face.

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_t ( ixx, iyy, n_cc_components + 2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost". 

! --------------------
! Unpack data from "Flux_recev_t" to "Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 0 - NGhost, ntp_zones  + NGhost 

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Efy_amr      ( ixx, iyy, izz, ii) =
     1  ( Efy_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_t ( ixx, iyy, n_cc_components + 3)) * 0.5

        ELSE                              ! Ghost face.

          Efy_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_t ( ixx, iyy, n_cc_components + 3)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost". 

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the bottom sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! bottom sibling patch.

! STEP 2.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of 
! "Flux_snd_t ( :, :, :, yr_patch_in_proc)", just as we would do 
! in a Cartesian mesh code. This is saved in "Flux_recev_b" on the 
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) Now that the data is in "Flux_recev_b" we can unpack it into 
! "Flux_z_amr, EntFlux_z_amr, Efx_amr, Efy_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! top sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access top sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for top sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 7)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 7)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 7)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "Flux_snd_t ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_b".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_t ( :, :, :, yr_patch_in_proc)" to "Flux_recev_b".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_ZT_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_b, num_of_reals_z, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,     ! remote data on yr.
     1               MPI_REAL8, Data_Flux_ZT_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_ZT_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) Unpack data from "Flux_recev_b" to
! "Flux_z_amr, EntFlux_z_amr, Efx_amr, Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = 0

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1  ( Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii)
     1  + Flux_recev_b ( ixx, iyy,      1: n_cc_components)) * 0.5

          EntFlux_z_amr ( ixx, iyy, izz, ii) = 
     1  ( EntFlux_z_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_b  ( ixx, iyy,      n_cc_components + 1)) * 0.5

        ELSE                              ! Ghost face.

          Flux_z_amr   ( ixx, iyy, izz, 1: n_cc_components, ii) = 
     1    Flux_recev_b ( ixx, iyy,      1: n_cc_components)

          EntFlux_z_amr ( ixx, iyy, izz, ii) = 
     1    Flux_recev_t  ( ixx, iyy,      n_cc_components + 1)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_b" to "Efx_amr".

      DO iyy = 0 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 1 - NGhost, ntp_zones  + NGhost 

        izz = 0

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones)) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_b ( ixx, iyy, n_cc_components + 2)) * 0.5

        ELSE                              ! Ghost face.

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_b ( ixx, iyy, n_cc_components + 2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost". 

! --------------------
! Unpack data from "Flux_recev_b" to "Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost                           
      DO ixx = 0 - NGhost, ntp_zones  + NGhost 

        izz = 0

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Efy_amr      ( ixx, iyy, izz, ii) =
     1  ( Efy_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_b ( ixx, iyy, n_cc_components + 3)) * 0.5

        ELSE                              ! Ghost face.

          Efy_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_b ( ixx, iyy, n_cc_components + 3)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".        
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost". 

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE UNPACK_DEFECT_FLUXES ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called before the synchronization of fluxes and 
! electric fields at the coinciding faces and edges, respectively.
! For that, the amr fluxes and electric fields at the coinciding faces and edges
! are populated in this subroutine from the "Flux_snd_X_amr" arrays.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        izt, izb,

     1        ixs, iys, izs

! ----------------------------------------------------------------------
! Synchronize on entry.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Assign Defect fluxes and electric fields at the coinciding faces and edges.

      DO ii = 1, npatch_per_proc

! --------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz )                             &

!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components,                                              &
!$OMP&    Defect_Flux_x_amr, Flux_x_amr, Defect_Efx_amr, Efx_amr,       &
!$OMP&    Defect_Flux_y_amr, Flux_y_amr, Defect_Efy_amr, Efy_amr,       &
!$OMP&    Defect_Flux_z_amr, Flux_z_amr, Defect_Efz_amr, Efz_amr,       &
!$OMP&    Defect_EntFlux_x_amr, EntFlux_x_amr,                          &
!$OMP&    Defect_EntFlux_y_amr, EntFlux_y_amr,                          &
!$OMP&    Defect_EntFlux_z_amr, EntFlux_z_amr)                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0,  nr_zones
        DO iyy = 0, ntp_zones
        DO ixx = 0, ntp_zones

! --------------------
! x-face.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_x_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Flux_x_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Flux_x_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_x_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_x_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_x_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "x-face".

! --------------------
! y-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_y_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Flux_y_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Flux_y_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_y_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_y_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_y_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "y-face".

! --------------------
! z-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Flux_z_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Flux_z_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Flux_z_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_z_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_z_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "z-face".

! --------------------
! x-edge.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efx_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efx_amr ( ixx, iyy, izz, ii)
     1       + Efx_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "x-edge".

! --------------------
! y-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efy_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efy_amr ( ixx, iyy, izz, ii)
     1       + Efy_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "y-edge".

! --------------------
! z-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efz_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efz_amr ( ixx, iyy, izz, ii)
     1       + Efz_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "z-edge".

! --------------------

        END DO ! End of "DO ixx = 0, ntp_zones".
        END DO ! End of "DO iyy = 0, ntp_zones".
        END DO ! End of "DO izz = 0,  nr_zones".

!$OMP END PARALLEL DO

! --------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE UNPACK_DEFECT_FLUXES

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
!
! Written By: Sethupathy S.
!
! In this subroutine the data window handles are created for MHD variables
! and the Flux variables in Steps 1 and 2 respectively.

! STEP 1.
! Create windows for message passing for,
!   "Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
!    Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
!    Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle".
! These handles points to,
!   "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr, U_snd_d_amr,"
!   "U_snd_t_amr, U_snd_b_amr".
!
! STEP 2.
! Create windows for message passing for,
!   "Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
!    Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
!    Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle".
! These handles points to,
!   "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr,"
!   "Flux_snd_d_amr, Flux_snd_t_amr, Flux_snd_b_amr".
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: num_of_reals_xy, num_of_reals_z, ierr

! ----------------------------------------------------------------------

! Type of variable being called in bytes is stored in "type_mpi_real".

      CALL MPI_Type_size ( MPI_REAL8, type_mpi_real, ierr)

! ----------------------------------------------------------------------
! STEP 1. Create handles for MHD variables.
! ----------------------------------------------------------------------
! Evaluate number of real variables in the six handles.
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message passing data size is same for "x" and "y" axes
! and it is stored in "num_of_reals_xy". The message passing data size
! along "z" axis is different and it is stored in "num_of_reals_z".

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 ) * Max_patches_per_proc

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 ) * Max_patches_per_proc

! Calculate window size in bytes.

      window_size_xy = type_mpi_real * num_of_reals_xy
      window_size_z  = type_mpi_real * num_of_reals_z

! Create six windows for,
!   "Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
!    Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
!    Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle" to hold values.

! "MPI_INFO_NULL" is something that you just do in MPI.

! The six window handle names will later allow us to get handles on the
! relevant data arrays "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr, 
! U_snd_d_amr, U_snd_t_amr, U_snd_b_amr" respectively on the remote arrays.

! It is crucial to release it with "MPI_WIN_FREE".

! x-data
      CALL MPI_WIN_CREATE ( U_snd_r_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_l_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_XL_Msg_handle, ierr)

! y-data
      CALL MPI_WIN_CREATE ( U_snd_u_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_YU_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_d_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_YD_Msg_handle, ierr)

! z-data
      CALL MPI_WIN_CREATE ( U_snd_t_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_ZT_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_b_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_ZB_Msg_handle, ierr)

! ----------------------------------------------------------------------
! STEP 2. Create handles for Flux variables.
! ----------------------------------------------------------------------
! Evaluate number of real variables in the six handles.
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message passing data size is same for "x" and "y" axes
! and it is stored in "num_of_reals_xy". The message passing data size
! along "z" axis is different and it is stored in "num_of_reals_z".

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 ) * Max_patches_per_proc

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( n_cc_components + 3 ) * Max_patches_per_proc

! Calculate window size in bytes.

      window_size_xy = type_mpi_real * num_of_reals_xy
      window_size_z  = type_mpi_real * num_of_reals_z

! Create six windows for,
!   "Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
!    Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
!    Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle" to hold values.

! "MPI_INFO_NULL" is something that you just do in MPI.

! The six window handle names will later allow us to get handles on the
! relevant data arrays "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr, 
! Flux_snd_d_amr, Flux_snd_t_amr, Flux_snd_b_amr" respectively on the 
! remote arrays.

! It is crucial to release it with "MPI_WIN_FREE".

! x-data
      CALL MPI_WIN_CREATE ( Flux_snd_r_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_l_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_XL_Msg_handle, ierr)

! y-data
      CALL MPI_WIN_CREATE ( Flux_snd_u_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_d_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_YD_Msg_handle, ierr)

! z-data
      CALL MPI_WIN_CREATE ( Flux_snd_t_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_ZT_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_b_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_ZB_Msg_handle, ierr)

! ----------------------------------------------------------------------

      END SUBROUTINE WINCREATE_MESSAGING

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE WINFREE_MESSAGING ()

! ----------------------------------------------------------------------
!
! Written By: Sethupathy S.
!
! Free up all the MHD handles and Flux handles. 
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: ierr

! ----------------------------------------------------------------------

! Must deallocate. Otherwise expect memory leaks.

! MHD handles.
      CALL MPI_WIN_FREE ( Data_U_XR_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_XL_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_YU_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_YD_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_ZT_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_ZB_Msg_handle, ierr)

! Flux handles.
      CALL MPI_WIN_FREE ( Data_Flux_XR_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_XL_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_YU_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_YD_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_ZT_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_ZB_Msg_handle, ierr)

! ----------------------------------------------------------------------

      END SUBROUTINE WINFREE_MESSAGING

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc, 
     1                 time, dtcur, rmin, rmax )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.                                              
!                                                                        
! This test subroutine computes the error in "U_amr, entropygr_amr,
! Bxf_amr, Byf_amr, Bzf_amr". 
! 
! ----------------------------------------------------------------------

      IMPLICIT NONE                                                      

      INTEGER iproblem, iproc

      REAL time, dtcur, rmin, rmax
                                                                         
      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ix_patch, iy_patch, iz_patch, isector,
     1        itri_patch_LD, itri_patch_LU, 
     1        itri_patch_RD, itri_patch_RU
                                                                         
      INTEGER ixx, iyy, izz, i, ii, izb, ize, im1, ip1

      REAL xc, yc, zc

      REAL err_U ( n_cc_components),
     1     err_Bxf, err_Byf, err_Bzf, 
     1     err_Fx, err_Fy, err_Fz,                                
     1     err_EntFx, err_EntFy, err_EntFz,                       
     1     err_Efx, err_Efy, err_Efz,

! Global error from all processors.
     1     gl_err_U ( n_cc_components),
     1     gl_err_Bxf, gl_err_Byf, gl_err_Bzf,
     1     gl_err_Fx, gl_err_Fy, gl_err_Fz, 
     1     gl_err_EntFx, gl_err_EntFy, gl_err_EntFz,   
     1     gl_err_Efx, gl_err_Efy, gl_err_Efz

      INTEGER :: write_data = 0
                                                                         
! ---------------------------------------------------------------------- 
                                                                         
! Synchronize on entry.                                                  
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
                                                                         
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)                
                                                                         
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)                
      me = my_rank + 1                                                   
                                                                         
      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages      
                                                                         
! ---------------------------------------------------------------------- 

      IF ( me == 1 ) THEN

      IF ( iproblem .NE. 0 ) THEN                                 
        WRITE ( 6, *)"Stopping at TEST_MESSAGING_ERROR_CALCULATION"  
        WRITE ( 6, *)"Expecting iproblem == 0"                    
        WRITE ( 6, *)"iproblem = ",iproblem
        STOP                                                             
      END IF

      END IF

! Initialize the error variables to zero.
      err_U ( :) = 0.0;   gl_err_U ( :) = 0.0
     
      err_Bxf = 0.0;      err_Byf = 0.0;      err_Bzf = 0.0
      gl_err_Bxf = 0.0;   gl_err_Byf = 0.0;   gl_err_Bzf = 0.0

      err_Fx = 0.0;      err_Fy = 0.0;     err_Fz = 0.0
      err_EntFx = 0.0;   err_EntFy = 0.0;  err_EntFz = 0.0
      err_Efx = 0.0;     err_Efy = 0.0;    err_Efz = 0.0

      gl_err_Fx = 0.0;      gl_err_Fy = 0.0;     gl_err_Fz = 0.0
      gl_err_EntFx = 0.0;   gl_err_EntFy = 0.0;  gl_err_EntFz = 0.0
      gl_err_Efx = 0.0;     gl_err_Efy = 0.0;    gl_err_Efz = 0.0

! ---------------------------------------------------------------------- 
! ---------------------------------------------------------------------- 

      DO ii = 1, npatch_per_proc

! ---------------------------------------------------------------------- 
! ---------------------------------------------------------------------- 

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ---------------------------------------------------------------------- 
! sector number and patch location in a sector.

      isector  = isector_ilevel_amr ( 1, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! ---------------------------------------------------------------------- 

      vert_coords_patch( :, :, :, :) = vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) = imask_active_amr ( :, :, :, ii)

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! ---------------------------------------------------------------------- 

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".
! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ---------------------------------------------------------------------- 
! Compute error in "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------
! Change the ghost zone inclusion along "z-axis", if the patch is at inner or 
! outer r-boundary of the simulation, where there is no MPI messaging.
      izb = 1 - NGhost
      ize = nr_zones + NGhost

      IF ( iz_patch == 1         ) izb = 1         ! patch at inner r-boundary.
      IF ( iz_patch == nr_patches) ize = nr_zones  ! patch at outer r-boundary.

! ---------------------------------------------------------------------- 

!      IF ( me == iproc ) write_data = 1

      U ( :, :, :, :, 1) = U_amr ( :, :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)

      Flux_x ( :, :, :, :) = Flux_x_amr ( :, :, :, :, ii)
      Flux_y ( :, :, :, :) = Flux_y_amr ( :, :, :, :, ii)
      Flux_z ( :, :, :, :) = Flux_z_amr ( :, :, :, :, ii)

      EntFlux_x ( :, :, :) = EntFlux_x_amr ( :, :, :, ii)
      EntFlux_y ( :, :, :) = EntFlux_y_amr ( :, :, :, ii)
      EntFlux_z ( :, :, :) = EntFlux_z_amr ( :, :, :, ii)

      Efx ( :, :, :) = Efx_amr ( :, :, :, ii)
      Efy ( :, :, :) = Efy_amr ( :, :, :, ii)
      Efz ( :, :, :) = Efz_amr ( :, :, :, ii)

      CALL TEST_MPI_ERROR_IN_B_XYZ_FACES
                                                                         
     1         ( nr_zones, ntp_zones, NGhost,                            
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,     
     1           n_fc_modes, n_cc_bfld_modes,                            
                                                                         
     1           iproblem, time, dtcur,                                  
     1           rmin, rmax,                                             
                                                                         
     1           vert_coords_patch, imask_active_patch,                  
     1           v_scale_length,                                         
                                                                         
     1           xf_area_patch, yf_area_patch, zf_area_patch,            
     1           v_quadpts_patch, v_quadwt_patch,                        
                                                                         
     1           U, V, entropygr, Bxf, Byf, Bzf,                         

     1           Flux_x, Flux_y, Flux_z, 
     1           EntFlux_x, EntFlux_y, EntFlux_z, Efx, Efy, Efz,
                                                                         
     1           gamma_eos_gr,                                           
                                                                         
     1           izb, ize, write_data,

     1           err_U, err_Bxf, err_Byf, err_Bzf,
     1           err_Fx, err_Fy, err_Fz,                                
     1           err_EntFx, err_EntFy, err_EntFz,                       
     1           err_Efx, err_Efy, err_Efz)   

! ----------------------------------------------------------------------
! ---------------------------------------------------------------------- 

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ---------------------------------------------------------------------- 

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

! ----------------------------------------
! Reduce and get max error for "err_U" from all processors. 
! Then broadcast the global maximum "gl_err_U" to all processors from 
! my_rank == 0  processor.

      CALL MPI_REDUCE ( err_U, gl_err_U, n_cc_components, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
      CALL MPI_Bcast   ( gl_err_U, n_cc_components, MPI_REAL8, 0,
     1                   MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_Bxf, err_Byf, err_Bzf" from all processors. 
! Then broadcast the global maximum "gl_err_Bxf, gl_err_Byf, gl_err_Bzf" 
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Bxf, gl_err_Bxf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Byf, gl_err_Byf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Bzf, gl_err_Bzf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

      CALL MPI_Bcast ( gl_err_Bxf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Byf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Bzf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

! ----------------------------------------
! Reduce and get max error for "err_Fx, err_Fy, err_Fz" from all processors. 
! Then broadcast the global maximum "gl_err_Fx, gl_err_Fy, gl_err_Fz" 
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Fx, gl_err_Fx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fy, gl_err_Fy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fz, gl_err_Fz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

      CALL MPI_Bcast ( gl_err_Fx, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fy, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fz, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

! ----------------------------------------
! Reduce and get max error for "err_EntFx, err_EntFy, err_EntFz" from all processors. 
! Then broadcast the global maximum "gl_err_EntFx, gl_err_EntFy, gl_err_EntFz" 
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_EntFx, gl_err_EntFx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_EntFy, gl_err_EntFy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_EntFz, gl_err_EntFz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

      CALL MPI_Bcast ( gl_err_EntFx, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_EntFy, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_EntFz, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

! ----------------------------------------
! Reduce and get max error for "err_Efx, err_Efy, err_Efz" from all processors. 
! Then broadcast the global maximum "gl_err_Efx, gl_err_Efy, gl_err_Efz" 
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Efx, gl_err_Efx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Efy, gl_err_Efy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Efz, gl_err_Efz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

      CALL MPI_Bcast ( gl_err_Efx, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Efy, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Efz, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           

! ---------------------------------------------------------------------- 

      IF ( me == iproc ) THEN

      WRITE ( 6, *)"err_U ( 2) = ",err_U ( 2)
      WRITE ( 6, *)"err_U ( 3) = ",err_U ( 3)
      WRITE ( 6, *)"err_U ( 4) = ",err_U ( 4)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_U ( 2) = ",gl_err_U ( 2)
      WRITE ( 6, *)"gl_err_U ( 3) = ",gl_err_U ( 3)
      WRITE ( 6, *)"gl_err_U ( 4) = ",gl_err_U ( 4)
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Bxf = ",err_Bxf
      WRITE ( 6, *)"err_Byf = ",err_Byf
      WRITE ( 6, *)"err_Bzf = ",err_Bzf
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Bxf = ",gl_err_Bxf
      WRITE ( 6, *)"gl_err_Byf = ",gl_err_Byf
      WRITE ( 6, *)"gl_err_Bzf = ",gl_err_Bzf
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Fx = ",err_Fx
      WRITE ( 6, *)"err_Fy = ",err_Fy
      WRITE ( 6, *)"err_Fz = ",err_Fz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Fx = ",gl_err_Fx
      WRITE ( 6, *)"gl_err_Fy = ",gl_err_Fy
      WRITE ( 6, *)"gl_err_Fz = ",gl_err_Fz
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_EntFx = ",err_EntFx
      WRITE ( 6, *)"err_EntFy = ",err_EntFy
      WRITE ( 6, *)"err_EntFz = ",err_EntFz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_EntFx = ",gl_err_EntFx
      WRITE ( 6, *)"gl_err_EntFy = ",gl_err_EntFy
      WRITE ( 6, *)"gl_err_EntFz = ",gl_err_EntFz
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Efx = ",err_Efx
      WRITE ( 6, *)"err_Efy = ",err_Efy
      WRITE ( 6, *)"err_Efz = ",err_Efz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Efx = ",gl_err_Efx
      WRITE ( 6, *)"gl_err_Efy = ",gl_err_Efy
      WRITE ( 6, *)"gl_err_Efz = ",gl_err_Efz
      WRITE ( 6, *)"------------------------------------------------"

      END IF

! ----------------------------------------------------------------------
                                                                         
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
! Synchronize before exit.                                               
                                                                         
! ----------------------------------------------------------------------

      END SUBROUTINE TEST_MESSAGING_ERROR_CALCULATION

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE TEST_XY_MESSAGING

! ---------------------------------------------------------------------- 
! Written By: Sethupathy S.  
!
! This test subroutine so far prints-out the Fluid variables
! assuming 1 patch/sector. It can be modified to adopt for different cases
! as well.
!
! STEPS:
! ------
! - Initialize the active zones and faces in "U_amr, entropygr_amr, 
!   Bxf_amr, Byf_amr, Bzf_amr".
! - Prints out the values before messaging
! - Do messaging
! - Prints out the values after messaging
!
! ---------------------------------------------------------------------- 

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc

      INTEGER ixx, iyy, i, ii

! ---------------------------------------------------------------------- 
                                                                         
! Synchronize on entry.                                                  
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
                                                                         
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)                
                                                                         
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)                
      me = my_rank + 1                                                   
                                                                         
      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages      
                                                                         
! ---------------------------------------------------------------------- 
! Fill in the values in U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, :, :) = 0.0      
      entropygr_amr ( :, :, :, :) = 0.0

      Bxf_amr ( :, :, :, :) = 0.0
      Byf_amr ( :, :, :, :) = 0.0
      Bzf_amr ( :, :, :, :) = 0.0

      DO iyy = 0, ntp_zones
        DO ixx = 0, ntp_zones

          IF( ixx > 0 .AND. iyy > 0 ) THEN
            U_amr ( ixx, iyy, :,   :, :) = iyy + 1
            entropygr_amr ( ixx, iyy, :,  :) = iyy + 1

            Bzf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

          IF( iyy > 0 ) THEN
            Bxf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

          IF( ixx > 0 ) THEN
            Byf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

        END DO
      END DO

! ----------------------------------------------------------------------

      i = 1
      CALL PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( i)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
! Enable/disable the function calls as needed for testing

      CALL LOAD_LOCAL_BUFFERS_UXY ()

      CALL FETCH_BUFFERS_FILL_IN_UX ()

!      CALL LOAD_LOCAL_BUFFERS_UXY ()

!      CALL FETCH_BUFFERS_FILL_IN_UY ()

!      CALL LOAD_LOCAL_BUFFERS_UXY ()
 
!      CALL FETCH_BUFFERS_FILL_IN_UY ()

! ----------------------------------------------------------------------

      CALL WINFREE_MESSAGING ()

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( i)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
                                                                         
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
! Synchronize before exit.                                               
                                                                         
! ----------------------------------------------------------------------

      END SUBROUTINE TEST_XY_MESSAGING

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( isector )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.

! Prints the values of zones assuming 1 patch per sector.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER isector

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ipatch, ii, i, j, ixx, iyy, izz, icomp 

       REAL U_1D ( 1 - NGhost: ntp_zones + NGhost),
     1      Bxf_1D ( 0 - NGhost: ntp_zones + NGhost),
     1      Byf_1D ( 0 - NGhost: ntp_zones + NGhost),
     1      Bzf_1D ( 0 - NGhost: ntp_zones + NGhost)


! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------

      IF ( me == isector) THEN

! ----------------------------------------------------------------------

      IF ( npatch_per_proc .NE. 1 ) THEN
        WRITE ( 6, *)"Stopping at PRINT_VARIABLES_ONE_PATCH_PER_SECTOR"
        WRITE ( 6, *)"Expecting npatch_per_proc == 1"
        WRITE ( 6, *)"npatch_per_proc = ",npatch_per_proc
        STOP
      END IF

      IF ( nimages .NE. 6 ) THEN
        WRITE ( 6, *)"Stopping at PRINT_VARIABLES_ONE_PATCH_PER_SECTOR"
        WRITE ( 6, *)"Expecting nimages == 6"
        WRITE ( 6, *)"nimages = ",nimages
        STOP
      END IF

      ii = 1

! ----------------------------------------------------------------------

      IF ( isector == iactive_FromProc ( ii, 1)) THEN

! ------------------------------------------------------------

      izz = 1
      icomp = 1

! ----------

      WRITE ( 6, *)"U_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        U_1D ( :) = U_amr ( 1 - NGhost: ntp_zones + NGhost,        
     1                      iyy, izz, icomp, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost, 
     1             U_1D ( 1 - NGhost: ntp_zones + NGhost), "I")

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Bxf_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        Bxf_1D  ( 0 - NGhost: ntp_zones + NGhost) = 
     1  Bxf_amr ( 0 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 0 - NGhost, ntp_zones + NGhost, Bxf_1D, "R" )

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Byf_amr"

      DO iyy = ntp_zones + NGhost, 0 - NGhost, -1

        Byf_1D  ( 1 - NGhost: ntp_zones + NGhost) = 
     1  Byf_amr ( 1 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost, 
     1           Byf_1D ( 1 - NGhost: ntp_zones + NGhost), "R" )

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Bzf_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        Bzf_1D  ( 1 - NGhost: ntp_zones + NGhost) = 
     1  Bzf_amr ( 1 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost,
     1           Bzf_1D ( 1 - NGhost: ntp_zones + NGhost), "I")

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ------------------------------------------------------------

      END IF ! End of "IF ( isector == iactive_FromProc ( ii, 1)) THEN".

! ----------------------------------------------------------------------

      END IF ! End of "IF ( me == isector ) THEN".

! ----------------------------------------------------------------------
                                                                         
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
! Synchronize before exit.                                               
                                                                         
! ----------------------------------------------------------------------

      END SUBROUTINE PRINT_VARIABLES_ONE_PATCH_PER_SECTOR

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( ACCURACY_ANALYSIS == 1)

      SUBROUTINE MAKE_GLOBAL_ACCURACY_CS ( iproblem, 
     1                                     rmin, rmax, time, dtcur )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.
!
! Compute the accuracy of the simulation. It prints out the error after
! each time step.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER iproblem

      REAL rmin, rmax, time, dtcur

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, isector, ix_patch, iy_patch, iz_patch

      INTEGER nRK_stage, irk, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum, ibdy_inner, ibdy_outer

      INTEGER, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,               
     1                     1 - NGhost: ntp_zones + NGhost, 2) :: bcarr_r

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Now we loop over all "npatch_per_proc" patches on each processor and
! build fluxes and electric fields.

      V_min_reduce ( :) = 1.0e32
      V_max_reduce ( :) = -1.0e32

      U_L1_reduce      ( :) = 0.0
      U_Linf_reduce    ( :) = 0.0
      U_zc_Linf_reduce ( :) = 0.0

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! STEP 1) First reverse all that we wrote into the big "*_amr" arrays during
! initialization.

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! ------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! ------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! ------------------------------

      U ( :, :, :, 1:n_cc_components, 1) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr ( :, :, :, 1) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)      

      flattengr ( :, :, :) = flatten_amr ( :, :, :, ii)

! ----------------------------------------------------------------------

!$acc update device ( 

!$acc&       time, dtcur, cfl_coef,                                     &

!$acc&       isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch,     &
!$acc&       vert_coords_patch, imask_active_patch,                     &

!$acc&       itri_patch_LD, itri_patch_RD, itri_patch_RU,               &
!$acc&       itri_patch_LU, istencil_patch,                             &

!$acc&       U, entropygr, Bxf, Byf, Bzf, flattengr )


! ----------------------------------------------------------------------
! STEP 2) Now build all the rest of the geometry information for this patch.

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! --------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes, n_cc_modes_one_higher, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! STEP 3) With geometry built, now do limiting, predictor step, flux
! and field evaluation.

! --------------------------------------------------

! Pad the boundaries in the radial direction if this is the innermost
! radial patch or the outermost radial patch.

! ------------------------------

      IF ( iz_patch .EQ. 1) THEN  

! --------------------
! Inner radial boundary needs to the padded. 

! Most of the geomesh simulations require the re-initialization of boundary
! conditions at the inner boundary. The ADD_BOUNDARY subroutine does that.

      ibdy_inner = 1   ! 1 == Impose Add_boundary at inner; 0 == No Add_Bdy
      ibdy_outer = 1   ! 1 == Impose Add_boundary at outer; 0 == No Add_Bdy

!$acc update device ( ibdy_inner, ibdy_outer)

      CALL ADD_BOUNDARY

     1         ( ibdy_inner, ibdy_outer,

     1           nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_centroid_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           gamma_eos_gr,

     1           V_min_local, V_max_local)

! --------------------

      END IF

! ------------------------------

      IF ( iz_patch .EQ. nr_patches) THEN
! Outer radial boundary needs to the padded.

! Most of the geomesh simulations require to implement the outflow  boundary
! conditions at the outer boundary. The PAD_BOUNDARY subroutine does that.

      bcarr_r ( :, :, 1)  = 17   ! No inner boundary enforcement here.
      bcarr_r ( :, :, 2)  = 17    ! Enforce outflow conditions at Outer boundary.

!$acc update device ( bcarr_r)

      CALL PAD_BOUNDARY_R

     1    ( nr_zones, ntp_zones, NGhost,
     1      n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1      n_fc_modes, n_cc_bfld_modes,

     1      vert_coords_patch, imask_active_patch,

     1      xf_area_patch, yf_area_patch, zf_area_patch,
     1      v_volume_patch,

     1      v_centroid_patch,

     1      U, entropygr, Bxf, Byf, Bzf,

     1      bcarr_r)

      END IF

! --------------------------------------------------

! Apply limiter to MHD variables.

      CALL LIMIT_ALL_SLOPES_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           iz_patch, nr_patches,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! Store local min and max for each patch that resides on this processor.
      V_min_amr ( :, ii) = V_min_local ( :)
      V_max_amr ( :, ii) = V_max_local ( :)

! For all the patches that reside on this processor, make a reduction of their
! min and max values.
      V_min_reduce = AMIN1 ( V_min_reduce, V_min_local)
      V_max_reduce = AMAX1 ( V_max_reduce, V_max_local)

! --------------------------------------------------
! Compute the error.

      CALL ACCURACY_CALCULATION_PATCH_CS

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           ii_mode, jj_mode, kk_mode,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_quadpts_patch, v_quadwt_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           Flux_x, Flux_y, Flux_z,
     1           EntFlux_x, EntFlux_y, EntFlux_z, Efx, Efy, Efz,
     1           gamma_eos_gr,

     1           U_L1, U_Linf, U_zc_Linf)

      U_L1_reduce ( :) = U_L1_reduce ( :) + U_L1 ( :)

      U_Linf_reduce ( :) = MAX ( U_Linf_reduce ( :), U_Linf ( :))

      U_zc_Linf_reduce ( :) = MAX ( U_zc_Linf_reduce( :), U_zc_Linf( :))

! --------------------
#if 0
      IF ( me == 1) THEN

        WRITE ( 6, *)"U_L1       = ", U_L1
        WRITE ( 6, *)"U_Linf     = ", U_Linf
        WRITE ( 6, *)"U_zc_Linf  = ", U_zc_Linf

      END IF
#endif
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      U_L1_reduce ( :) = U_L1_reduce ( :) / REAL ( npatch_per_proc)

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

! Now make a reduction across all processors for "V_min_reduce, V_max_reduce".

      U_L1      ( :) = U_L1_reduce      ( :)
      U_Linf    ( :) = U_Linf_reduce    ( :)
      U_zc_Linf ( :) = U_zc_Linf_reduce ( :)


      CALL MPI_ALLREDUCE ( U_L1,              ! Input for reudction
     1                     U_L1_reduce,       ! Global sum
     1                     n_cc_components,   ! length of reduction variable   
     1                     MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      U_L1_reduce ( :) = U_L1_reduce ( :) / REAL ( nimages)


      CALL MPI_ALLREDUCE ( U_Linf,            ! Input for reudction
     1                     U_Linf_reduce,     ! Global maximum
     1                     n_cc_components,   ! length of reduction variable   
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE ( U_zc_Linf,         ! Input for reudction
     1                     U_zc_Linf_reduce,  ! Global maximum
     1                     n_cc_components,   ! length of reduction variable   
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)


! Write out the global min and max of all variables from processor 0.

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" ---------------------- Error --------------- "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_L1      = ", U_L1_reduce ( :)
      WRITE ( 6, *)"U_Linf    = ", U_Linf_reduce ( :)
      WRITE ( 6, *)"U_zc_Linf = ", U_zc_Linf_reduce ( :)
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE MAKE_GLOBAL_ACCURACY_CS

#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE WRITE_SILO_ONE_BLOCK ( filenum, time, dtcur )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.

! Tests the silo mesh write out for one block.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER filenum
      
      REAL time, dtcur

! ----------------------------------------------------------------------

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ii

      INTEGER block_index              ! global number for the patch.
       
      INTEGER itri_patch ( 4)          ! extracted from itri_patch_amr

! ---------------------------------------------------------------------- 
                                                                         
! Synchronize on entry.                                                  
!     CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)                           
                                                                         
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)                
                                                                         
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)                
      me = my_rank + 1                                                   
                                                                         
      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages      
                                                                         
! ---------------------------------------------------------------------- 

      DO ii = 1, npatch_per_proc

! ---------------------------------------------------------------------- 

! Global patch location, counted first along the processor direction,
      block_index = me + ( ii - 1) * nimages

      itri_patch ( :) = itri_patch_amr ( :, ii)

      vert_coords_patch( :, :, :, :) = vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) = imask_active_amr ( :, :, :, ii)

! Call the tester subroutine for only one patch.
      IF ( ii == 1 ) THEN
#if 0
      CALL WRITE_ONE_BLOCK ( block_index, 
     1                       nr_zones,  NGhost,
     1                       ntp_zones, NGhost,
     1                       itri_patch,
     1                       vert_coords_patch )
#endif
      END IF

! ----------------------------------------------------------------------

      END DO     ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      END SUBROUTINE WRITE_SILO_ONE_BLOCK

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END MODULE RIEMANN_3D_CS_AMR_SUPPORT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FILL_FLUX_SND_PATCH 

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes,

     1           Flux_x, Flux_y, Flux_z,
     1           EntFlux_x, EntFlux_y, EntFlux_z,
     1           Efx, Efy, Efz,

     1           Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,
     1           Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch)

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after building the fluxes and electric fields,

! In this subroutine, only the boundary data that is stored in the 
! "Flux_xyz, Ent_Flux_xyz, Efxyz" arrays are loaded into the 
! "Flux_snd_r_patch", "Flux_snd_l_patch", "Flux_snd_u_patch",
! "Flux_snd_d_patch", "Flux_snd_t_patch", "Flux_snd_b_patch" arrays. 

! After that the "Flux_snd_X_patch" arrays will be moved to the CPU and
! then they will fill their respective "_amr" arrays on the CPU.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1        n_fc_modes, n_cc_bfld_modes

! ----------

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,                  
     1                  0 - NGhost: ntp_zones + NGhost, 
     1                  0 - NGhost:  nr_zones + NGhost, 
     1                  n_cc_components) :: 
     1      Flux_x, Flux_y, Flux_z
                                                                         
      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,                  
     1                  0 - NGhost: ntp_zones + NGhost, 
     1                  0 - NGhost:  nr_zones + NGhost)::  
     1      EntFlux_x, EntFlux_y, EntFlux_z                              
                                                                         
      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,                  
     1                  0 - NGhost: ntp_zones + NGhost, 
     1                  0 - NGhost:  nr_zones + NGhost):: Efx, Efy, Efz

! ----------

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components + 3) ::
     1   Flux_snd_r_patch, Flux_snd_l_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components + 3) ::
     1   Flux_snd_u_patch, Flux_snd_d_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  n_cc_components + 3) ::
     1   Flux_snd_t_patch, Flux_snd_b_patch

! ----------------------------------------------------------------------
! Variables that are declared below are internal to this subroutine.

      INTEGER ixx, iyy, izz, icomp,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        izt, izb,

     1        ixs, iys, izs

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

!$acc kernels
      Flux_snd_r_patch ( :, :, :) = 0.0
      Flux_snd_l_patch ( :, :, :) = 0.0
      Flux_snd_u_patch ( :, :, :) = 0.0
      Flux_snd_d_patch ( :, :, :) = 0.0
      Flux_snd_t_patch ( :, :, :) = 0.0
      Flux_snd_b_patch ( :, :, :) = 0.0
!$acc end kernels

! ----------------------------------------------------------------------
! First store right and left x-data in "Flux_snd_r_patch, Flux_snd_l_patch".
! Then store up and down y-data in "Flux_snd_u_patch, Flux_snd_d_patch".

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_patch, Flux_snd_l_patch". First for the face centered fluxes,
! "Flux_x, EntFlux_x" then for the edge centered variables,
! "Efy, Efz".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_x, EntFlux_x" data in
! "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch,          &
!$OMP&    Flux_x, EntFlux_x )                                           &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, Flux_x, EntFlux_x)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice
          ixr = ntp_zones      ! location of right slice

! ----------

          Flux_snd_r_patch (      iyy, izz, 1: n_cc_components) =
     1              Flux_x ( ixr, iyy, izz, 1: n_cc_components) 

          Flux_snd_l_patch (      iyy, izz, 1: n_cc_components) =
     1              Flux_x ( ixl, iyy, izz, 1: n_cc_components) 

! ----------

          Flux_snd_r_patch (      iyy, izz, n_cc_components + 1) =
     1           EntFlux_x ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, n_cc_components + 1) =
     1           EntFlux_x ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy" data in "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch, Efy )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, Efy )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl )

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Flux_snd_r_patch (      iyy, izz, n_cc_components + 2) =
     1                 Efy ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, n_cc_components + 2) =
     1                 Efy ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz" data in "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch, Efz )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, Efz )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.          

          Flux_snd_r_patch (      iyy, izz, n_cc_components + 3) =
     1                 Efz ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, n_cc_components + 3) =
     1                 Efz ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "Flux_snd_u_patch, Flux_snd_d_patch". First for the face centered fluxes,
! "Flux_y, EntFlux_y" then for the edge centered variables,
! "Efz, Efx".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_y, EntFlux_y" data in
! "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch, Flux_y,  &
!$OMP&    EntFlux_y )                                                   &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, Flux_y, EntFlux_y)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice
          iyu = ntp_zones      ! location of up   slice

! ----------

          Flux_snd_u_patch ( ixx,      izz, 1: n_cc_components) =
     1              Flux_y ( ixx, iyu, izz, 1: n_cc_components) 

          Flux_snd_d_patch ( ixx,      izz, 1: n_cc_components) =
     1              Flux_y ( ixx, iyd, izz, 1: n_cc_components) 

! ----------

          Flux_snd_u_patch ( ixx,      izz, n_cc_components + 1) =
     1           EntFlux_y ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, n_cc_components + 1) =
     1           EntFlux_y ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz" data in "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch, Efz )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, Efz )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Flux_snd_u_patch ( ixx,      izz, n_cc_components + 2) =
     1                 Efz ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, n_cc_components + 2) =
     1                 Efz ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx" data in "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch, Efx )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, Efx )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.          

          Flux_snd_u_patch ( ixx,      izz, n_cc_components + 3) =
     1                 Efx ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, n_cc_components + 3) =
     1                 Efx ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Store top and bottom z--data in "Flux_snd_t_patch, Flux_snd_b_patch".

! ----------------------------------------------------------------------
! Loading "Flux_snd_t_patch, Flux_snd_b_patch". First for the face centered fluxes,
! "Flux_z, EntFlux_z" then for the edge centered variables,
! "Efx, Efy".

! ------------------------------------------------------------
! Storing face centered fluxes "Flux_z, EntFlux_z" data in
! "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch, Flux_z,  &
!$OMP&    EntFlux_z )                                                   &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, Flux_z, EntFlux_z)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

! ----------

          Flux_snd_t_patch ( ixx, iyy,      1: n_cc_components) =
     1              Flux_z ( ixx, iyy, izt, 1: n_cc_components)

          Flux_snd_b_patch ( ixx, iyy,      1: n_cc_components) =
     1              Flux_z ( ixx, iyy, izb, 1: n_cc_components)

! ----------

          Flux_snd_t_patch ( ixx, iyy,      n_cc_components + 1) =
     1           EntFlux_z ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      n_cc_components + 1) =
     1           EntFlux_z ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx" data in "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch, Efx )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, Efx )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_patch ( ixx, iyy,      n_cc_components + 2) =
     1                 Efx ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      n_cc_components + 2) =
     1                 Efx ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy" data in "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch, Efy )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, Efy )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_patch ( ixx, iyy,      n_cc_components + 3) =
     1                 Efy ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      n_cc_components + 3) =
     1                 Efy ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE FILL_FLUX_SND_PATCH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if (1 == 1)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM RIEMANN_GEOMESH_CS_MHD

! ----------------------------------------------------------------------
! Written by : Dinshaw Balsara and Sethupathy S.

      USE RIEMANN_3D_CS_AMR_SUPPORT

      IMPLICIT NONE

#include "riemann.com"

#include "cool.com"

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, isector, ix_patch, iy_patch, iz_patch

      INTEGER nRK_stage, irk, iproblem, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum, ibdy_inner, ibdy_outer,
     1        nr_geom_zones

      REAL rmin, rmax, rratio, time, dtcur, dtcurp1, 
     1     a, b, dt_defcor, time_defcor,
     1     timeimage, dtimage, timestop,
     1     timing_write1, timing_write2,
     1     timing_timestep1, timing_timestep2,

     1     geom_dr1

!$acc declare create ( time, dtcur, dtcurp1, ibdy_inner, ibdy_outer,    &
!$acc&                 dt_defcor, time_defcor )

      REAL sector_vertices ( 4, 3, 6)

! ----------------------------------------------------------------------
! File names used for imaging.

      CHARACTER*13 silofil

      CHARACTER*8 rhofil, vlxfil, vlyfil, vlzfil, prsfil,
     1            bfxfil, bfyfil, bfzfil, dumpfil,
     1            speciefil

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)::
     1   r4gr

! Enforcing radial boundary conditions. Used for PAD_BOUNDARY subroutine.

      INTEGER, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost, 2) :: bcarr_r

! ----------------------------------------------------------------------

      REAL timing_info1, timing_info2, alph_cak

! ----------------------------------------------------------------------
! Variables used in testing the MPI messaging.

      INTEGER iproc, nsilo_images

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "riemann_geomesh.out")

! ----------------------------------------------------------------------
! Initialize MPI here.

      CALL MPI_INIT (ierr)

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Fix up the total number of processors "nimages" and my number in that
! space of processors, "me".

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

! Create windows for message passing. Establishes handles for MPI-3
! one-sided get calls.
      CALL WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

! Also create fluid variable arrays that are needed by the algorithm and only need
! to reside on the GPU. These arrays are not directly used for the update
! and will never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(                                                      &
!$acc&           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch, &
!$acc&           xf_area_patch, xf_centroid_patch,                      &

!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &

!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &

!$acc&           v_scale_length,                                        &

!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &

!$acc&           isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch, &

!$acc&           Coef_Basis_Moments,                                    &
!$acc&           itri_patch_LD, itri_patch_RD, itri_patch_RU,           &
!$acc&           itri_patch_LU, istencil_patch,                         &

!$acc&         U, entropygr, Bxf, Byf, Bzf,                             &
!$acc&         V, B1_zc, B2_zc, B3_zc, Source,                          &
!$acc&         U_skinny, entropygr_skinny,                              &
!$acc&         Bxf_skinny, Byf_skinny, Bzf_skinny,                      &

!$acc&         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,       &
!$acc&         divvelgr, msonicgr, flattengr, workgr,                   &
!$acc&         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                   &

!$acc&         V_min_local, V_max_local,                                &

!$acc&         Flux_x, Flux_y, Flux_z,                                  &
!$acc&         EntFlux_x, EntFlux_y, EntFlux_z,                         &
!$acc&         Efx, Efy, Efz,                                           &
!$acc&         bcarr_r,                                                 &
!$acc&         dU_dt_predgr, dEnt_dt_predgr, U_TAvg_modes_gr,           &
!$acc&         flag_entprs, flag_Tfloor,                                &

!$acc&         Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,    &
!$acc&         Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch,    &

!$acc&         iproblem, rmin, rmax )

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Initializing prs and rho floors, inner and outer radius and fluid variables.

      ilevel = 1  ! By default, we work on 1 level.

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

! Initialize "rmin, rmax" to reasonable values.
      rmin = 2.0
      rmax = 3.5

! If the ratio = 1.0 then there is no ratioing of the radial mesh, when
! the geometric mesh is choosen.
      radial_mesh_ratio = 1.0

      base_grid_rmin = rmin
      base_grid_rmax = rmax

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4           ! Polytropic index for the gas.
      gasconst_euler = 1.0        ! Set up scaled gas constant.
      molewt_euler = 1.4          ! Set up mole wt. of the gas.

      isospeed = 1.0              ! Isothermal sound speed; only for isomhd problems.
      barotropic_stiffness = 0.0  ! Stiffens the isothermal EOS if non-zero.

      tempfloor = 1.0e-5          ! used only when the cooling fn is called.

! Note: The following numbers to be changed for suitable problems.
! If "Mstar" is non-zero then gravity will be switched on in the ADER subroutine.
! If "omega_passive" is non-zero then non-inertial terms will be switched
! on in the ADER subroutine.

      newtonsc = 6.674E-08        ! Newton's constant; in cgs unit
      Mstar = 0.0                 ! Mass of Sun; in grams  =  1.989E33
      omega_passive ( :) = 0.0    ! Rotation axis of star.
      omega_passive ( 3) = 0.0

      omega_active ( :) = 0.0

! ----------
! For multipole simulations.

      iflag_npoles ( 1:4) = 0
      iflag_npoles ( 1)   = 1     ! Enable only dipole by default.

      multipole_Bmag ( 1:4) = 0.0
      multipole_Bmag ( 1)   = 0.0

! ----------
! Define the dipole and qudrapole magnetic field strengths.

      a = pi / 6.0

! monopole.
      iflag_monopole  = 0    ! 0 - no monopole field, 1 - monpole, 2 - split-monopole.
      monopole_u0     = 0.0  ! radial velocity, used in parker spiral solution
      monopole_B ( :) = 0.0

! dipole and quadrupole.
      dipole_B ( 1) = 0.0 * SIN ( a)
      dipole_B ( 2) = 0.0
      dipole_B ( 3) = 0.0 * COS ( a)

      qpole_B ( 1:3, 1:3) = 0.0

! ----------
! Initial time is always zero.
      time = 0.0  

! ----------
! Decide which physical problem you want to initialize with "iproblem". 
! iproblem = 0  for testing the MPI messaging or Regular testing
! iproblem = 2  MHD Blast;
! iproblem = 7  Parker wind test;
! iproblem = 9  SOD test;                iproblem = 10 Brio-Wu test 
! iproblem = 13 VF slow solar wind;

      iproblem = 7

!$acc update device( iproblem, rmin, rmax )

! ----------------------------------------------------------------------

! These parameters control the numerics.

#if ( NUMORDER == 2)
      nRK_stage = 2  ! Two stage RK at second order.
#elif ( NUMORDER == 3)
      nRK_stage = 3  ! Three stage RK at third order.
#elif ( NUMORDER == 4)
      nRK_stage = 5  ! Four stage RK at fourth order.
#endif

! Courant number
      cfl_coef = 0.2

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1   ! Set fractional change needed for entropy fix.(deprecated)

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

! Using 100% multid flux can sometimes mean we take flux contributions from
! a stencil that does not depend on the zones that contribute to that face.
! For that reason, we blend multid RS with 1d RS associated with that face.
! Blending in more multidimensionality also seems to help out with pressure
! positivity without any artificial fixes.
! This value should range between [0.4, 0.6] for good performance.

! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
! flux from each edge.
! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
! The electric field is always fully evaluated by the multid. RS.

      blend_multid_flux = 0.1666

! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
! (almost or fully) stationary contacts. To prevent that we blend in a small
! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]

! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
! reverse in smooth flow.

      blend_mdhll_to_mdhllc = 0.0

! ----------------------------------------------------------------------

! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 5000000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc. (not used)
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Time at which you want the simulation to stop.
      timestop = 1.0

! Choose the number of silo images that are expected to be written-out. 
      nsilo_images = 1

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Populate iproblem specific values of magnetic field, angular velocity,
! sound speed etc., The subroutine which is called here, is present in 
! the applications file, specific to the application.
! This part is, switched OFF in the main version of the code and
! Switched ON when running a specific application; if needed.

#if( CHOOSE_APPL_CAK == 1)

! For CAK winds problem:-

      CALL GET_CAK_PARAMETERS(
     1     iproblem, nsilo_images,
     1     rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive, isospeed, dipole_B, alph_cak,
     1     Mstar, newtonsc, tempfloor, iflag_npoles, multipole_Bmag)

! Fill in the cooling table. different cooling functions can be read in by making
! small modifications to this subroutine.
      CALL FILL_CAK_COOLTABLE ( )

!$acc update device  (/COOL_COM/)

#endif

! ----------------------------------------------------------------------

! update "rmin, rmax".
      base_grid_rmin = rmin
      base_grid_rmax = rmax

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.

      dtimage = timestop / REAL ( nsilo_images)

!$acc update device( iproblem, rmin, rmax )

! ----------------------------------------------------------------------

! Set default multispecies values, just to make sure they are not NANS.

      gamma_specie ( :) = gamma
      molewt_specie ( :) = molewt
      heat_of_formation ( :) = 0.0
      cp_specie ( :) = gasconst * gamma_specie ( :)
     1            / ( ( gamma_specie ( :) - 1.0) * molewt_specie ( :) )
      cv_specie ( :) = cp_specie ( :) / gamma_specie ( :)

! ----------

      gamma_specie_euler ( :) = gamma_euler
      molewt_specie_euler ( :) = molewt_euler
      heat_of_formation_euler ( :) = 0.0
      cp_specie_euler ( :) = gasconst_euler * gamma_specie_euler ( :)
     1 / ( ( gamma_specie_euler ( :) - 1.0) * molewt_specie_euler ( :))
      cv_specie_euler ( :) = cp_specie_euler ( :)
     1                     / gamma_specie_euler ( :)

! ----------
! Make a default initialization of big fluid arrays. This is still done
! on CPU.

      gamma_eos_gr ( :, :, :) = gamma
      gamma_soundspeed_gr ( :, :, :) = gamma
      mean_mol_wt_gr ( :, :, :) = molewt

! ----------------------------------------------------------------------

      divgr_amr ( :, :, :, :) = 0.0

      flag_entprs_amr ( :, :, :, :) = 0.0
      flag_Tfloor_amr  ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! By this time we have updated a lot of variables on the CPU. Let us
! migrate that data to the GPU. We still want to continue working
! on CPU, so we set "do_on_gpu = .FALSE.".

      do_on_gpu = .FALSE.

!$acc update device  (/RIEMANN_COM/)

! ----------------------------------------------------------------------

! Set the variable names used in the silo file.

      silo_allvars ( 1) = "Density "
      silo_allvars ( 2) = "XMomentum "
      silo_allvars ( 3) = "YMomentum "
      silo_allvars ( 4) = "ZMomentum "
      silo_allvars ( 5) = "Energy "
      silo_allvars ( 6) = "XMagneticField "
      silo_allvars ( 7) = "YMagneticField "
      silo_allvars ( 8) = "ZMagneticField "

#if( NFLUID == 1)
      silo_allvars ( NFIELD + 1 + 1) = "Tracer_1 "
#endif

#if( USE_B0B1_FORMULATION == 1)
      silo_allvars ( NFIELD + 1 + NFLUID + 1) = "XTotalMagnField "
      silo_allvars ( NFIELD + 1 + NFLUID + 2) = "YTotalMagnField "
      silo_allvars ( NFIELD + 1 + NFLUID + 3) = "ZTotalMagnField "
#endif

! ----------------------------------------
! Visualizing the extra internal variables.

! --------------------
#if( USE_B0B1_FORMULATION == 1)
!-----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NFIELD + 1 + NFLUID + 3 + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NFIELD + 1 + NFLUID + 3 + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NFIELD + 1 + NFLUID + 3 + 3) = "xflag_entflag "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NFIELD + 1 + NFLUID + 3 + 4) = "xflag_tempfloor "
#endif

! ----------
#else
! ----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NFIELD + 1 + NFLUID + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NFIELD + 1 + NFLUID + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NFIELD + 1 + NFLUID + 3) = "xflag_entprs "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NFIELD + 1 + NFLUID + 4) = "xflag_tempfloor "
#endif

! ----------
#endif
! --------------------


#if ( SILO_EXTRA_DATA >= 5)
      IF ( me == 1 ) THEN
        WRITE( 6, *)" "
        WRITE( 6, *)"Stopping at - Set the silo variable names"
        WRITE( 6, *)"SILO_EXTRA_DATA is expected to be <= 4"
        WRITE( 6, *)"SILO_EXTRA_DATA = ", SILO_EXTRA_DATA
        WRITE( 6, *)" "
        CALL FLUSH( 6)
        ERROR STOP
      END IF ! End of "IF ( me == 1 ) THEN".
#endif

! ----------------------------------------
! Required in "WRITE_SILO" routine.

      filenum = 1

! Set silo file name.

      silofil = "pltf0001.silo"

! ------------------------------
! The below filenames are deprecated relative to silo variables.

! Set file names for the image files. Note that this code just dumps
! out zone-centered variables for imaging. Can't restart from these files.

      rhofil = "rhoa0001"
      vlxfil = "vlxa0001"
      vlyfil = "vlya0001"
      vlzfil = "vlza0001"
      prsfil = "prsa0001"

      bfxfil = "bfxa0001"
      bfyfil = "bfya0001"
      bfzfil = "bfza0001"

      speciefil = "sp010001"

! Set file name for the first dump file.

      dumpfil = "dump0001"

! ----------------------------------------------------------------------

      IF ( me == 1 ) THEN

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Chose physical problem given by iproblem = ",
     1              iproblem

      WRITE( 6, *)" "
      WRITE( 6, *)"iproblem, nsilo_images = ", iproblem, nsilo_images   
      WRITE( 6, *)"rmin, rmax = ", rmin, rmax
      WRITE( 6, *)"timestop = ", timestop
      WRITE( 6, *)" "

      WRITE( 6, *)"omega_passive = ", omega_passive
      WRITE( 6, *)"dipole_B = ", dipole_B
      WRITE( 6, *)"isospeed, Mstar, newtonsc, tempfloor = ",
     1             isospeed, Mstar, newtonsc, tempfloor
      WRITE( 6, *)"iflag_npoles   = ", iflag_npoles
      WRITE( 6, *)"multipole_Bmag = ", multipole_Bmag
      WRITE( 6, *)" "

      CALL FLUSH ( 6)

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"Starting timestep; Documenting constants in the module"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Order of accuracy = ", NUMORDER
      WRITE ( 6, *)"Number of zones theta-phi, Number of ghost cells= ",
     1             ntp_zones, NGhost
      WRITE ( 6, *)"Number of zones in radial direction = ", nr_zones
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1       "n_cc_components, n_cc_modes, n_cc_modes_one_higher = ",
     1        n_cc_components, n_cc_modes, n_cc_modes_one_higher
      WRITE ( 6, *)"n_fc_modes = ", n_fc_modes
      WRITE ( 6, *)" "
      WRITE ( 6, *)"n_cc_bfld_modes = ", n_cc_bfld_modes
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"pi, rhofloor, prsfloor, smallnum = ",
     1              pi, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)"gamma, gasconst, molewt = ",gamma, gasconst, molewt
      WRITE ( 6, *)" "
      WRITE ( 6, *)"cfl_coef = ", cfl_coef
      WRITE ( 6, *)" "

      WRITE ( 6, *)"flatten_coef, divratio_einfeldt_rs = ",
     1              flatten_coef, divratio_einfeldt_rs
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lapidus_threshold, lapidus_coef = ",
     1              lapidus_threshold, lapidus_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"blend_multid_flux, blend_mdhll_to_mdhllc = ",
     1              blend_multid_flux, blend_mdhll_to_mdhllc

      WRITE ( 6, *)"ntstep, ndumpstep, nimgstep, time, timestop = ",
     1              ntstep, ndumpstep, nimgstep, time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of processors used = ", nimages
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      CALL FLUSH ( 6)

      END IF ! End of "IF ( me == 1 ) THEN".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! The subroutine sequence ""BUILD_BASE_PATCHES, GENERATE_GLOBAL_RADIAL_LOG_MESH,
! BUILD_MODE_INDICES, Mesh_Sector_General, FILL_GHOST2D_GLOBAL" initializes
! the global mesh.

! These subroutines only depend on "rmin, rmax, smallnum,
! ntp_zones, nr_zones, NGhost, ntp_patches, nr_patches".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Here we build sibling info. for the base level mesh.

      CALL BUILD_BASE_PATCHES ()

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction. This calls a generic on-processor
! serial code. Therefore, synchronize after exit.

#if ( CHOOSE_RADIAL_MESH == 0)

      CALL GENERATE_GLOBAL_RADIAL_UNIFORM_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

#elif ( CHOOSE_RADIAL_MESH == 1)

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

#elif ( CHOOSE_RADIAL_MESH == 2)

      CALL GENERATE_GLOBAL_RADIAL_STRETCHED_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

#elif ( CHOOSE_RADIAL_MESH == 3)

!      radial_mesh_ratio = 1.23606      ! Geometric ratioing.
!     geom_dr1          = 0.003        ! Size of first zone 
!      nr_geom_zones     = 8            ! number of radial patches with geometric mesh
!      IF ( ( nr_global .NE. 160) .OR. ( NINT( rmax) .NE. 7)) THEN

!      radial_mesh_ratio = 1.021783     ! Geometric ratioing.
!      geom_dr1          = 0.003        ! Size of first zone 
!      nr_geom_zones     = 100            ! number of radial patches with geometric mesh
!      IF ( ( nr_global .NE. 200) .OR. ( NINT( rmax) .NE. 7)) THEN

      radial_mesh_ratio = 1.024813     ! Geometric ratioing.
      geom_dr1          = 0.003        ! Size of first zone 
      nr_geom_zones     = 100            ! number of radial patches with geometric mesh
      IF ( ( nr_global .NE. 200) .OR. ( NINT( rmax) .NE. 10)) THEN
        WRITE( 6, *)" "
        WRITE( 6, *)"Stopping before GENERATE_GLOBAL_RADIAL_HYBRID_MESH"
        WRITE( 6, *)"Mismatch in nr_zones and nr_global or rmax"
        WRITE( 6, *)" "
        CALL FLUSH( 6)
        ERROR STOP
      END IF

      CALL GENERATE_GLOBAL_RADIAL_HYBRID_MESH

     1       ( nr_zones, nr_global, NGhost_global, nr_geom_zones,
     1         rmin, rmax, smallnum, radial_mesh_ratio, geom_dr1,
     1         rb_save_global, dr_save_global, rc_save_global)

#endif

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes, n_cc_modes_one_higher,
     1           n_cc_bfld_modes,

     1           ii_mode, jj_mode, kk_mode,

     1           ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1           ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1           ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz )

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
! Copy in gamma, molewt grid variables and then the modes of U and magnetic field.

!$acc enter data copyin(                                                &
!$acc&    vert_coords_patch, imask_active_patch,                        &

!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,            &

!$acc&    ii_mode, jj_mode, kk_mode,                                    &
!$acc&    ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,               &
!$acc&    ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,               &
!$acc&    ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz )

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Loop over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! Since we are dealing with only one level mesh, each processor will work on
! "npatch_per_proc" number of patches. This variable has already been built.
! We loop over patches here.

! Note that "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! and "iactive_FromProc" to be properly initialized on each processor. This
! was already done via calll to "BUILD_BASE_PATCHES".
! ----------------------------------------------------------------------
! START INITIALIZING GEOMETRY AND PROBLEM  ;  START INITIALIZING GEOMETRY AND PROBLEM

! The following subroutine sequence extracts patch coordinates from the
! global mesh, "EXTRACT_VERT_COORDS_PATCH, STENCIL_INDEXING_TRICORNERS".
! This extraction only needs to be done once and fills in
! "vert_coords_patch, imask_active_patch". Saved in "vert_coords_amr, imask_active_amr".

! From unpacking those two arrays, we can rebuild the entire geometry with the
! sequence, "COORDS_TO_WT_QUADPTS_FACES_ISO2, COORDS_TO_WT_QUADPTS_VOLUME_ISO2,
! INIT_BASIS_FUNCTION_COEFS".

! MHD Data is initialized after geometry with "INIT_PATCH_CUBE_SPHERE".

! ----------------------------------------------------------------------

      dtcurp1_reduce = 1.0e32  ! Set up to globally minimize the first timestep.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

! From "ipatch" (which is the logical patch number in "isector") we can
! extract "ix_patch, iy_patch, iz_patch".

      jj = ipatch
      iz_patch = 1 + (jj - 1) / (ntp_patches**2)

      jj = jj - (iz_patch - 1) * ntp_patches**2
      iy_patch = 1 + (jj - 1) / ntp_patches

      jj = jj - (iy_patch - 1) * ntp_patches
      ix_patch = jj

! ------------------------------

! Fill in "isector" and level information. Also fill in "ix_patch, iy_patch, iz_patch"
! information.
      isector_ilevel_amr ( 1, ii) = isector
      isector_ilevel_amr ( 2, ii) = 1  ! For now, only one level in the AMR hierarchy.

! Fill in "ix_patch, iy_patch, iz_patch" for this patch
      ixyz_patch_amr ( 1, ii) = ix_patch
      ixyz_patch_amr ( 2, ii) = iy_patch
      ixyz_patch_amr ( 3, ii) = iz_patch


! ----------------------------------------------------------------------
! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! ------------------------------
! Save "vert_coords_patch, imask_active_patch" for later reuse. These are
! also needed for imaging.

      vert_coords_amr ( :, :, :, :, ii) =
     1                              vert_coords_patch ( :, :, :, :)

      imask_active_amr ( :, :, :, ii) = imask_active_patch ( :, :, :)

! ------------------------------

! Understand that we are working on a single patch on each processor.
! We have just filled in "vert_coords_patch, imask_active_patch" which is
! the minimal amount of information needed for regenating geometry,
! quadrature points, weights etc.

! ----------------------------------------------------------------------
! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ------------------------------
! Save information on whether corners are tri-corners or quad-corners.
! Also, save stencil information because for cube-sphere this is only 2D
! and very lightweight.

      itri_patch_amr ( 1, ii) = itri_patch_LD
      itri_patch_amr ( 2, ii) = itri_patch_RD
      itri_patch_amr ( 3, ii) = itri_patch_RU
      itri_patch_amr ( 4, ii) = itri_patch_LU

      istencil_amr ( :, :, :, :, ii) = istencil_patch ( :, :, :, :)

! ----------------------------------------------------------------------
! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".
! ----------------------------------------------------------------------
! BEGIN BUILDING GEOMETRY  ;  BEGIN BUILDING GEOMETRY  ;  BEGIN BUILDING GEOMETRY
! ----------------------------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

! ----------------------------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ----------------------------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes, n_cc_modes_one_higher, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------

! Now we save the mimimal amount of geometry information from "_patch" variables
! to "_amr" variables. These are just the variables that are needed for making
! a time-update and timestep evaluation. They are designed to have very
! low storage requirements.

      xf_quadnorml_amr ( :, :, :, :, ii) = 
     1                                xf_quadnorml_patch ( :, :, :, :)

      yf_quadnorml_amr ( :, :, :, :, ii) = 
     1                                yf_quadnorml_patch ( :, :, :, :)

      zf_quadnorml_amr ( :, :, :, :, :, ii) = 
     1                             zf_quadnorml_patch ( :, :, :, :, :)

      xf_area_amr ( :, :, :, ii) = xf_area_patch ( :, :, :)
      yf_area_amr ( :, :, :, ii) = yf_area_patch ( :, :, :)
      zf_area_amr ( :, :, :, ii) = zf_area_patch ( :, :, :)


      v_scale_length_amr ( :, :, :, ii) = v_scale_length ( :, :, :)

      v_volume_amr ( :, :, :, ii) = v_volume_patch ( :, :, :)

      v_centroid_amr ( :, :, :, :, ii) = v_centroid_patch ( :, :, :, :)
      
! ----------------------------------------------------------------------
! END BUILDING GEOMETRY  ;  END BUILDING GEOMETRY  ;  END BUILDING GEOMETRY
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! START INITIALIZE ; START INITIALIZE ; START INITIALIZE ; START INITIALIZE
! ----------------------------------------------------------------------
! By this point, all geometry is built. Therefore, we are now in a position
! to initialize the MHD variables for this patch.

! ------------------------------
! The restart condition is removed, so as to fill the ghost zone information
! of the "B0" field.

!#if ( RESTART == 0)
! ------------------------------

      a = cfl_coef
      cfl_coef = 0.05

!$acc update device ( cfl_coef)

      CALL INIT_PATCH_CUBE_SPHERE

     1         ( nr_zones, ntp_zones, NGhost,                            
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,     
     1           n_fc_modes, n_cc_bfld_modes,                            
                                                                         
     1           iproblem, time, dtcur,                                  
     1           rmin, rmax,                                             
                                                                         
     1           vert_coords_patch, imask_active_patch,                  
     1           v_scale_length,                                         

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch,                        
                                                                         
     1           U, V, entropygr, Bxf, Byf, Bzf,                         
     1           B0field,
                                                                         
     1           gamma_eos_gr,                                           
                                                                         
     1           V_min_local, V_max_local)

      cfl_coef = a  ! Reset the CFL number to its original value.

! Store local "dtcur" for each patch that resides on this processor.
      dtcurp1_amr ( ii) = dtcur

! For all the patches that reside on this processor, make a reduction of their
! "dtcur" value
      dtcurp1_reduce = AMIN1 ( dtcurp1_reduce, dtcur)

! ----------------------------------------------------------------------
! By this point, all "U, V, entropygr, Bxf, Byf, Bzf, gamma_eos_gr" for this
! "ipatch" are initialized in the interior of the mesh. But each processor
! can have multiple patches, which are all stored in
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

      B0field_amr ( :, :, :, 1:3, ii) = B0field ( :, :, :, 1:3)

! ------------------------------
!#endif
! ------------------------------

! ----------------------------------------------------------------------
! END INITIALIZE ; END INITIALIZE ; END INITIALIZE ; END INITIALIZE ;
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
#if ( RESTART == 0)
! ------------------------------

! Now make a reduction across all processors for "dtcurp1_reduce".

      CALL MPI_ALLREDUCE ( dtcurp1_reduce,   ! Input for reudction
     1                     dtcur,            ! Global minimum
     1                     1,                ! length of reduction variable   
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Initialized all patches, dtcur = ", dtcur
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! ------------------------------
#endif
! ------------------------------

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! END INITIALIZING GEOMETRY AND PROBLEM  ;  END INITIALIZING GEOMETRY AND PROBLEM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! End of looping over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! Build, "is_bottom_or_top" information for each patch.
! Used in reading and writing of silo files.

      is_bottom_or_top ( :, :) = 0.0

      DO ii = 1, npatch_per_proc

        iz_patch = ixyz_patch_amr ( 3, ii)

        IF ( iz_patch == 1)           is_bottom_or_top ( 1, ii) = 1
        IF ( iz_patch == nr_patches)  is_bottom_or_top ( 2, ii) = 1

      END DO

! ----------------------------------------------------------------------
! START READING IN RESTART DUMP ; START READING IN RESTART DUMP ;
! ----------------------------------------------------------------------

! Read in a prior dump file if one has been written out. Increment dump
! file name.

! ------------------------------
#if ( RESTART == 1)
! ------------------------------

      IF (me == 1) THEN

      OPEN ( unit = 5, file = "riemann.in")

      READ ( 5, *) time, dtcur, jj
      READ ( 5, *) filenum
      CLOSE ( 5)

      a = time
      b = dtcur

      END IF

! Now broadcast "time, dtcur, filenum" to all processors.
      CALL MPI_Bcast ( time,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast ( dtcur,   1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast ( filenum, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Now read the restart dump.

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Reading dumpfile with filenum = ", filenum
      WRITE ( 6, *)" "
      END IF

! --------------------------------------------------

#if( USE_B0B1_FORMULATION == 1 )

      CALL READDATA ( npatch_per_proc,
     1                iactive_FromProc ( :, 2),
     1                nr_zones, ntp_zones, NGhost,
     1                is_bottom_or_top,
     1                itri_patch_amr,
     1                vert_coords_amr,
     1                n_cc_components + 3 + SILO_EXTRA_DATA,
     1                SILO_VARLEN, silo_allvars,
     1                U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                filenum, time
     1              )

      DO ii = 1, npatch_per_proc

            U_amr ( :, :, :, 1:n_cc_components, ii) =
     1 U_silo_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 2)
      entropygr_amr ( :, :, :, ii) = 
     1   U_silo_amr ( :, :, :, n_cc_components + 3 + 2, ii)
#endif

      END DO

#else

      CALL READDATA ( npatch_per_proc,
     1                iactive_FromProc ( :, 2),
     1                nr_zones, ntp_zones, NGhost,
     1                is_bottom_or_top,
     1                itri_patch_amr,
     1                vert_coords_amr,
     1                n_cc_components + SILO_EXTRA_DATA,
     1                SILO_VARLEN, silo_allvars,
     1                U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                filenum, time
     1              )

      DO ii = 1, npatch_per_proc

            U_amr ( :, :, :, 1:n_cc_components, ii) =
     1 U_silo_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 2)
      entropygr_amr ( :, :, :, ii) = 
     1   U_silo_amr ( :, :, :, n_cc_components + 2, ii)
#endif

      END DO

#endif

! --------------------------------------------------

      IF ( me == 1) THEN

        IF ( ABS ( time - a) >= smallnum) THEN
        WRITE ( 6, *)
     1  "time does not match between the restart dump and riemann.in",
     1  time, a
        CALL FLUSH ( 6)
        STOP
        END IF

        IF ( ABS ( dtcur - b) >= smallnum) THEN
        WRITE ( 6, *)
     1  "dtcur does not match between the restart dump and riemann.in",
     1  dtcur, b
        CALL FLUSH ( 6)
        STOP
        END IF

        IF ( jj .NE. n_cc_components) THEN
        WRITE ( 6, *)"Stopping because : jj .NE. n_cc_components ",
     1                jj, n_cc_components
        CALL FLUSH ( 6)
        STOP
        END IF

      END IF

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Reading of restart file successful !!!"
      WRITE ( 6, *)"Entering time step loop with time, dtcur:",
     1              time, dtcur
      WRITE ( 6, *)" "
      END IF

      timeimage = time

      filenum = filenum + 1

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
#endif
! ----------------------------------------------------------------------
! END READING IN RESTART DUMP ; END READING IN RESTART DUMP ;
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Now that we have set up the initial conditions we put "do_on_gpu = .TRUE.". 
! This allows all the Open ACC loops in all ACC-enabled subroutines to 
! execute on the GPU.   
                                                                         
! So far, we have been completely serial on CPU. Therefore, do one more  
! update of all "(/RIEMANN_COM/)" variables on GPU.                      
                                                                         
      do_on_gpu = .TRUE.                                                 
                                                                         
!$acc update device (/RIEMANN_COM/)
!$acc update device ( time, dtcur, cfl_coef)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! START TIMESTEP LOOP  ;  START TIMESTEP LOOP  ;  START TIMESTEP LOOP
! ----------------------------------------------------------------------

      DO istep = 1, ntstep

      CALL CPU_TIME ( timing_timestep1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Synchronize all ghost zone data for "U_amr, entropygr_amr, Bxf_amr,Byf_amr,Bzf_amr".
! All these subroutines have "CALL MPI_Barrier" at begin and end.
! All these MPI calls are on CPU.

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UX ()

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UY ()

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UY_TOP ()

      CALL LOAD_LOCAL_BUFFERS_UZ ()
      CALL FETCH_BUFFERS_FILL_IN_UZ ()

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! If we get here, we write out a dump. The ghost zones in the theta-phi
! directions are filled in. Not so in the r-direction.

      IF ( MOD ( istep, 2) == 1 ) THEN
      IF ( ( time .GE. timeimage) .OR. ( time .GE. timestop) ) THEN

! ----------------------------------------------------------------------

      IF ( me == 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Writing image file with filenum = ", filenum
      WRITE ( 6, *)" "
      WRITE ( 6, *)"time, dtcur = ", time, dtcur
      WRITE ( 6, *)" "

      CALL FLUSH ( 6)


      OPEN ( unit = 5, file = "riemann.in")

      WRITE ( 5, *) time, dtcur, n_cc_components
      WRITE ( 5, *) filenum

      CALL FLUSH ( 5)
      CLOSE ( 5)

      END IF

! --------------------------------------------------

      CALL CPU_TIME ( timing_write1)

#if( USE_B0B1_FORMULATION == 1 )

      DO ii = 1, npatch_per_proc

        U_silo_amr ( :, :, :, 1:n_cc_components, ii) =
     1       U_amr ( :, :, :, 1:n_cc_components, ii)

        U_silo_amr ( :, :, :, n_cc_components+1:
     1                        n_cc_components+3, ii)=
     1 B0field_amr ( :, :, :, 1:3, ii) + U_silo_amr ( :, :, :, 6:8, ii)


#if ( SILO_EXTRA_DATA >= 1)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 1, ii) =
     1                        divgr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 2)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 2, ii) =
     1                         entropygr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 3)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 3, ii) =
     1                     flag_entprs_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 4)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 4, ii) =
     1                      flag_Tfloor_amr ( :, :, :, ii)
#endif

      END DO

      CALL WRITEDATA ( npatch_per_proc,
     1                 iactive_FromProc ( :, 2),
     1                 nr_zones, ntp_zones, NGhost,
     1                 is_bottom_or_top,
     1                 itri_patch_amr,
     1                 vert_coords_amr,
     1                 n_cc_components + 3 + SILO_EXTRA_DATA,
     1                 SILO_VARLEN, silo_allvars,
     1                 U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                 filenum, time
     1                )

#else

      DO ii = 1, npatch_per_proc

        U_silo_amr ( :, :, :, 1:n_cc_components, ii) =
     1       U_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 1)
        U_silo_amr ( :, :, :, n_cc_components + 1, ii) =
     1                        divgr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 2)
        U_silo_amr ( :, :, :, n_cc_components + 2, ii) =
     1                         entropygr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 3)
        U_silo_amr ( :, :, :, n_cc_components + 3, ii) =
     1                     flag_entprs_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 4)
        U_silo_amr ( :, :, :, n_cc_components + 4, ii) =
     1                      flag_Tfloor_amr ( :, :, :, ii)
#endif

      END DO

      CALL WRITEDATA ( npatch_per_proc,
     1                 iactive_FromProc ( :, 2),
     1                 nr_zones, ntp_zones, NGhost,
     1                 is_bottom_or_top,
     1                 itri_patch_amr,
     1                 vert_coords_amr,
     1                 n_cc_components + SILO_EXTRA_DATA,
     1                 SILO_VARLEN, silo_allvars,
     1                 U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                 filenum, time
     1                )

#endif

      CALL CPU_TIME ( timing_write2)

! --------------------------------------------------

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed writing image:", timing_write2 - timing_write1
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)
      END IF

! ----------------------------------------------------------------------

      filenum = filenum + 1

      timeimage = timeimage + dtimage

      END IF  ! End of "IF ( time .GE. timeimage)"
      END IF  ! End of "IF ( MOD ( istep, 2) == 0 )"

! ----------------------------------------------------------------------
! This is the principal stopping point in this code.

      IF ( time .GE. timestop)  EXIT

! ----------------------------------------------------------------------
! BEGIN LOOP OVER PATCHES ON PROCESSOR   ;   BEGIN LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------
! Now we loop over all "npatch_per_proc" patches on each processor and
! build fluxes and electric fields.

      dtcurp1_reduce = 1.0e32

      V_min_reduce ( :) = 1.0e32
      V_max_reduce ( :) = -1.0e32

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 1) First reverse all that we wrote into the big "*_amr" arrays during
! initialization.

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! --------------------------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! --------------------------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! --------------------------------------------------
! To be efficient on GPUs, we need to store in variables like "U_skinny" and then
! move them to GPU.

      U_skinny ( :, :, :, 1:n_cc_components) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr_skinny ( :, :, :) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf_skinny ( :, :, :) = Bxf_amr ( :, :, :, ii)
      Byf_skinny ( :, :, :) = Byf_amr ( :, :, :, ii)
      Bzf_skinny ( :, :, :) = Bzf_amr ( :, :, :, ii)      

      flattengr ( :, :, :) = flatten_amr ( :, :, :, ii)

! --------------------------------------------------
! Now that all the "_patch" and "_skinny" type data has been updated on CPU, 
! move it to GPU.

!$acc update device ( U_skinny, entropygr_skinny, Bxf_skinny,Byf_skinny,&
!$acc&       Bzf_skinny, flattengr,                                     &

!$acc&       vert_coords_patch, imask_active_patch, istencil_patch,     &

!$acc&       isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch,     &
!$acc&       itri_patch_LD, itri_patch_RD, itri_patch_RU,               &
!$acc&       itri_patch_LU)

! --------------------------------------------------

! Now we copy all the "_skinny" variables back to the actual variables on the GPU

!$acc kernels
      U ( :, :, :, 1:n_cc_components, 1) =
     1       U_skinny ( :, :, :, 1:n_cc_components)
      entropygr ( :, :, :, 1) = entropygr_skinny ( :, :, :)

      Bxf ( :, :, :, 1) = Bxf_skinny ( :, :, :)
      Byf ( :, :, :, 1) = Byf_skinny ( :, :, :)
      Bzf ( :, :, :, 1) = Bzf_skinny ( :, :, :)
!$acc end kernels

! ----------------------------------------------------------------------
! STEP 2) Now build all the rest of the geometry information for this patch.

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! --------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes, n_cc_modes_one_higher, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! STEP 3) With geometry built, now do limiting, predictor step, flux
! and field evaluation.

! --------------------------------------------------
! Pad the boundaries in the radial direction if this is the innermost
! radial patch or the outermost radial patch.

! ------------------------------
! Most of the geomesh simulations require to implement the outflow  boundary
! conditions at the outer boundary. The PAD_BOUNDARY subroutine does that.
! "2" - outflow boundary condition. "17" - no boundary condition imposed.

      bcarr_r ( :, :, 1) = 17   ! Inner boundary.
      bcarr_r ( :, :, 2) = 17   ! Outer boundary.

      IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17     ! Inner boundary.
      IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 17     ! Outer boundary.


#if( CHOOSE_APPL_CAK == 1)
      IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17     ! Inner boundary.
      IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 2      ! Outer boundary.
#endif

!$acc update device ( bcarr_r)

      IF ( ( iz_patch .EQ. 1) .OR. ( iz_patch .EQ. nr_patches)) THEN

      CALL PAD_BOUNDARY_R

     1    ( nr_zones, ntp_zones, NGhost,
     1      n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1      n_fc_modes, n_cc_bfld_modes,

     1      vert_coords_patch, imask_active_patch,

     1      xf_area_patch, yf_area_patch, zf_area_patch,
     1      v_volume_patch,

     1      v_centroid_patch,

     1      U, entropygr, Bxf, Byf, Bzf,

     1      bcarr_r)

      END IF

! ------------------------------
! Most of the geomesh simulations require the re-initialization of boundary
! conditions at the inner boundary. The ADD_BOUNDARY subroutine does that.

      ibdy_inner = 0
      ibdy_outer = 0

      IF ( iz_patch .EQ. 1         )  ibdy_inner = 1   ! 1 == Impose Add_boundary
      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 1   ! 1 == Impose Add_boundary


#if( CHOOSE_APPL_CAK == 1)
      IF ( iz_patch .EQ. 1         )  ibdy_inner = 1   ! 1 == Impose Add_boundary
      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 0   ! 1 == Impose Add_boundary
#endif

!$acc update device ( ibdy_inner, ibdy_outer )

      IF ( ( iz_patch .EQ. 1) .OR. ( iz_patch .EQ. nr_patches)) THEN

      CALL ADD_BOUNDARY

     1         ( ibdy_inner, ibdy_outer,

     1           nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_centroid_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           gamma_eos_gr,

     1           V_min_local, V_max_local)

      END IF

! ----------------------------------------------------------------------
! Part of STEP 3) also as STEP 5) for even time steps.
! Evolve the internal energy for the cooling function source term;
! using the operator-splitting. This is a stiff source term. It is evolved 
! using the adaptive-stepsize RK scheme.

! For odd-steps, call the operator split-source after the MHD evolution.

#if ( USE_COOLING_FUNCTION == 1)

      IF ( MOD ( istep, 2) == 1) THEN

      CALL OPERATOR_SPLIT_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,      

     1           U, entropygr, Bxf, Byf, Bzf, flag_Tfloor)

#if ( SILO_EXTRA_DATA >= 4)
!$acc update host( flag_Tfloor )
      flag_Tfloor_amr ( :, :, :, ii) = flag_Tfloor ( :, :, :)
#endif

      END IF

#endif

! --------------------------------------------------

! Apply limiter to MHD variables.

      CALL LIMIT_ALL_SLOPES_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           iz_patch, nr_patches,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! Store local min and max for each patch that resides on this processor.
      V_min_amr ( :, ii) = V_min_local ( :)
      V_max_amr ( :, ii) = V_max_local ( :)

! For all the patches that reside on this processor, make a reduction of their
! min and max values.
      V_min_reduce = AMIN1 ( V_min_reduce, V_min_local)
      V_max_reduce = AMAX1 ( V_max_reduce, V_max_local)

! ----------
! divergence, stored in "workgr"

#if ( SILO_EXTRA_DATA >= 1)

!$acc update host( workgr )
      divgr_amr ( :, :, :, ii) = workgr ( :, :, :)

#endif

! --------------------------------------------------

!$acc kernels
      dU_dt_predgr    = 0.0
      dEnt_dt_predgr  = 0.0
      U_TAvg_modes_gr = 0.0
!$acc end kernels

! Call ADER predictor step.

      CALL PREDICTOR_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           iproblem, time, dtcur,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc, Source,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         dU_dt_predgr, dEnt_dt_predgr, U_TAvg_modes_gr)


! --------------------------------------------------

! Call routines that build the fluxes and electric fields.

      CALL MAKE_FACIAL_FLUXES_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         dU_dt_predgr, dEnt_dt_predgr, U_TAvg_modes_gr,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,

     1         Flux_x, Flux_y, Flux_z,
     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz)

! --------------------------------------------------

! Call routine that builds the electric fields via 2Ddimensional Riemann solver.

#if ( USE_FACIAL_EFIELDS != 1)

      CALL MAKE_EDGE_EFIELDS_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch, 
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, 
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, 
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU, 
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         dU_dt_predgr, dEnt_dt_predgr, U_TAvg_modes_gr,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,

     1         Flux_x, Flux_y, Flux_z,
     1         Efx, Efy, Efz)

#endif

! ----------------------------------------------------------------------
! STEP 4) With geometry built, and with all fluxes in-hand, carry out
! the update. Note that these fluxes (or electric fields) have not been
! synchronized across patch boundaries. So the update is not perfect and we
! willl have to make it perfect with a defect correction step.

      dt_defcor = dtcur         ! backup the dtcur value for defect-corrector step.
      time_defcor = time        ! backup the time  value for defect-corrector step.

!$acc update device ( time, dtcur, dt_defcor, time_defcor )

      CALL UPDATE_RK_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadnorml_patch, yf_quadnorml_patch,
     1           zf_quadnorml_patch,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, Source,

     1         Flux_x, Flux_y, Flux_z, EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz, flag_entprs)


#if ( SILO_EXTRA_DATA >= 3)
!$acc update host( flag_entprs )
      flag_entprs_amr ( :, :, :, ii) = flag_entprs ( :, :, :)
#endif

! ----------------------------------------------------------------------
! STEP 5) Evolve the internal energy for the cooling function source term;
! using the operator-splitting. This is a stiff source term. It is evolved 
! using the adaptive-stepsize RK scheme.

! For even-steps, call the operator split-source after the MHD evolution.

#if ( USE_COOLING_FUNCTION == 1)

      IF ( MOD ( istep, 2) == 0) THEN

      CALL OPERATOR_SPLIT_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,    
     1           n_fc_modes, n_cc_bfld_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,      

     1           U, entropygr, Bxf, Byf, Bzf, flag_Tfloor)

#if ( SILO_EXTRA_DATA >= 4)
!$acc update host( flag_Tfloor )
      flag_Tfloor_amr ( :, :, :, ii) = flag_Tfloor ( :, :, :)
#endif

      END IF

#endif

! ----------------------------------------------------------------------
! STEP 6) Evaluate a new timestep for this patch. Store "dtcurp1" so that we can
! globally reduce it after this loop is done on ALL processors.

      CALL EVALUATE_TIMESTEP

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           time, dtcur, dtcurp1,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr)

! Store local "dtcurp1" for each patch that resides on this processor.
! Inside the timestep subroutine, "dtcurp1" has been moved back from GPU to CPU.
      dtcurp1_amr ( ii) = dtcurp1

! For all the patches that reside on this processor, make a reduction of their
! "dtcurp1" value.
      dtcurp1_reduce = AMIN1 ( dtcurp1_reduce, dtcurp1)


! ----------------------------------------------------------------------
! STEP 7) "U, entropygr, Bxf, Byf, Bzf" have been almost perfectly updated on GPU.
! So, copy them back to CPU. The only fluxes we need are the ones that fill in
! the faces for the defect correction. Therefore, we have an extra subroutine
! on the GPU that copies "Flux_snd_r, l, u, d, t, b_patch" back to CPU.

! First, we load all the "_skinny" variables; Second, use the flux and 
! electric field arrays on GPU to load the "Flux_snd_r_patch" etc. on the GPU. 
! Only then can we efficiently bring them back to CPU.

!$acc kernels
      U_skinny ( :, :, :, 1:n_cc_components) = 
     1       U ( :, :, :, 1:n_cc_components, 1)

      entropygr_skinny ( :, :, :) = entropygr ( :, :, :, 1)

      Bxf_skinny ( :, :, :) = Bxf ( :, :, :, 1)
      Byf_skinny ( :, :, :) = Byf ( :, :, :, 1)
      Bzf_skinny ( :, :, :) = Bzf ( :, :, :, 1)
!$acc end kernels


      CALL FILL_FLUX_SND_PATCH

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           Flux_x, Flux_y, Flux_z,
     1           EntFlux_x, EntFlux_y, EntFlux_z,
     1           Efx, Efy, Efz,

     1           Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,
     1           Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch)


! ----------------------------------------------------------------------
!$acc update host ( U_skinny, entropygr_skinny,                         &
!$acc&              Bxf_skinny, Byf_skinny, Bzf_skinny,                 &
!$acc&      Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,       &
!$acc&      Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch)

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1           U_skinny ( :, :, :, 1:n_cc_components)

      entropygr_amr ( :, :, :, ii) = entropygr_skinny ( :, :, :)

      Bxf_amr ( :, :, :, ii) = Bxf_skinny ( :, :, :)
      Byf_amr ( :, :, :, ii) = Byf_skinny ( :, :, :)
      Bzf_amr ( :, :, :, ii) = Bzf_skinny ( :, :, :)

      Flux_snd_r_amr ( :, :, :, ii) = Flux_snd_r_patch ( :, :, :)
      Flux_snd_l_amr ( :, :, :, ii) = Flux_snd_l_patch ( :, :, :)

      Flux_snd_u_amr ( :, :, :, ii) = Flux_snd_u_patch ( :, :, :)
      Flux_snd_d_amr ( :, :, :, ii) = Flux_snd_d_patch ( :, :, :)

      Flux_snd_t_amr ( :, :, :, ii) = Flux_snd_t_patch ( :, :, :)
      Flux_snd_b_amr ( :, :, :, ii) = Flux_snd_b_patch ( :, :, :)

! ----------------------------------------------------------------------
! We may never need this, but it helps to have the lines to fill in the fluxes.
! If bandwidth to main memory ever improves for GPUs, this would be an excellent step.

#if ( 0 == 1)

!$acc update host ( Flux_x, Flux_y, Flux_z, EntFlux_x, EntFlux_y,       &
!$acc&              EntFlux_z, Efx, Efy, Efz)

      Flux_x_amr ( :, :, :, :, ii) = Flux_x ( :, :, :, :)
      Flux_y_amr ( :, :, :, :, ii) = Flux_y ( :, :, :, :)
      Flux_z_amr ( :, :, :, :, ii) = Flux_z ( :, :, :, :)

      EntFlux_x_amr ( :, :, :, ii) = EntFlux_x ( :, :, :)
      EntFlux_y_amr ( :, :, :, ii) = EntFlux_y ( :, :, :)
      EntFlux_z_amr ( :, :, :, ii) = EntFlux_z ( :, :, :)

      Efx_amr ( :, :, :, ii) = Efx ( :, :, :)
      Efy_amr ( :, :, :, ii) = Efy ( :, :, :)
      Efz_amr ( :, :, :, ii) = Efz ( :, :, :)

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Let all patches synchronize before variable min, max reduction.
! That way, they have all updated their "V_min_reduce, V_max_reduce".


      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! END LOOP OVER PATCHES ON PROCESSOR   ;   END LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------

! STEP 8) Now make a reduction across all processors for "V_min_reduce, V_max_reduce".

      V_min_local = V_min_reduce
      V_max_local = V_max_reduce

      CALL MPI_ALLREDUCE ( V_min_local,          ! Input for reudction
     1                     V_min_reduce,         ! Global minimum
     1                     n_cc_components + 1,  ! length of reduction variable   
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE ( V_max_local,          ! Input for reudction
     1                     V_max_reduce,         ! Global maximum
     1                     n_cc_components + 1,  ! length of reduction variable   
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

! Write out the global min and max of all variables from processor 0.

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Step # = ", istep," time = ", time,
     1             " dtcur = ", dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V_min = ", V_min_reduce ( :)
      WRITE ( 6, *)"V_max = ", V_max_reduce ( :)
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! STEP 9) Obtain globally minimized "dtcurp1_reduce" and broadcast it to
! all "dtcurp1". Also update "time".

! Let all patches synchronize before timestep reduction.
! That way, they have all updated their "dtcurp1_reduce".
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Now make a reduction across all processors for "dtcurp1_reduce".

      CALL MPI_ALLREDUCE ( dtcurp1_reduce,   ! Input for reudction
     1                     dtcurp1,          ! Global minimum
     1                     1,                ! length of reduction variable
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      CALL CPU_TIME ( timing_timestep2)

      time = time + dtcur

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Step # = ", istep," found dtcurp1 = ", dtcurp1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed this timestep = ", timing_timestep2
     1                                     - timing_timestep1
      WRITE ( 6, *)" "

      END IF

! ----------------------------------------

#if ( ACCURACY_ANALYSIS == 1)
      CALL MAKE_GLOBAL_ACCURACY_CS ( iproblem, rmin, rmax, time, dtcur)
#endif

! ----------------------------------------

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "

      END IF

! ----------------------------------------------------------------------
! Make the code to stop at the "time == timestop".
      IF ( (time + dtcurp1) .GT. timestop)
     1   dtcurp1 = MAX ( ( timestop - time), smallnum / 100.0)

      dtcur = dtcurp1

!$acc update device ( time, dtcur, dtcurp1)

! ----------------------------------------------------------------------

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 10) Now recall that all the fluxes are slightly defective. To make
! them perfect, we have to synchronize flux and e-field data across patch boundaries.

! Synchronize all flux & e-field data for 
! "Flux_x_amr, Flux_y_amr, Flux_z_amr, EntFlux_x_amr, EntFlux_y_amr, EntFlux_z_amr,
! Efx_amr, Efy_amr, Efz_amr".
! All these subroutines have "CALL MPI_Barrier" at begin and end.

! Once we come here, all defects should start off as zeros, even if they are
! filled in at appropriate locations in the course of this step.

      Defect_Flux_x_amr = 0.0
      Defect_Flux_y_amr = 0.0
      Defect_Flux_z_amr = 0.0

      Defect_EntFlux_x_amr = 0.0
      Defect_EntFlux_y_amr = 0.0
      Defect_EntFlux_z_amr = 0.0

      Defect_Efx_amr = 0.0
      Defect_Efy_amr = 0.0
      Defect_Efz_amr = 0.0


! --------------------------------------------------
! STEP 10.a) Fill the boundary values of "Flux_xyz_amr, Defect_Flux_xyz_amr"

      CALL LOAD_AMR_FLUXES_AT_BDY ( )


! --------------------------------------------------
! STEP 10.a) Synchronize all flux & e-field data

      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()

      CALL FETCH_BUFFERS_FILL_IN_Flux_X ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Y ()

      CALL LOAD_LOCAL_BUFFERS_Flux_Z ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Z ()

      CALL SYNCHRONIZE_EFZ_QUAD_CORNER ()


! --------------------------------------------------
! STEP 10.c) Unpack the defect fluxes and defect electric fields.

      CALL UNPACK_DEFECT_FLUXES ( )


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 11) The defects have been built in the previous step.
! This is the defect correction step. It operates just like a timestep.
! The defects have to be stored as if they are fluxes, which has already
! been done in the previous step. The defects only affect the boundary zones.
! This is such a low cost operation that it can be done entirely on CPU.

! ----------------------------------------------------------------------
! BEGIN LOOP OVER PATCHES ON PROCESSOR   ;   BEGIN LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------
! Now we loop over all "npatch_per_proc" patches on each processor and
! just correct the defects.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! --------------------------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! --------------------------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! --------------------------------------------------
! This is entirely on CPU, so we can copy from "U_amr" to "U" and so on.
! Also copy in the minimal amount of geometry variables, because we will
! not be truly rebuilding geometry here.

      U ( :, :, :, 1:n_cc_components, 1) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr ( :, :, :, 1) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)


      Source ( :, :, :, :) = 0.0     ! NO Source for defect-correction step.

      xf_quadnorml_patch ( :, :, :, :) =
     1                              xf_quadnorml_amr ( :, :, :, :, ii)

      yf_quadnorml_patch ( :, :, :, :) =
     1                              yf_quadnorml_amr ( :, :, :, :, ii)

      zf_quadnorml_patch ( :, :, :, :, :) =
     1                           zf_quadnorml_amr ( :, :, :, :, :, ii)


      xf_area_patch ( :, :, :) = xf_area_amr ( :, :, :, ii)
      yf_area_patch ( :, :, :) = yf_area_amr ( :, :, :, ii)
      zf_area_patch ( :, :, :) = zf_area_amr ( :, :, :, ii)


      v_scale_length ( :, :, :) = v_scale_length_amr ( :, :, :, ii)

      v_volume_patch ( :, :, :) = v_volume_amr ( :, :, :, ii)

      v_centroid_patch ( :, :, :, :) = v_centroid_amr ( :, :, :, :, ii)

! --------------------------------------------------
! The defects are used in the update step just like fluxes. Therefore,
! copy in "Defect_Flux_x_amr" etc. into "Flux_x" etc. Likewise, copy in
! "Defect_Efx_amr" etc. into "Efx" etc.

      Flux_x ( :, :, :, :) = Defect_Flux_x_amr ( :, :, :, :, ii)
      Flux_y ( :, :, :, :) = Defect_Flux_y_amr ( :, :, :, :, ii)
      Flux_z ( :, :, :, :) = Defect_Flux_z_amr ( :, :, :, :, ii)

      EntFlux_x ( :, :, :) = Defect_EntFlux_x_amr ( :, :, :, ii)
      EntFlux_y ( :, :, :) = Defect_EntFlux_y_amr ( :, :, :, ii)
      EntFlux_z ( :, :, :) = Defect_EntFlux_z_amr ( :, :, :, ii)

      Efx ( :, :, :) = Defect_Efx_amr ( :, :, :, ii)
      Efy ( :, :, :) = Defect_Efy_amr ( :, :, :, ii)
      Efz ( :, :, :) = Defect_Efz_amr ( :, :, :, ii)


! --------------------------------------------------
! Defect correction step. Here, use the previously stored "time_defcor, dt_defcor".

      CALL UPDATE_DEFECT

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes, time_defcor, dt_defcor,  ! note _defcor.

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadnorml_patch, yf_quadnorml_patch,
     1           zf_quadnorml_patch,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, Source,

     1         Flux_x, Flux_y, Flux_z,
     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz)

! --------------------------------------------------
! Now reverse the process; i.e. copy "U" into "U_amr" etc.

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! Let all patches synchronize before starting the next timestep.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! END LOOP OVER PATCHES ON PROCESSOR   ;   END LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO istep = 1, ntstep".

! ----------------------------------------------------------------------
! END TIMESTEP LOOP  ;  END TIMESTEP LOOP  ;  END TIMESTEP LOOP 
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! We created a lot of temporary arrays that are used for processing a timestep
! on the GPU. Since we are at the end of the timestep, we can get rid of 
! these arrays on the GPU.                                               
                                                                         
!$acc end data 
! --------------------------------------------------------------------- 
! All the messaging is done; about to stop code. Now free up all the handles.

      CALL WINFREE_MESSAGING ()

      CALL MPI_FINALIZE (ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END PROGRAM RIEMANN_GEOMESH_CS_MHD

#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------



#if (1 == 0)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM TEST_RM_CS_MODULE

! ----------------------------------------------------------------------
! Written by : Dinshaw Balsara and Sethupathy S.

      USE RIEMANN_3D_CS_AMR_SUPPORT

      IMPLICIT NONE

#include "riemann.com"

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, isector, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, ix_patch, iy_patch, iz_patch

      INTEGER nRK_stage, irk, iproblem, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum

      REAL rmin, rmax, time, dtcur, dtcurp1, dtmove, a,
     1     timeimage, dtimage, timestop

!$acc declare create ( time, dtcur, dtcurp1, dtmove )

      REAL sector_vertices ( 4, 3, 6)

! ----------------------------------------------------------------------
! File names used for imaging.

      CHARACTER*8 rhofil, vlxfil, vlyfil, vlzfil, prsfil,
     1            bfxfil, bfyfil, bfzfil, dumpfil,
     1            speciefil

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)::
     1   r4gr

! ----------------------------------------------------------------------

      REAL timing_info1, timing_info2

! ----------------------------------------------------------------------
! Variables used in testing the MPI messaging.

      INTEGER iproc

! ----------------------------------------------------------------------
! Initialize MPI here.

      CALL MPI_INIT (ierr)

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

! Create windows for message passing.
      CALL WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

! Also create fluid variable arrays that are needed by the algorithm and only need
! to reside on the GPU. These arrays are not directly used for the update
! and will never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,&
!$acc&           xf_area_patch, xf_centroid_patch,                      &

!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &

!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &

!$acc&           v_scale_length,                                        &

!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &

!$acc&           Coef_Basis_Moments,                                    &
!$acc&           itri_patch_LD, itri_patch_RD, itri_patch_RU,           &
!$acc&           itri_patch_LU, istencil_patch,                         &

!$acc&         U, entropygr, Bxf, Byf, Bzf,                             &
!$acc&         V, B1_zc, B2_zc, B3_zc,                                  &

!$acc&         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,       &
!$acc&         divvelgr, msonicgr, flattengr, workgr,                   &
!$acc&         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                   &

!$acc&         V_min_local, V_max_local,                                &

!$acc&         Flux_x, Flux_y, Flux_z,                                  &
!$acc&         EntFlux_x, EntFlux_y, EntFlux_z,                         &
!$acc&         Efx, Efy, Efz)

! ----------------------------------------------------------------------
! Initializing prs and rho floors, inner and outer radius and fluid variables.

      ilevel = 1  ! By default, we work on 1 level.

      OPEN ( UNIT = 6, FILE = "riemann_cs.out")

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

      ix_patch = 2
      iy_patch = 2
      iz_patch = 1

      isector = 1

      rmin = 10.0
      rmax = 15.23

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4     ! Polytropic index for the gas.
      gasconst_euler = 1.0  ! Set up scaled gas constant.
      molewt_euler = 1.4    ! Set up mole wt. of the gas.

      time = 0.0  ! Initial time is always zero.

! Decide which physical problem you want to initialize with "iproblem".
! iproblem = 0 for testing the MPI messaging;
! iproblem = 1 for method of manufactored solution; iproblem = 2 is blast;
! iproblem = 3 is hydro vortex; iproblem = 4 is mhd vortex;
! iproblem = 5 is accuracy test; iproblem = 6 is rotor test;
! iproblem = 7 is Parker wind test; iproblem = 8 hydrodynamic manufactored solution;
! iproblem = 9 is SOD test; iproblem = 10 is Brio-Wu test.

      iproblem = 0

! ----------------------------------------------------------------------

! These parameters control the numerics.

#if ( NUMORDER == 2)
      nRK_stage = 2  ! Two stage RK at second order.
#elif ( NUMORDER == 3)
      nRK_stage = 3  ! Three stage RK at third order.
#elif ( NUMORDER == 4)
      nRK_stage = 5  ! Four stage RK at fourth order.
#endif

! Courant number
      cfl_coef = 0.25

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

! Using 100% multid flux can sometimes mean we take flux contributions from
! a stencil that does not depend on the zones that contribute to that face.
! For that reason, we blend multid RS with 1d RS associated with that face.
! Blending in more multidimensionality also seems to help out with pressure
! positivity without any artificial fixes.
! This value should range between [0.4, 0.6] for good performance.

! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
! flux from each edge.
! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
! The electric field is always fully evaluated by the multid. RS.

      blend_multid_flux = 0.1666

! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
! (almost or fully) stationary contacts. To prevent that we blend in a small
! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]

! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
! reverse in smooth flow.

      blend_mdhll_to_mdhllc = 0.0

! ----------------------------------------------------------------------

! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 50000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.
      dtimage = 0.2 / 40.0

! Time at which you want the simulation to stop.
      timestop = 0.5

! ----------------------------------------------------------------------
! By this time we have updated a lot of variables on the CPU. Let us
! migrate that data to the GPU. We still want to continue working
! on CPU, so we set "do_on_gpu = .FALSE.".

      do_on_gpu = .FALSE.

!$acc update device  (/RIEMANN_COM/)

! ----------------------------------------------------------------------

! Set default multispecies values, just to make sure they are not NANS.

      gamma_specie ( :) = gamma
      molewt_specie ( :) = molewt
      heat_of_formation ( :) = 0.0
      cp_specie ( :) = gasconst * gamma_specie ( :)
     1            / ( ( gamma_specie ( :) - 1.0) * molewt_specie ( :) )
      cv_specie ( :) = cp_specie ( :) / gamma_specie ( :)

! ----------

      gamma_specie_euler ( :) = gamma_euler
      molewt_specie_euler ( :) = molewt_euler
      heat_of_formation_euler ( :) = 0.0
      cp_specie_euler ( :) = gasconst_euler * gamma_specie_euler ( :)
     1 / ( ( gamma_specie_euler ( :) - 1.0) * molewt_specie_euler ( :))
      cv_specie_euler ( :) = cp_specie_euler ( :)
     1                     / gamma_specie_euler ( :)

! ----------
! Make a default initialization of big fluid arrays. This is still done
! on CPU.

      gamma_eos_gr ( :, :, :) = gamma
      gamma_soundspeed_gr ( :, :, :) = gamma
      mean_mol_wt_gr ( :, :, :) = molewt

! ----------------------------------------------------------------------

! Set the variable names used in the silo file. Not done yet.

! Set file names for the image files. Note that this code just dumps
! out zone-centered variables for imaging. Can't restart from these files.

      rhofil = "rhoa0001"
      vlxfil = "vlxa0001"
      vlyfil = "vlya0001"
      vlzfil = "vlza0001"
      prsfil = "prsa0001"

      bfxfil = "bfxa0001"
      bfyfil = "bfya0001"
      bfzfil = "bfza0001"

      speciefil = "sp010001"

! Set file name for the first dump file.

      dumpfil = "dump0001"

! ----------------------------------------------------------------------

      IF ( me == 1 ) THEN

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Chose physical problem given by iproblem = ",
     1              iproblem
      WRITE ( 6, *)" "


      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"Starting timestep; Documenting constants in the module"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Order of accuracy = ", NUMORDER
      WRITE ( 6, *)"Number of zones theta-phi, Number of ghost cells= ",
     1             ntp_zones, NGhost
      WRITE ( 6, *)"Number of zones in radial direction = ", nr_zones
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Picking one of the patches from the sector"
      WRITE ( 6, *)"isector = ", isector
      WRITE ( 6, *)"ix_patch, iy_patch, iz_patch = ",
     1              ix_patch, iy_patch, iz_patch
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1       "n_cc_components, n_cc_modes, n_cc_modes_one_higher = ",
     1        n_cc_components, n_cc_modes, n_cc_modes_one_higher
      WRITE ( 6, *)"n_fc_modes = ", n_fc_modes
      WRITE ( 6, *)" "
      WRITE ( 6, *)"n_cc_bfld_modes = ", n_cc_bfld_modes
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"pi, rhofloor, prsfloor, smallnum = ",
     1              pi, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)"gamma, gasconst, molewt = ",gamma, gasconst, molewt
      WRITE ( 6, *)" "
      WRITE ( 6, *)"cfl_coef = ", cfl_coef
      WRITE ( 6, *)" "

      WRITE ( 6, *)"flatten_coef, divratio_einfeldt_rs = ",
     1              flatten_coef, divratio_einfeldt_rs
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lapidus_threshold, lapidus_coef = ",
     1              lapidus_threshold, lapidus_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"blend_multid_flux, blend_mdhll_to_mdhllc = ",
     1              blend_multid_flux, blend_mdhll_to_mdhllc

      WRITE ( 6, *)"ntstep, ndumpstep, nimgstep, time, timestop = ",
     1              ntstep, ndumpstep, nimgstep, time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of RK stages, nRK_stage = ", nRK_stage
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

      END IF ! End of "IF ( me == 1 ) THEN".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Here we build sibling info. for the base level mesh.

      CALL BUILD_BASE_PATCHES ()

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Write out the iactive, isibling_* data.

!      CALL TEST_1_PRINT_IACTIVE_ISIBLING ( )

! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction. This calls a generic on-processor
! serial code. Therefore, synchronize after exit.

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         rb_save_global, dr_save_global, rc_save_global)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes, n_cc_modes_one_higher,
     1           n_cc_bfld_modes,

     1           ii_mode, jj_mode, kk_mode,

     1           ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1           ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1           ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz )

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Loop over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! Since we are dealing with only one level mesh, each processor will work on
! "npatch_per_proc" number of patches. This variable has already been built.
! We loop over patches here.

! Note that "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! and "iactive_FromProc" to be properly initialized on each processor. This
! was already done via calll to "BUILD_BASE_PATCHES".
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

! From "ipatch" (which is the logical patch number in "isector") we can
! extract "ix_patch, iy_patch, iz_patch".

      jj = ipatch
      iz_patch = 1 + (jj - 1) / (ntp_patches**2)

      jj = jj - (iz_patch - 1) * ntp_patches**2
      iy_patch = 1 + (jj - 1) / ntp_patches

      jj = jj - (iy_patch - 1) * ntp_patches
      ix_patch = jj

! ------------------------------

! Fill in "isector" and level information. Also fill in "ix_patch, iy_patch, iz_patch"
! information.
      isector_ilevel_amr ( 1, ii) = isector
      isector_ilevel_amr ( 2, ii) = 1  ! For now, only one level in the AMR hierarchy.

! Fill in "ix_patch, iy_patch, iz_patch" for this patch
      ixyz_patch_amr ( 1, ii) = ix_patch
      ixyz_patch_amr ( 2, ii) = iy_patch
      ixyz_patch_amr ( 3, ii) = iz_patch

! ------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! ------------------------------
! Save "vert_coords_patch, imask_active_patch" for later reuse. These are
! also needed for imaging.

      vert_coords_amr ( :, :, :, :, ii) =
     1                              vert_coords_patch ( :, :, :, :)

      imask_active_amr ( :, :, :, ii) = imask_active_patch ( :, :, :)

! ------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Understand that we are working on a single patch on each processor.
! We have just filled in "vert_coords_patch, imask_active_patch" which is
! the minimal amount of information needed for regenating geometry,
! quadrature points, weights etc.

! ----------------------------------------------------------------------
! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ------------------------------
! Save information on whether corners are tri-corners or quad-corners.
! Also, save stencil information because for cube-sphere this is only 2D
! and very lightweight.

      itri_patch_amr ( 1, ii) = itri_patch_LD
      itri_patch_amr ( 2, ii) = itri_patch_RD
      itri_patch_amr ( 3, ii) = itri_patch_RU
      itri_patch_amr ( 4, ii) = itri_patch_LU

      istencil_amr ( :, :, :, :, ii) = istencil_patch ( :, :, :, :)

! ------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".
! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes, n_cc_modes_one_higher, n_cc_bfld_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,
     1         ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,
     1         ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! By this point, all geometry is built. Therefore, we are now in a position
! to initialize the MHD variables for this patch.

      a = cfl_coef
      cfl_coef = 0.05

!$acc update device ( cfl_coef)

      CALL TEST_MPI_INIT_PATCH_CUBE_SPHERE

     1         ( nr_zones, ntp_zones, NGhost,                            
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,     
     1           n_fc_modes, n_cc_bfld_modes,                            
                                                                         
     1           iproblem, time, dtcur,                                  
     1           rmin, rmax,                                             
                                                                         
     1           vert_coords_patch, imask_active_patch,                  
     1           v_scale_length,                                         
                                                                         
     1           xf_area_patch, yf_area_patch, zf_area_patch,            
     1           v_quadpts_patch, v_quadwt_patch,                        
                                                                         
     1           U, V, entropygr, Bxf, Byf, Bzf,                         

     1           Flux_x, Flux_y, Flux_z, 
     1           EntFlux_x, EntFlux_y, EntFlux_z, Efx, Efy, Efz,
                                                                         
     1           gamma_eos_gr,                                           
                                                                         
     1           V_min_local, V_max_local)

      cfl_coef = a  ! Reset the CFL number to its original value.

! ----------------------------------------------------------------------
! By this point, all "U, V, entropygr, Bxf, Byf, Bzf, gamma_eos_gr" for this
! "ipatch" are initialized in the interior of the mesh. But each processor
! can have multiple patches, which are all stored in
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

      Flux_x_amr ( :, :, :, :, ii) = Flux_x ( :, :, :, :)
      Flux_y_amr ( :, :, :, :, ii) = Flux_y ( :, :, :, :)
      Flux_z_amr ( :, :, :, :, ii) = Flux_z ( :, :, :, :)

      EntFlux_x_amr ( :, :, :, ii) = EntFlux_x ( :, :, :)
      EntFlux_y_amr ( :, :, :, ii) = EntFlux_y ( :, :, :)
      EntFlux_z_amr ( :, :, :, ii) = EntFlux_z ( :, :, :)

      Efx_amr ( :, :, :, ii) = Efx ( :, :, :)
      Efy_amr ( :, :, :, ii) = Efy ( :, :, :)
      Efz_amr ( :, :, :, ii) = Efz ( :, :, :)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! End of looping over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Test the XY messaging
!       CALL TEST_XY_MESSAGING ()   

! Calculate error in Ghost zones before MPI messaging.

      iproc = 1

      iproc = MAX ( iproc, 1)
      iproc = MIN ( iproc, nimages)

      IF ( me == iproc ) THEN
        WRITE ( 6, *)"iproc, nimages, npatch_per_proc = ",
     1                iproc, nimages, npatch_per_proc
        WRITE ( 6, *)"----------------------------------------------"
        WRITE ( 6, *)"Before MPI messaging:"
        WRITE ( 6, *)" "
      END IF

      CALL TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc,     
     1                 time, dtcur, rmin, rmax )

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Do messaging and fill in Ghost zones.

! --------------------
! Fill ghoszt zones and faces in X and Y directions; i.e., along the
! spherical surface.

! Fill send buffers in both X and Y directions.
      CALL LOAD_LOCAL_BUFFERS_UXY ()                                   
                                                                         
! Fetch the X-data from left and right patches.
      CALL FETCH_BUFFERS_FILL_IN_UX () 

! Fill send buffers in both X and Y directions.
      CALL LOAD_LOCAL_BUFFERS_UXY ()                                   
                                                                         
! Fetch the Y-data from up and down patches.
      CALL FETCH_BUFFERS_FILL_IN_UY ()                                 

! ----------
! Till above all the ghost zones are filled for interior and the bottom
! patches in the sector. For the top patches, we need to call the filling,
! fecthing subroutines in Y-direction again.

      CALL LOAD_LOCAL_BUFFERS_UXY () 
                                                                         
      CALL FETCH_BUFFERS_FILL_IN_UY_TOP ()

! --------------------
! Now message in the Z direction. This completes the messaging.

! Fill send buffers in Z direction.
      CALL LOAD_LOCAL_BUFFERS_UZ ()

! Fetch the Z-data from top and bottom patches.
      CALL FETCH_BUFFERS_FILL_IN_UZ ()

! ----------------------------------------------------------------------
! Synchronize Fluxes and E-fields in the same sequence as of MHD variables
! written above.

      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()                                 
      CALL FETCH_BUFFERS_FILL_IN_Flux_X () 

      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()                                 
      CALL FETCH_BUFFERS_FILL_IN_Flux_Y ()                               
                                                                         
      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()                                 
      CALL FETCH_BUFFERS_FILL_IN_Flux_Y_TOP ()                           

      CALL LOAD_LOCAL_BUFFERS_Flux_Z ()                                  
      CALL FETCH_BUFFERS_FILL_IN_Flux_Z () 

! ----------------------------------------------------------------------
! Assert the ghost zones and faces have received the correct values
! for the case of "iproblem = 0".

      IF ( me == iproc ) THEN
        WRITE ( 6, *)"----------------------------------------------"
        WRITE ( 6, *)"After MPI messaging: "
        WRITE ( 6, *)" "
      END IF

      CALL TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc,     
     1                 time, dtcur, rmin, rmax )

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! All the ghost zones and faces are filled. Now free up all the handles.
      CALL WINFREE_MESSAGING ()

! ----------------------------------------------------------------------

      is_bottom_or_top ( :, :) = 0.0

      DO ii = 1, npatch_per_proc

        iz_patch = ixyz_patch_amr ( 3, ii)

        IF ( iz_patch == 1)           is_bottom_or_top ( 1, ii) = 1
        IF ( iz_patch == nr_patches)  is_bottom_or_top ( 2, ii) = 1

      END DO

      filenum = 1;  time = 0.0

      CALL  WRITEDATA ( npatch_per_proc, 
     1                  iactive_FromProc ( :, 2),
     1                  nr_zones, ntp_zones, NGhost,
     1                  is_bottom_or_top,
     1                  itri_patch_amr,
     1                  vert_coords_amr,
     1                  n_cc_components, 
     1                  SILO_VARLEN, silo_allvars, 
     1                  U_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                  filenum, time
     1                )

! ----------------------------------------------------------------------

      CALL MPI_FINALIZE (ierr)

! ----------------------------------------------------------------------

      END PROGRAM TEST_RM_CS_MODULE

#endif
