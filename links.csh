#!/bin/csh

# Note extra '../' in the paths so that this has to be used in a subdirectory of em_quarter_ss
#
set rundir = ../../run

echo linking to some physics data files in ${rundir} directory

             ln -sf ${rundir}/ETAMPNEW_DATA . ;                     \
             ln -sf ${rundir}/ETAMPNEW_DATA.expanded_rain . ;       \
             ln -sf ${rundir}/RRTM_DATA . ;                         \
             ln -sf ${rundir}/RRTMG_LW_DATA . ;                     \
             ln -sf ${rundir}/RRTMG_SW_DATA . ;                     \
             ln -sf ${rundir}/CAM_ABS_DATA . ;                      \
             ln -sf ${rundir}/CAM_AEROPT_DATA . ;                   \
             ln -sf ${rundir}/CAMtr_volume_mixing_ratio.RCP4.5 . ;  \
             ln -sf ${rundir}/CAMtr_volume_mixing_ratio.RCP6   . ;  \
             ln -sf ${rundir}/CAMtr_volume_mixing_ratio.RCP8.5 CAMtr_volume_mixing_ratio ;   \
             ln -sf ${rundir}/CAMtr_volume_mixing_ratio.A1B    . ;  \
             ln -sf ${rundir}/CAMtr_volume_mixing_ratio.A2     . ;  \
             ln -sf ${rundir}/CLM_ALB_ICE_DFS_DATA . ;              \
             ln -sf ${rundir}/CLM_ALB_ICE_DRC_DATA . ;              \
             ln -sf ${rundir}/CLM_ASM_ICE_DFS_DATA . ;              \
             ln -sf ${rundir}/CLM_ASM_ICE_DRC_DATA . ;              \
             ln -sf ${rundir}/CLM_DRDSDT0_DATA . ;                  \
             ln -sf ${rundir}/CLM_EXT_ICE_DFS_DATA . ;              \
             ln -sf ${rundir}/CLM_EXT_ICE_DRC_DATA . ;              \
             ln -sf ${rundir}/CLM_KAPPA_DATA . ;                    \
             ln -sf ${rundir}/CLM_TAU_DATA . ;                      \
             ln -sf ${rundir}/ozone.formatted . ;                   \
             ln -sf ${rundir}/ozone_lat.formatted . ;               \
             ln -sf ${rundir}/ozone_plev.formatted . ;              \
             ln -sf ${rundir}/aerosol.formatted . ;                 \
             ln -sf ${rundir}/aerosol_lat.formatted . ;             \
             ln -sf ${rundir}/aerosol_lon.formatted . ;             \
             ln -sf ${rundir}/aerosol_plev.formatted . ;            \
             ln -sf ${rundir}/capacity.asc . ;            \
             ln -sf ${rundir}/coeff_p.asc . ;            \
             ln -sf ${rundir}/coeff_q.asc . ;            \
             ln -sf ${rundir}/constants.asc . ;            \
             ln -sf ${rundir}/masses.asc . ;            \
             ln -sf ${rundir}/termvels.asc . ;            \
             ln -sf ${rundir}/kernels.asc_s_0_03_0_9 . ;            \
             ln -sf ${rundir}/kernels_z.asc . ;            \
             ln -sf ${rundir}/bulkdens.asc_s_0_03_0_9 . ;            \
             ln -sf ${rundir}/bulkradii.asc_s_0_03_0_9 . ;            \
             ln -sf ${rundir}/CCN_ACTIVATE.BIN . ;                   \
             ln -sf ${rundir}/p3_lookup_table_1.dat-v2.8.2 . ;                   \
             ln -sf ${rundir}/p3_lookup_table_2.dat-v2.8.2 . 
             
ln -s ${rundir}/GENPARM.TBL . 
ln -s ${rundir}/LANDUSE.TBL . 
 ln -s ${rundir}/SOILPARM.TBL .
ln -s ${rundir}/URBPARM.TBL .
ln -s ${rundir}/VEGPARM.TBL .
ln -s ${rundir}/MPTABLE.TBL .
ln -s ${rundir}/tr49t67 
ln -s ${rundir}/tr49t85 .
ln -s ${rundir}/tr67t85 .

ln -s ../../../main/ideal.exe .
ln -s ../../../main/wrf.exe .

