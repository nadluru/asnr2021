#!/bin/bash
# Concatenation of AP PA 75 76 and regridding to acquired resolution in plane and denoising
source CondorFunctions.sh
subjid=$1
parallel echo ADCP_${subjid}_3T_DWI_dir{2}_{3}.nii.gz ::: 75 76 ::: AP PA | parallel --plus -j1 -n4 --colsep " " 'mrcat {1} {3} {2} {4} -axis 3 - | mrgrid - regrid -scale 0.5,0.5,1.0 - -datatype int16le -interp sinc -force | dwidenoise -noise {..}_noiseSigma.nii.gz -datatype float64 - {..}_denoised.nii.gz' :::: -

# DeGibbs.
mrdegibbs ${subjid}_denoised.nii.gz {.}_degibbs.mif -nthreads 8 -datatype int16le

# eddy+topup
export bvecs=AP_PA_DV26_flipx.bvecs
export bvals=AP_PA_DV26_flipx.bvals
export rotime=0.110088
dwifslpreproc {.}_degibbs.mif {.}_degibbs_fsl_preprocessed.mif -pe_dir AP -rpe_all -readout_time $rotime -eddyqc_all {.}_eddyqcdir -fslgrad $bvecs $bvals -scratch {.}_scratch -force -nocleanup -nthreads 12 ::: $(pwd)/*denoised.mif > FSLPreproc.dag

# b1 bias correction
parallel -j1 'export executable=$(which dwibiascorrect);export job={/.};export numCPUs=1;export RAM="4 Gb";export initialDir=$(pwd)/B1biascorrect;mkdir -p $initialDir;export args="-force ants {} {.}_b1bc.mif -nthreads 1";CondorEcho' ::: $(pwd)/*preprocessed.mif > B1BiasCorrect.dag

# rician bias correction
parallel -j12 --dry-run --plus mrcalc -force {} -finite {} 0 -if {..}_lowb.mif -force ::: *_noiseSigma.nii.gz
parallel -j12 --dry-run --rpl '{i} s:_denoised.*::' 'mrcalc -force {} 2 -pow {i}_noiseSigma_lowb.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if {.}_rc.mif -force' ::: *b1bc.mif

# brain masking
parallel -j12 --dry-run 'dwiextract {} -bzero - | mrmath - mean {.}_meanB0.nii.gz -axis 3' ::: *rc.mif
parallel -j1 --plus 'export executable=$(which bet);export job={/..};export RAM="4 Gb";export numCPUs=1;export initialDir=$(pwd)/BET;mkdir -p $initialDir;export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' ::: $(pwd)/*meanB0.nii.gz > BET.dag

parallel --dry-run mrgrid {} regrid -vox 1.25 {.}_finalRes.mif -interp sinc -datatype int16le -force ::: *rc.mif
parallel -j6 --bar --rpl '{i} s:_AP.*::' 'dwiextract {} -bzero - | mrmath - mean {i}_meanB0_finalRes.nii.gz -axis 3 -force' ::: *finalRes.mif
parallel -j1 --plus 'export executable=$(which bet);export job={/..};export initialDir=$(pwd)/BETB0FinalRes;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' ::: $(pwd)/*_meanB0_finalRes.nii.gz

parallel -j1 --plus --rpl '{i} s:_meanB0.*::' 'export executable=$(which epi_reg);export args="-v --epi={i}_meanB0_t1wRes_bet.nii.gz --t1={i}_T1w_acpc_dc_restore.nii.gz --t1brain={i}_T1w_acpc_dc_restore_brain.nii.gz --out={i}_epi2mprage";export job={/..};export initialDir=$(pwd)/epi_reg;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";CondorEcho' ::: $(pwd)/*_meanB0_t1wRes.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/epi_reg_August182020.dag
parallel --dry-run convert_xfm -omat {.}_inverse.mat -inverse {} ::: *epi2mprage.mat

parallel --dry-run --rpl '{i} s:_T1.*::' flirt -in {} -ref {i}_meanB0_t1wRes_bet.nii.gz -applyxfm -init {i}_epi2mprage_inverse.mat -out {i}_mprageInDWI_t1wRes.nii.gz ::: *brain.nii.gz


# region Module 6
parallel --plus --dry-run 'export executable=$(which 5ttgen);export args="fsl {} {..}_5TT.mif -premasked -force -nocrop";export job={/..};export numCPUs=1;export RAM="4 Gb";export initialDir=$(pwd)/5TT;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*_mprageInDWI_t1wRes.nii.gz
parallel --bar -j6  mrgrid {} regrid -vox 1.25 {.}_finalRes.mif -force ::: *_mprageInDWI_t1wRes_5TT.mif
# endregion

# region module 7
parallel --dry-run --rpl '{i} s:_AP.*::' dwi2response -mask {i}_meanB0_finalRes_bet_mask.nii.gz dhollander {} {i}_sfwm.txt {i}_gm.txt {i}_csf.txt -force ::: *finalRes.mif
# fod estimation
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which dwi2fod);export args="-mask {i}_meanB0_finalRes_bet_mask.nii.gz -lmax 10,0,0 msmt_csd {} {i}_sfwm.txt {i}_wmfod.mif {i}_gm.txt {i}_gmfod.mif {i}_csf.txt {i}_csffod.mif -force";export job={/.};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/dwi2fod;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*rc_finalRes.mif
# fod normalization
parallel -j6 --dry-run --rpl '{i} s:_AP.*::' mtnormalise {i}_wmfod.mif {i}_wmfod_norm.mif {i}_gmfod.mif {i}_gmfod_norm.mif {i}_csffod.mif {i}_csffod_norm.mif -mask {i}_meanB0_finalRes_bet_mask.nii.gz -force ::: *rc_finalRes.mif
# mtnormalise double check
1086
2084
2095
2106
2114
# brain masking seems good.
# ignoring gm normalize seems to fix it.
# b0 thresholding (>0) is the best solution.

# Decided to use a script instead of tckgen and tcksift2 separately due to space constraints etc.
#parallel -j1 --rpl '{i} s:_mprage.*::' 'export executable=$(which tckgen);export args="{i}_wmfod_norm.mif {i}_tractogram.tck -act {} -backtrack -crop_at_gmwmi -maxlength 250 -power 0.33 -select 15M -seed_dynamic {i}_wmfod_norm.mif";export job={/.};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/tckgen;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*5TT_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/tckgen_August192020.dag
#parallel -j1 --rpl '{i} s:_mprage.*::' 'export executable=$(which tcksift2);export args="{i}_tractogram.tck {i}_wmfod_norm.mif {i}_weights.txt -act {} -out_mu {i}_mu.txt -fd_scale_gm";export job={/.};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/tcksift2;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*5TT_finalRes.mif
# endregion

# region T1w organizing
cd /study/utaut2/T1WIAnalysisNA/ADCP_T1w_nagesh
comm -13 <(ls ADCP_*/*brain.nii.gz | sed 's:/.*::;s:.*_::' | sort) <(ls ADCP_* -d | sed 's:.*_::g;s:/::' | sort)
1067
1068
2031
2086
parallel --bar cp {} {//}_{/} ::: ADCP_*/*brain.nii.gz
parallel --bar cp {} {//}_{/} ::: ADCP_*/*brain_1mm.nii.gz
parallel --bar cp {} {//}_{/} ::: ADCP_*/*restore.nii.gz
parallel --bar cp {} {//}_{/} ::: ADCP_*/*restore_1mm.nii.gz

# 8/20/2020
ls ADCP_????/*brain.nii.gz | sed 's:*::' | parallel --bar cp {} ADCP_DTI_nagesh/{//}_{/}
ls ADCP_????/*brain_1mm.nii.gz | sed 's:*::' | parallel --bar cp {} ADCP_DTI_nagesh/{//}_{/}
ls ADCP_????/*restore_1mm.nii.gz | sed 's:*::' | parallel --bar cp {} ADCP_DTI_nagesh/{//}_{/}
ls ADCP_????/*restore.nii.gz | sed 's:*::' | parallel --bar cp {} ADCP_DTI_nagesh/{//}_{/}
# endregion

# region processing missed subjects
# 8/20/2020
missedSubjs=$(comm -23 <(comm -13 <(ls *meanB0_t1wRes.nii.gz | sed 's:_mean.*::;s:.*_::' | sort) <(ls *brain.nii.gz | sed 's:_T1.*::;s:.*_::' | sort)) <(comm -13 <(ls *rc.mif | sed 's:_AP.*::;s:.*_::' | sort) <(ls *brain.nii.gz | sed 's:_T1.*::;s:.*_::' | sort)))

parallel echo *{}*brain.nii.gz ::: $missedSubjs | parallel --dry-run -j5 --plus --rpl '{i} s:_T1.*::g' 'mrgrid {i}_AP_PA_NoZFI_regrid_denoised_degibbs_fsl_preprocessed_b1bc_rc.mif regrid -vox $(mrinfo {} -spacing | awk "{print \$1}") - | dwiextract - -bzero - | mrmath - mean {i}_meanB0_t1wRes.nii.gz -axis 3 -force' :::: -

missedSubjs="1022 1038 1078 2059 2066"
parallel echo $(pwd)/*{}*_meanB0_t1wRes.nii.gz ::: $missedSubjs | parallel -j1 --plus 'export executable=$(which bet);export job={/..};export initialDir=$(pwd)/BETB0T1wRes;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' :::: - > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/BETB0T1wResMissedSubjs_August202020.dag

missedSubjs="1022 1038 1078 2059 2066"
parallel echo $(pwd)/*{}*_meanB0_t1wRes.nii.gz ::: $missedSubjs | parallel -j1 --plus --rpl '{i} s:_meanB0.*::' 'export executable=$(which epi_reg);export args="-v --epi={i}_meanB0_t1wRes_bet.nii.gz --t1={i}_T1w_acpc_dc_restore.nii.gz --t1brain={i}_T1w_acpc_dc_restore_brain.nii.gz --out={i}_epi2mprage";export job={/..};export initialDir=$(pwd)/epi_reg;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";CondorEcho' :::: - > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/epi_reg_August202020.dag

missedSubjs="1022 1038 1078 2059 2066"
parallel echo *{}*epi2mprage.mat ::: $missedSubjs | parallel --dry-run convert_xfm -omat {.}_inverse.mat -inverse {}

missedSubjs="1022 1038 1078 2059 2066"
parallel echo *{}*brain.nii.gz ::: $missedSubjs | parallel --dry-run --rpl '{i} s:_T1.*::' flirt -in {} -ref {i}_meanB0_t1wRes_bet.nii.gz -applyxfm -init {i}_epi2mprage_inverse.mat -out {i}_mprageInDWI_t1wRes.nii.gz

missedSubjs="1022 1038 1078 2059 2066"
parallel echo $(pwd)/*{}*_mprageInDWI_t1wRes.nii.gz ::: $missedSubjs | parallel --plus -j1 'export executable=$(which 5ttgen);export args="fsl {} {..}_5TT.mif -premasked -force -nocrop";export job={/..};export numCPUs=1;export RAM="4 Gb";export initialDir=$(pwd)/5TT;mkdir -p $initialDir;CondorEcho' > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/5TT_August202020.dag

missedSubjs="1022 1038 1078 2059 2066"
parallel echo *{}*_mprageInDWI_t1wRes_5TT.mif ::: $missedSubjs | parallel --dry-run -j5  mrgrid {} regrid -vox 1.25 {.}_finalRes.mif -force

missedSubjs="1022 1038 1078 2059 2066"
parallel echo *{}*_mprageInDWI_t1wRes.nii.gz ::: $missedSubjs | sed 's:_mp.*::g' | parallel -j1 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/AlignIITAtlas.sh;export args={};export job=IITReg_{};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/IITRegCondorLogs;mkdir -p $initialDir;CondorEcho' > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/IITReg_August202020.dag
# endregion

# region align atlas
ls *mprage*t1wRes.nii.gz | sed 's:_mp.*::g' | parallel -j1 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/AlignIITAtlas.sh;export args={};export job=IITReg_{};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/IITRegCondorLogs;mkdir -p $initialDir;CondorEcho' > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/IITReg_August192020.dag

parallel --dry-run --plus mrgrid {} regrid -vox 1.25 {..}_finalRes.nii.gz -interp nearest ::: IITReg_*/*_IIT_RANLRegistered_regions_IIT_GM_Desikan_atlas_mrtrix3.nii.gz
# endregion 

# region organizing 2059 and 2117
mv ADCP_2117_noiseSigma.nii.gz ADCP_2117_AP_PA_NoZFI_regrid_noiseSigma.nii.gz
mv ADCP_2059_noiseSigma.nii.gz ADCP_2059_AP_PA_NoZFI_regrid_noiseSigma.nii.gz
# endregion


# 5TT should be done with the highest resolution T1w available
# IIT to T1 at 1mm since IIT is at 1 mm.
# DWI at 1.25 mm
# DWI regridding > 0 for now working on this after eddy+topup. For ECP we should try this even with the first regrid.

# region fod with b0 thresholding
# 8/21/2020
parallel --dry-run -j6 'mrgrid {} regrid -vox 1.25 -interp sinc -datatype int16le - | mrcalc - 0 -lt 0 - -if {.}_finalRes.mif -force' ::: *rc.mif
parallel -j6 --bar --rpl '{i} s:_AP.*::' 'dwiextract {} -bzero - | mrmath - mean {i}_meanB0_finalRes.nii.gz -axis 3 -force' ::: *rc_finalRes.mif
parallel -j1 --plus 'export executable=$(which bet);export job={/..};export initialDir=$(pwd)/BETB0FinalResV2;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' ::: $(pwd)/*_meanB0_finalRes.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/BETB01FinalRes_August222020.dag

# condor seems to be slow here..
# parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which dwi2response);export args="-mask {i}_meanB0_finalRes_bet_mask.nii.gz dhollander {} {i}_sfwm.txt {i}_gm.txt {i}_csf.txt -force";export job={/.};export initialDir=$(pwd)/Response;mkdir -p $initialDir;export numCPUs=1;export RAM="8 Gb";CondorEcho' ::: $(pwd)/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/Response_August222020.dag


parallel --bar -j12 --rpl '{i} s:_AP.*::' dwi2response -mask {i}_meanB0_finalRes_bet_mask.nii.gz dhollander {} {i}_sfwm.txt {i}_gm.txt {i}_csf.txt -force ::: *rc_finalRes.mif
# fod estimation
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which dwi2fod);export args="-mask {i}_meanB0_finalRes_bet_mask.nii.gz -lmax 10,0,0 msmt_csd {} {i}_sfwm.txt {i}_wmfod.mif {i}_gm.txt {i}_gmfod.mif {i}_csf.txt {i}_csffod.mif -force";export job={/.};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/dwi2fodV2;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/FOD_August222020.dag
# fod normalization
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which mtnormalise);export args="{i}_wmfod.mif {i}_wmfod_norm.mif {i}_gmfod.mif {i}_gmfod_norm.mif {i}_csffod.mif {i}_csffod_norm.mif -mask {i}_meanB0_finalRes_bet_mask.nii.gz -force";export numCPUs=1;export job={/.};export RAM="4 Gb";export initialDir=$(pwd)/mtnormaliseCondor;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/FODNormalize_August232020.dag

# no gm normalization for now.
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which mtnormalise);export args="{i}_wmfod.mif {i}_wmfod_norm.mif {i}_csffod.mif {i}_csffod_norm.mif -mask {i}_meanB0_finalRes_bet_mask.nii.gz -force";export numCPUs=1;export job={/.};export RAM="4 Gb";export initialDir=$(pwd)/mtnormaliseCondor;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/FODNormalizeNoGM_August232020.dag
# endregion

# region t1w to dwi at 1mm for atlas to t1 and also remember to threshold the regridded DWI to have positive values.
parallel --dry-run -j6 --plus --rpl '{i} s:_T1.*::g' 'mrgrid {i}_AP_PA_NoZFI_regrid_denoised_degibbs_fsl_preprocessed_b1bc_rc.mif regrid -vox $(mrinfo {} -spacing | awk "{print \$1}") - | mrcalc - 0 -lt 0 - -if - | dwiextract - -bzero - | mrmath - mean {i}_meanB0_t1wRes_1mm.nii.gz -axis 3 -force' ::: *brain_1mm.nii.gz

parallel -j1 --plus 'export executable=$(which bet);export job={/..};export initialDir=$(pwd)/BETB0T1wRes_1mm;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' ::: $(pwd)/*_meanB0_t1wRes_1mm.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/BETB01TwRes1mm_August212020.dag

parallel -j1 --plus --rpl '{i} s:_meanB0.*::' 'export executable=$(which epi_reg);export args="-v --epi={i}_meanB0_t1wRes_1mm_bet.nii.gz --t1={i}_T1w_acpc_dc_restore_1mm.nii.gz --t1brain={i}_T1w_acpc_dc_restore_brain_1mm.nii.gz --out={i}_epi2mprage_1mm";export job={/..};export initialDir=$(pwd)/epi_reg_1mm;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";CondorEcho' ::: $(pwd)/*_meanB0_t1wRes_1mm.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/epi_reg_1mm_August212020.dag

parallel --dry-run convert_xfm -omat {.}_inverse.mat -inverse {} ::: *epi2mprage_1mm.mat
parallel --dry-run --rpl '{i} s:_T1.*::' flirt -in {} -ref {i}_meanB0_t1wRes_1mm_bet.nii.gz -applyxfm -init {i}_epi2mprage_1mm_inverse.mat -out {i}_mprageInDWI_t1wRes_1mm.nii.gz ::: *brain_1mm.nii.gz


ls *mprageInDWI_t1wRes_1mm.nii.gz | sed 's:_mp.*::g' | parallel -j1 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/AlignIITAtlas.sh;export args={};export job=IITReg_{};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/IITRegCondorLogs1mm;mkdir -p $initialDir;CondorEcho' > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/IITReg1mm_August222020_8Gb.dag

parallel -j1 --plus 'export executable=$(which mrgrid);export args="{} regrid -vox 1.25 {..}_finalRes.nii.gz -interp nearest";export job={/..};export initialDir=$(pwd)/iitmrgrid;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";CondorEcho' ::: $(pwd)/IITReg_*/*_IIT_RANLRegistered_regions_IIT_GM_Desikan_atlas_mrtrix3.nii.gz

parallel -j1 --plus --rpl '{i} s:_IIT.*::;s:.*/::' 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/StructuralConnectome.sh;export args={i};export job=SC_{i};export numCPUs=2;export RAM="8 Gb";export initialDir=$(pwd)/SC;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/IITReg_*/*Desikan_*_finalRes.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/SC_August242020.dag
# endregion

# 8/25/2020 epi2mprage, atlas align all at 1.25 mm isotropic since mrgrid of the atlas align step revealed that I don't fully understand mrgrid in such settings. regrid for DWI is fine since I have seen it being used by experts like Thjis. 5ttgen might be okay too.

# region re-process 8/25/2020
# finalRes DWI
parallel --dry-run -j6 'mrgrid {} regrid -vox 1.25 -interp sinc -datatype int16le - | mrcalc - 0 -lt 0 - -if {.}_finalRes.mif -force' ::: *rc.mif
parallel -j6 --dry-run --rpl '{i} s:_AP.*::' 'dwiextract {} -bzero - | mrmath - mean {i}_meanB0_finalRes.nii.gz -axis 3 -force' ::: *rc_finalRes.mif
parallel -j1 --plus 'export executable=$(which bet);export job={/..};export initialDir=$(pwd)/BETB0FinalResV3;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";export args="{} {..}_bet -f 0.3 -R -m";CondorEcho' ::: $(pwd)/*_meanB0_finalRes.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/BETB01FinalRes_August252020.dag

# epi2mprage
parallel -j1 --plus --rpl '{i} s:_meanB0.*::' 'export executable=$(which epi_reg);export args="-v --epi={i}_meanB0_finalRes_bet.nii.gz --t1={i}_T1w_acpc_dc_restore_1mm.nii.gz --t1brain={i}_T1w_acpc_dc_restore_brain_1mm.nii.gz --out={i}_epi2mprage_finalRes";export job={/..};export initialDir=$(pwd)/epi_reg_finalRes;mkdir -p $initialDir;export numCPUs=1;export RAM="4 Gb";CondorEcho' ::: $(pwd)/*_meanB0_finalRes.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/epi_reg_finalRes_August252020.dag
parallel --dry-run convert_xfm -omat {.}_inverse.mat -inverse {} ::: *epi2mprage_finalRes.mat
parallel --dry-run --rpl '{i} s:_T1.*::' flirt -in {} -ref {i}_meanB0_finalRes_bet.nii.gz -applyxfm -init {i}_epi2mprage_finalRes_inverse.mat -out {i}_mprageInDWI_finalRes.nii.gz ::: *brain_1mm.nii.gz

# atlas align
ls *_mprageInDWI_finalRes.nii.gz | sed 's:_mp.*::g' | parallel -j1 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/AlignIITAtlas.sh;export args={};export job=IITReg_{};export numCPUs=1;export RAM="8 Gb";export initialDir=$(pwd)/IITRegCondorLogsFinalRes;mkdir -p $initialDir;CondorEcho' > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/IITReg1mm_August252020_8Gb.dag

# structural connectome
parallel -j1 --plus --rpl '{i} s:_IIT.*::;s:.*/::' 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/StructuralConnectome.sh;export args={i};export job=SC_{i};export numCPUs=2;export RAM="8 Gb";export initialDir=$(pwd)/SCV2;mkdir -p $initialDir;CondorEcho' ::: $(pwd)/IITReg_FinalRes_*/*Desikan_*_mrtrix3.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/SC_August262020.dag
# endregion

# region reading graphs
# https://stackoverflow.com/questions/12751497/how-to-convert-csv-file-containing-network-data-into-gml
library(igraph)
my.data <- read.delim(url("http://dl.dropbox.com/u/22681355/email.csv"), sep = '')
my.graph <- graph.data.frame(d = my.data, directed = FALSE)
write.graph(graph = my.graph, file = 'email.gml', format = 'gml')
# endregion


# region fod spatial normalization
# respone averaging (on bender)
cd /mounts/data/preprocessed/modalities/dti/adcp_adluru/ADCP_DTI/FullSet_July2020/ADCP_DTI_nagesh/WaismanProcessed/AfterProcessing148Complete/ADCP_DTI_nagesh

responsemean *sfwm.txt GroupAverageResponse_sfwm.txt
responsemean *csf.txt GroupAverageResponse_csf.txt
responsemean *gm.txt GroupAverageResponse_gm.txt

# because these were not saved for some reason (on bender)
parallel --dry-run -j6 'mrgrid {} regrid -vox 1.25 -interp sinc -datatype int16le - | mrcalc - 0 -lt 0 - -if {.}_finalRes.mif -force' ::: *rc.mif

# fod estimation (on Waisman condor)
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which dwi2fod);export args="-mask {i}_meanB0_finalRes_bet_mask.nii.gz -lmax 10,0,0 msmt_csd {} {//}/GroupAverageResponse_sfwm.txt {i}_Group_wmfod.mif {//}/GroupAverageResponse_gm.txt {i}_Group_gmfod.mif {//}/GroupAverageResponse_csf.txt {i}_Group_csffod.mif -force -nthreads 2";export job={/.};export numCPUs=1;export RAM="8 Gb";export initialDir=/scratch/adcp/dwi2fodPop;mkdir -p $initialDir;CondorEcho' ::: /scratch/adcp/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/GroupFOD_November1022020.dag

# no gm normalization for now (on Waisman Condor)
parallel -j1 --rpl '{i} s:_AP.*::' 'export executable=$(which mtnormalise);export args="{i}_Group_wmfod.mif {i}_Group_wmfod_norm.mif {i}_Group_csffod.mif {i}_Group_csffod_norm.mif -mask {i}_meanB0_finalRes_bet_mask.nii.gz -force -nthreads 2";export numCPUs=1;export job={/.};export RAM="4 Gb";export initialDir=/scratch/adcp/mtnormalisePop;mkdir -p $initialDir;CondorEcho' ::: /scratch/adcp/*rc_finalRes.mif > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/GroupFODNormalizeNoGM_November112020.dag

# 02/01/2021
/study/utaut2/T1WIAnalysisNA/Tmp/grouprespfodsnorm/*_Group_wmfod_norm.mif
/study/utaut2/T1WIAnalysisNA/Tmp/IITReg
/study/utaut2/T1WIAnalysisNA/Tmp/5TT
parallel -j1 --plus --rpl '{i} s:_IIT.*::;s:.*/::' 'export executable=/home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/StructuralConnectome.sh;export args={i};export job=SC_{i};export numCPUs=2;export RAM="8 Gb";export initialDir=$(pwd)/SCV2;mkdir -p $initialDir;CondorEcho' ::: /study/utaut2/T1WIAnalysisNA/Tmp/IITReg/*Desikan_*_mrtrix3.nii.gz > /home/adluru/ADRCWRAPCPCPBiostatsProjects/adcp/SC_February012021.dag