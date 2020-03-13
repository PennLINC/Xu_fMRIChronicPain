% Renders image using blue as the 2nd color
se_render_imageCol(fullfile(pwd,'ALE','Foci','fslCustom_aberrantFoci.nii'), spm_invNcdf(0,0,1),20,[1 0 0; 0 0 1])
print('-dpng',fullfile(pwd,'ALE','Foci','fslCustom_aberrantFoci.png'))
se_ImageCut(fullfile(pwd,'ALE','Foci','fslCustom_aberrantFoci.png'),'X');
delete(fullfile(pwd,'ALE','Foci','fslCustom_aberrantFoci.png'));