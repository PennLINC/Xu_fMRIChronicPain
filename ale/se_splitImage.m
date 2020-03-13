global defaults

mkdir(fullfile(pwd,'VOIs','VOIfiles'))


Vi = spm_vol(spm_select(1,'image','Select image'));

pm  = spm_input(['premultiply by'],'+0','r',1,1);

u  = spm_input(['Hight threshold (0 for none)'],'+0','r',0,1);
k  = spm_input(['Extend threshold (0 for none)'],'+0','r',0,1);

dat = spm_read_vols(Vi)*pm;
ind = find(dat>u);
[X Y Z] = ind2sub(Vi.dim,ind);
XYZ = [X Y Z]';
Z = dat(ind);
clear X Y dat ind;


A = spm_clusters(XYZ);

for i=1:numel(unique(A))
  
  if sum(A==i)>=k
    nXYZ        = XYZ(:,A==i);
    mXYZ   = mean(Vi.mat * [nXYZ; ones(1,size(nXYZ,2))],2); mXYZ   = mXYZ(1:3)';
    
    fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
    spm_orthviews('Image', spm_vol(fullfile(spm('dir'),'canonical','single_subj_T1.nii')), [0.0 0.22 1 .8]);
    spm_orthviews('addcolouredblobs',1,XYZ,Z,Vi(1).mat,[.5 0 0]);
    spm_orthviews('addcolouredblobs',1,nXYZ,Z(A==i),Vi(1).mat,[1 1 0]);
    spm_orthviews('reposition',mXYZ)
    
    ROIname = spm_input('VOI name [0 skips]',1,'s',[spm_str_manip(Vi.fname,'rt') '_']);
    
    if ~strcmp(ROIname,'0')
      Vo = Vi;
      Vo.fname = fullfile(pwd,'VOIs','VOIfiles',[ROIname '.nii']);
      
      Vo = spm_create_vol(Vo);
      for p = 1:Vi.dim(3)
        Q = find(nXYZ(3,:) == p);
        Vo = spm_write_plane(Vo, full(sparse(nXYZ(1,Q),nXYZ(2,Q),ones(sum(Q>0),1)',Vi(1).dim(1),Vi(1).dim(2))), p);
      end
    end
  end
end

