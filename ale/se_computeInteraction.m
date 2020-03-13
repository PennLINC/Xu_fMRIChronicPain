function se_computeInteraction(IAExperiments, names, weights, doPrint)


doPrint     = strcmp(doPrint,'+');
useParallel = se_getParallelInfo;

TEMPLATE = spm_vol('./MaskenEtc/Grey10.nii');

IAname = [int2str(weights(1)) names{1} '_' int2str(weights(2)) names{2} '_' int2str(weights(3)) names{3} '_' int2str(weights(4)) names{4}];

if numel(dir(fullfile(pwd,'ALE','Interactions',[ IAname  '_P95.nii'])))==0
    
    for i=1:4
        fx = dir(fullfile(pwd,'ALE','VolumesZ',[names{i} '.nii']));
        dat(:,:,:,i) = spm_read_vols(spm_vol(fullfile(pwd,'ALE','VolumesZ',fx(1).name)));
    end
    dat = max(dat,[],4);
    ind = find(dat>2.4);
    
    Anz     = sum(cellfun('length',IAExperiments));
    VP      = nan(Anz, numel(ind));
    weight  = nan(Anz,2);
    xleer   = zeros(TEMPLATE.dim+[30 30 30]);
    
    cnt = 1;
    for xi=1:4
        for i=1:numel(IAExperiments{xi})
            data   = xleer; Experiments1 = IAExperiments{xi};
            for ii = 1:Experiments1(i).Peaks
                data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30),Experiments1(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); VP(cnt,:) = data(ind);
            weight(cnt,:) = [weights(xi) xi];
            cnt = cnt+1;
        end
    end
    VP = 1-VP;
    
    benchmark = (weights(1) * (1-prod(VP(weight(:,2)==1,:)))) + ...
        (weights(2) * (1-prod(VP(weight(:,2)==2,:)))) + ...
        (weights(3) * (1-prod(VP(weight(:,2)==3,:)))) + ...
        (weights(4) * (1-prod(VP(weight(:,2)==4,:))));
    
    for i=1:4
        ALE(i,:) = 1-prod(VP(weight(:,2)==i,:));
    end
    
    
    
    repeats = 25000;
    
    null = nan(repeats,numel(benchmark));
    
    null(1,:) = benchmark;
    if useParallel == 1
        parfor repeat = 2:repeats
            e = randperm(Anz); now = weight(e,:);
            null(repeat,:) = (weights(1) * (1-prod(VP(now(:,2)==1,:)))) + ...
                (weights(2) * (1-prod(VP(now(:,2)==2,:)))) + ...
                (weights(3) * (1-prod(VP(now(:,2)==3,:)))) + ...
                (weights(4) * (1-prod(VP(now(:,2)==4,:))));
        end;
    else
        for repeat = 2:repeats
            e = randperm(Anz); now = weight(e,:);
            null(repeat,:) = (weights(1) * (1-prod(VP(now(:,2)==1,:)))) + ...
                (weights(2) * (1-prod(VP(now(:,2)==2,:)))) + ...
                (weights(3) * (1-prod(VP(now(:,2)==3,:)))) + ...
                (weights(4) * (1-prod(VP(now(:,2)==4,:))));
        end;
    end
    
    
    p = sum(null>=repmat(benchmark,repeats,1))/repeats;
    
    Z = spm_invNcdf(1-p,0,1);
    
    [tX tY tZ] = ind2sub(TEMPLATE.dim,ind);
    pXYZ = [tX tY tZ]';
    
    Z    = Z(p<=0.05);
    pXYZ = pXYZ(:,p<=0.05);
    
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Interactions'));
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Interactions','Images'));
    
    
    xVo        = TEMPLATE;
    xVo.dt     = [64 1];
    xVo        = rmfield(xVo,'pinfo');
    xVo.fname  = fullfile(pwd,'ALE','Interactions',[ IAname  '_P95.nii']);
    if numel(pXYZ)>0
        xVo        = spm_write_vol(xVo,accumarray(pXYZ(1:3,:)',Z,xVo.dim));
    else
        xVo = spm_write_vol(xVo,zeros(xVo.dim));
    end
    
    
end


if doPrint
    try
        se_render_imageCol(fullfile(pwd,'ALE','Interactions',[ IAname  '_P95.nii']),0.1,50,NaN)
        print('-dpng',fullfile(pwd,'ALE','Interactions','Images',[ IAname '_P95.png']))
        se_ImageCut(fullfile(pwd,'ALE','Interactions','Images',[ IAname '_P95.png']),'X');
        delete(fullfile(pwd,'ALE','Interactions','Images',[ IAname '_P95.png']));
    end
    
    dat = spm_read_vols(spm_vol(fullfile(pwd,'ALE','Interactions',[ IAname  '_P95.nii'])));
    ind = find(dat>0);
    
    VP      = nan(sum(cellfun('length',IAExperiments)), numel(ind));
    weight  = nan(sum(cellfun('length',IAExperiments)),2);
    xleer   = zeros(TEMPLATE.dim+[30 30 30]);
    
    cnt = 1;
    for xi=1:4
        for i=1:numel(IAExperiments{xi})
            data   = xleer; Experiments1 = IAExperiments{xi};
            for ii = 1:Experiments1(i).Peaks
                data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30),Experiments1(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); VP(cnt,:) = data(ind);
            weight(cnt,:) = [weights(xi) xi];
            cnt = cnt+1;
        end
    end
    VP = 1-VP;
    
    clear ALE;  for i=1:4; ALE(i,:) = 1-prod(VP(weight(:,2)==i,:)); end
    
    [tX tY tZ] = ind2sub(TEMPLATE.dim,ind);
    XYZ = [tX tY tZ]';
    Z = dat(ind);
    A = spm_clusters(XYZ);
    
    cnt = 1;
    for cl = 1:max(A)
        if sum(A==cl)>49
            xyz = XYZ(:,A==cl);
            fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
            spm_orthviews('Image', spm_vol('./MaskenEtc/MNI152.nii'), [0.0 0.3 1 .8]);
            spm_orthviews('addblobs',1,xyz,Z(A==cl)'-(.5*min(Z(A==cl)')),TEMPLATE.mat,NaN);
            mXYZ = median([TEMPLATE.mat * [xyz; ones(1,size(xyz,2)) ]]');
            spm_orthviews('reposition',mXYZ(1:3));
            axes('Position',[0 .05 1 .25])
            bar(sum(ALE(:,A==cl),2))
            set(gca,'XTickLabel',strrep(names,'_',' '))
            print('-dpng',fullfile(pwd,'ALE','Interactions','Images',[ IAname '__C' int2str(cnt) '.png']));
            cnt = cnt+1;
            
        end
    end
    
    
    
end
