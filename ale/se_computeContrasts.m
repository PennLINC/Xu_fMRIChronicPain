function se_computeContrasts(Experiments1, Experiments2, study1, study2, col, doPrint)


doPrint     = strcmp(doPrint,'+');
useParallel = se_getParallelInfo;


u     = 0.05;


TEMPLATE = spm_vol('./MaskenEtc/Grey10.nii');

if numel(dir(fullfile(pwd,'ALE','Contrasts',[ study1 '--' study2  '_P95.nii'])))~=1
    
    
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Contrasts'));
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Contrasts','Images'));
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Conjunctions'));
    [~,~,~] =mkdir(fullfile(pwd,'ALE','Conjunctions','Images'));
    
    
    Anz1  = numel(Experiments1);
    Anz2  = numel(Experiments2);
    xleer = zeros(TEMPLATE.dim+[30 30 30]);
    repeats = 25000;
    
    
    fprintf(1,'%s','Positive contrast')
    fx = dir(fullfile(pwd,'ALE','Results',[study1 '*cFWE*.nii']));
    try
        dat = spm_read_vols(spm_vol(fullfile(pwd,'ALE','Results',fx(1).name)));
        % find wherever there's 'activation' in the first study
        ind = find(dat>0);
    catch
        ind = [];
    end
    
    if numel(ind)>0
        Z1  = dat(ind); clear dat
        VP  = nan(Anz1+Anz2, numel(ind));
        A   = zeros(size(ind))';
        
        for i=1:numel(Experiments1)
            data   = xleer;
            for ii = 1:Experiments1(i).Peaks
                data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30),Experiments1(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); 
            VP(i,:) = data(ind);
        end
        
        for i=1:numel(Experiments2)
            data   = xleer;
            for ii = 1:Experiments2(i).Peaks
                data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30),Experiments2(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); 
            VP(i+Anz1,:) = data(ind);
        end
        
        VP        = 1-VP;
        
        % where they differed
        benchmark = (1-prod(VP(1:Anz1,:))) - (1-prod(VP(Anz1+1:Anz1+Anz2,:)));
        % simulate a random difference
        fprintf(1,'%s','  randomising')
        if useParallel==0
            for repeat = 1:repeats
                e = randperm(Anz1+Anz2);
                draw = (1-prod(VP(e(1:Anz1),:)))-(1-prod(VP(e(Anz1+1:end),:)));
                A = A + double(draw>=benchmark);
                if rem(repeat,1000)==0; fprintf(1,'%s','.'); end
            end
        else
            parfor repeat = 1:repeats
                e = randperm(Anz1+Anz2);
                draw = (1-prod(VP(e(1:Anz1),:)))-(1-prod(VP(e(Anz1+1:end),:)));
                A = A + double(draw>=benchmark);
            end
        end
        
        fprintf(1,'\n')
        
        % get the z-scores
        A = spm_invNcdf(1-(A/repeats),0,1);
        A(isinf(A) & A>0) = spm_invNcdf(1-eps,0,1);
        Z1 = min(Z1',A); 
        ind1 = ind(Z1 > spm_invNcdf(1-u,0,1)); 
        Z1 = Z1(Z1>spm_invNcdf(1-u,0,1));
        
    else
        Z1 = []; ind1 = [];
    end
    
    % save Z1 ; save ind1
    
    
    fprintf(1,'%s','Negative contrast')
    fx = dir(fullfile(pwd,'ALE','Results',[study2 '*cFWE*.nii']));
    try
        dat = spm_read_vols(spm_vol(fullfile(pwd,'ALE','Results',fx(1).name)));
        ind = find(dat>0);
    catch
        ind = [];
    end
    
    if numel(ind)>0
        
        Z2  = dat(ind); clear dat
        VP  = nan(Anz1+Anz2, numel(ind));
        A   = zeros(size(ind))';
        
        for i=1:numel(Experiments1)
            data   = xleer;
            for ii = 1:Experiments1(i).Peaks
                data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30),Experiments1(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); VP(i,:) = data(ind);
        end
        
        for i=1:numel(Experiments2)
            data   = xleer;
            for ii = 1:Experiments2(i).Peaks
                data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30) = ...
                    max(data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30),Experiments2(i).Kernel);
            end
            data    = data(16:end-15,16:end-15,16:end-15); VP(i+Anz1,:) = data(ind);
        end
        VP        = 1-VP;
        benchmark = (1-prod(VP(Anz1+1:Anz1+Anz2,:))) - (1-prod(VP(1:Anz1,:)));
        
        fprintf(1,'%s','  randomising')
        if useParallel==0
            for repeat = 1:repeats
                e = randperm(Anz1+Anz2);
                draw = (1-prod(VP(e(Anz1+1:end),:))) - (1-prod(VP(e(1:Anz1),:)));
                A = A + double(draw>=benchmark);
                if rem(repeat,1000)==0; fprintf(1,'%s','.'); end
            end;
        else
            parfor repeat = 1:repeats
                e = randperm(Anz1+Anz2);
                draw = (1-prod(VP(e(Anz1+1:end),:))) - (1-prod(VP(e(1:Anz1),:)));
                A = A + double(draw>=benchmark);
            end
        end
        fprintf(1,'\n')
        
        A = spm_invNcdf(1-(A/repeats),0,1);
        A(isinf(A) & A>0) = spm_invNcdf(1-eps,0,1);
        Z2 = min(Z2',A);   ind2 = ind(Z2 > spm_invNcdf(1-u,0,1)); Z2 = Z2(Z2>spm_invNcdf(1-u,0,1));
        
    else
        Z2 = []; ind2 = [];
    end
    
    % save Z2, indices
    
    [tX tY tZ] = ind2sub(TEMPLATE.dim,[ind1; ind2]);
    pXYZ = [tX tY tZ]';
    
    
    fprintf(1,'\n%s\n',['Inference and printing'])
    xVo        = TEMPLATE;
    xVo.dt     = [64 1];
    xVo        = rmfield(xVo,'pinfo');
    xVo.fname  = fullfile(pwd,'ALE','Contrasts',[ study1 '--' study2  '_P95.nii']);
    if numel(pXYZ)>0
        xVo        = spm_write_vol(xVo,accumarray(pXYZ(1:3,:)',[Z1 -1*Z2],xVo.dim));
    else
        xVo = spm_write_vol(xVo,zeros(xVo.dim));
    end
    
end

if doPrint
try
    se_render_imageCol(fullfile(pwd,'ALE','Contrasts',[ study1 '--' study2  '_P95.nii']),spm_invNcdf(1-u,0,1),20,col)
    print('-dpng',fullfile(pwd,'ALE','Contrasts','Images',[ study1 '--' study2 '_P95.png']))
    se_ImageCut(fullfile(pwd,'ALE','Contrasts','Images',[ study1 '--' study2 '_P95.png']),'X');
    delete(fullfile(pwd,'ALE','Contrasts','Images',[ study1 '--' study2 '_P95.png']));
end
end

try
    fx1 = dir(fullfile(pwd,'ALE','Results',[study1 '_cFWE05*.nii']));
    fx2 = dir(fullfile(pwd,'ALE','Results',[study2 '_cFWE05*.nii']));
    dat     = min(spm_read_vols(spm_vol(fullfile(pwd,'ALE','Results',fx1.name))), spm_read_vols(spm_vol(fullfile(pwd,'ALE','Results',fx2.name))));
    ind     = find(dat);
    [X Y Z] = ind2sub(TEMPLATE.dim,ind);
    XYZ = [X Y Z]'; Z = dat(ind);
    
    A     = spm_clusters(XYZ);
    Q     = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= 50; Q = [Q j]; end
    end
    
    Z     = Z(Q);
    XYZ   = XYZ(:,Q);
    
    ind = sub2ind(TEMPLATE.dim,XYZ(1,:),XYZ(2,:),XYZ(3,:));
    dat = zeros(TEMPLATE.dim);
    dat(ind) = Z;
    
    xVo        = TEMPLATE;
    xVo.dt     = [64 1];
    xVo        = rmfield(xVo,'pinfo');
    xVo.fname  = fullfile(pwd,'ALE','Conjunctions',[ study1 '_AND_' study2  '_cFWE05.nii']);
    xVo        = spm_write_vol(xVo,dat.*(dat>spm_invNcdf(1-0.001,0,1)));
    
    if doPrint
    try
        se_render_imageCol(xVo.fname,spm_invNcdf(1-0.001,0,1),50,NaN)
        print('-dpng',fullfile(pwd,'ALE','Conjunctions','Images',[ study1 '_AND_' study2 '_cFWE05.png']))
        se_ImageCut(fullfile(pwd,'ALE','Conjunctions','Images',[ study1 '_AND_' study2 '_cFWE05.png']),'X');
        delete(fullfile(pwd,'ALE','Conjunctions','Images',[ study1 '_AND_' study2 '_cFWE05.png']))
    end
    end
end


