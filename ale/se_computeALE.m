function se_computeALE(Experiments,study,doPrint)

doPrint = strcmp(doPrint,'+');

TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));
prior    = spm_read_vols(TEMPLATE)>.1;


uc  = 0.001;
c = num2str(uc); c = c(3:end);



[~, ~, ~] = mkdir('./ALE/ALEvolumes');
[~, ~, ~] = mkdir('./ALE/NullDistributions');
[~, ~, ~] = mkdir('./ALE/ALEvolumesZ');
[~, ~, ~] = mkdir('./ALE/Results');
[~, ~, ~] = mkdir('./ALE/Images');
[~, ~, ~] = mkdir('./ALE/Images/Foci');
[~, ~, ~] = mkdir('./ALE/Images/ALE');
[~, ~, ~] = mkdir('./ALE/Foci');



AnzV     = numel(Experiments);

mB = 1;
for i=1:AnzV
    mB = mB*(1-max(Experiments(i).Kernel(:)));
end

bin = [0:.0001:(1-mB)+.001];
step  = 1/mean(diff(bin));


if exist(fullfile(pwd,'ALE','Foci',[study '.nii']),'file')==0
    fprintf(1,'%s\n',[study ' - illustrate Foci'])
    xleer      = zeros(TEMPLATE.dim);;
    
    for i=1:AnzV
        for ii = 1:Experiments(i).Peaks
            xleer(Experiments(i).XYZ(1,ii),Experiments(i).XYZ(2,ii),Experiments(i).XYZ(3,ii)) = ...
                xleer(Experiments(i).XYZ(1,ii),Experiments(i).XYZ(2,ii),Experiments(i).XYZ(3,ii)) +1;
        end
    end
    
    xleer(~prior) = NaN;
    
    
    Vo       = TEMPLATE;
    Vo.fname = fullfile(pwd,'ALE','Foci',[study '.nii']);
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo       = spm_write_vol(Vo,xleer);
    
end
if doPrint & numel(dir(fullfile(pwd,'ALE','Images','Foci',[study '_LTR.png'])))==0
    se_render_image(fullfile(pwd,'ALE','Foci',[study '.nii']),0,0)
    print('-dpng',fullfile(pwd,'ALE','Images','Foci',[study '.png']))
    se_ImageCut(fullfile(pwd,'ALE','Images','Foci',[study '.png']),'X');
    delete(fullfile(pwd,'ALE','Images','Foci',[study '.png']))
end




if exist(fullfile(pwd,'ALE','NullDistributions',[ study '.mat']),'file')==0
    
    fprintf(1,'%s\n',[study ' - computing ALE'])
    xleer      = zeros(TEMPLATE.dim+[30 30 30]);
    leerALE    = ones(TEMPLATE.dim);
    Hx         = zeros(AnzV,numel(bin));
    
    for i=1:AnzV
        data   = xleer;
        for ii = 1:Experiments(i).Peaks
            data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30) = ...
                max(data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30),Experiments(i).Kernel);
        end
        data    = data(16:end-15,16:end-15,16:end-15);
        Hx(i,:) = hist(data(prior),bin);
        leerALE = leerALE.*(1-data);
    end
    
    leerALE         = 1-leerALE;
    leerALE(~prior) = NaN;
end


if exist(fullfile(pwd,'ALE','ALEvolumes', [study '.nii']),'file')==0
    Vo       = TEMPLATE;
    Vo.fname = fullfile(pwd,'ALE','ALEvolumes', [study '.nii']);
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo       = spm_write_vol(Vo,leerALE);
end
if doPrint & numel(dir(fullfile(pwd,'ALE','Images','ALE',[study '_LTR.png'])))==0
    se_render_image(fullfile(pwd,'ALE','ALEvolumes', [study '.nii']),0,0,.5)
    print('-dpng',fullfile(pwd,'ALE','Images','ALE',[study '.png']))
    se_ImageCut(fullfile(pwd,'ALE','Images','ALE',[study '.png']),'X');
    delete(fullfile(pwd,'ALE','Images','ALE',[study '.png']))
end


        
if exist(fullfile(pwd,'ALE','NullDistributions',[ study '.mat']),'file')==0
    fprintf(1,'%s\n',[study ' - permutation-null PDF'])
    step  = 1/mean(diff(bin));
    ALE   = Hx(1,:);
    
    for iMap=2:AnzV
        V1 = ALE;
        V2 = Hx(iMap,:);
        
        da1 = find(V1);            da2 = find(V2);
        V1  = V1/sum(V1);          V2  = V2/sum(V2);
        
        ALE = zeros(1,numel(bin));
        for i2 = 1:numel(da2)
            p     = V2(da2(i2))*V1(da1);
            score = 1-(1-bin(da2(i2)))*(1-bin(da1));
            wohin = round(score*step)+1;
            ALE(wohin) = ALE(wohin)+p;
        end
    end
    
    lastUsed = find(ALE>0,1,'last');
    cNULL    = fliplr(cumsum(fliplr(ALE(1:lastUsed))));
    
    save(fullfile(pwd,'ALE','NullDistributions',[ study '.mat']),'ALE','cNULL','lastUsed')
end

clear ALE V1 V2 ALEdist



if exist(fullfile(pwd,'ALE','ALEvolumesZ',[study '.nii']),'file')==0
    
    fprintf(1,'%s\n',[study ' - computing p-values'])
    load(fullfile(pwd,'ALE','NullDistributions',[ study '.mat']))
    Vi  = spm_vol(fullfile(pwd,'ALE','ALEvolumes',[ study '.nii']));
    try
        ALE = leerALE;
    catch
        ALE = spm_read_vols(spm_vol(fullfile(pwd,'ALE','ALEvolumes', [study '.nii'])));
    end
    
    p   = ones(Vi.dim);
    mx = numel(cNULL);
    for z=1:Vi.dim(3)
        [x y] = find(prior(:,:,z)>0 & ALE(:,:,z)>0);
        for i = 1:numel(x)
            indx = min(round(ALE(x(i),y(i),z)*step)+1,mx);
            p(x(i),y(i),z) = cNULL(indx);
        end
    end
    Z = spm_invNcdf(1-p,0,1);
    Z(p<eps) = spm_invNcdf(1-eps,0,1)+(ALE(p<eps)*2);
    
    p(~prior) = NaN;
    Z(~prior) = NaN;
    
    Vo       = Vi;
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo.fname = fullfile(pwd,'ALE','ALEvolumesZ',[study '.nii']);
    spm_write_vol(Vo,Z);
    
end




useParallel = se_getParallelInfo;


load(fullfile(pwd,'ALE','NullDistributions',[ study '.mat']))
toRepeat = 10000;

try
    load(fullfile(pwd,'ALE','NullDistributions',[study '_clustP.mat']));
    cut2; done = 1;
catch
    fprintf(1,'%s\n',[study ' - simulating noise']); done = 0;
end


if done == 0
    if useParallel == 0        
        
        try
            load(fullfile(pwd,'ALE','NullDistributions',[study '_clustP.mat']));
            StartRunde = runde;
        catch
            NM       = nan(1,toRepeat);
            NN       = nan(1,toRepeat);
            clusters = 1;
            StartRunde     = 1;
        end
        
        
        load(fullfile(pwd,'MaskenEtc','permSpace5.mat'))
        xleer      = zeros(TEMPLATE.dim+[30 30 30]);
        Vi         = single(nan(AnzV,numel(indices)));
        
        tic
        for runde = StartRunde+1:toRepeat
            
            for i=1:AnzV
                nowXYZ = allXYZ(:,ceil(rand(1,Experiments(i).Peaks)*anzXYZ));
                data   = xleer;
                for ii = 1:Experiments(i).Peaks
                    data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30) = ...
                        max(data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30),Experiments(i).Kernel);
                end
                Vi(i,:) = (data(indices));
            end
            
            data  = 1-prod(1-Vi);
            p     = cNULL(round(data*step)+1);
            
            NM(runde) = max(data(:));
            A = spm_clusters(double(allXYZ(1:3,p<uc)));
            
            if isempty(A); Q = 0;
            else;          Q = hist(A,1:max(A));
            end
            
            NN(runde) = max(Q);
            clusters = clusters+numel(Q);
            if rem(runde,250)==0
                fprintf(1,'%s\n',['Iteration ' int2str(runde) ' reached after ' num2str(toc/60,'%3.2f') ' min']);
                save(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat']),'NM','NN','clusters','runde')
            end
        end
        
        NM = sort(NM(~isnan(NM)),'descend');
        cut = NM(ceil(numel(NM)*.05));
        
        NN = sort(NN(~isnan(NN)),'descend');
        cut2 = NN(ceil(numel(NN)*.05));
        save(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat']),'NM','NN','cut','cut2')
        
        
    else
        
        try; load(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat'])); end
        
        if exist(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat']),'file')~=2 | exist('clusters','var')>0
                        
            NM       = nan(1,toRepeat);
            NN       = nan(1,toRepeat);
            
            load(fullfile(pwd,'MaskenEtc','permSpace5.mat'))
            nowIndices = indices;
            totalVoxel = numel(indices);
            AnzXYZ = anzXYZ;
            
            xleer      = single(zeros(TEMPLATE.dim+[30 30 30]));
            
            Kernel = {Experiments.Kernel};
            Peaks     = [Experiments.Peaks];
            clear Experiments
            
            tic
            for xi=1:ceil(toRepeat/500)
                parfor runde = ((xi-1)*500)+1:min((xi*500),toRepeat)
                    
                    Vx         = single(nan(AnzV,totalVoxel));
                    for i=1:AnzV
                        %                     stamp = single(se_MemSmooth(kernel,Smoothing(i)));
                        nowXYZ = allXYZ(:,ceil(rand(1,Peaks(i))*AnzXYZ));
                        data   = xleer;
                        for ii = 1:Peaks(i)
                            data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30) = ...
                                max(data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30),Kernel{i});
                        end
                        Vx(i,:) = (data(nowIndices));
                    end
                    
                    data  = 1-prod(1-Vx);
                    p     = cNULL(round(data*step)+1);
                    
                    NM(runde) = max(data(:));
                    A = spm_clusters(double(allXYZ(1:3,p<uc)));
                    
                    if isempty(A)
                        Q = 0;
                    else
                        Q = hist(A,1:max(A));
                    end
                    
                    NN(runde) = max(Q);
                end
                fprintf(1,'%s\n',[int2str(xi*500) ' finished in ' num2str(toc/60,'%3.1f') ' minutes'])
            end
            
            NM = sort(NM(~isnan(NM)),'descend');
            cut = NM(ceil(numel(NM)*.05));
            
            NN = sort(NN(~isnan(NN)),'descend');
            cut2 = NN(ceil(numel(NN)*.05));
            save(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat']),'NM','NN','cut','cut2')
            
        end
        
        
    end
    
end


if numel(dir(fullfile(pwd,'ALE','Results',[study '_cFWE05*.nii']))) ==0
    
    fprintf(1,'%s\n',[study ' - inference and printing'])
    
    load(fullfile(pwd,'ALE','NullDistributions',[ study '_clustP.mat']))
    
    
    
    VZ = spm_vol(fullfile(pwd,'ALE','ALEvolumesZ',[ study '.nii']));          Z = spm_read_vols(VZ);
    Z(prior == 0) = 0;      Z(Z < spm_invNcdf(1-uc,0,1)) =   0;
    
    if any(Z(:)); XYZ = []; zz = [];
        for p = 1:VZ.dim(3); d = Z(:,:,p);
            if any(any(d)); [i,j] = find(d>0); XYZ = [XYZ [i'; j'; p*ones(1,size(j,1))]]; zz = [zz d(find(d>0))']; end
        end;
        
        A = spm_clusters(XYZ);  Q = zeros(1,max(A));  for i = 1:max(A); Q(i) = sum(A == i); end
        writeXYZ = []; writeZ   = []; mini = Inf;
        
        if max(Q)>10
            [N I] = sort(Q,'descend');
            
            for i=1:numel(I)
                if N(i)>=cut2
                    mini = min(mini,N(i));
                    writeXYZ    = [writeXYZ XYZ(:,A==I(i))];
                    writeZ      = [writeZ zz(A==I(i))];
                end
                
            end
            
            if numel(writeXYZ)>0
                Vo          = VZ;
                c = num2str(uc); c = c(3:end);
                Vo.fname    = fullfile(pwd,'ALE','Results',[ study '_cFWE05_' c '_' int2str(mini) '.nii']);
                Vo = spm_create_vol(Vo);
                for p = 1:Vo.dim(3)
                    Q = find(writeXYZ(3,:) == p);
                    Vo = spm_write_plane(Vo, full(sparse(writeXYZ(1,Q),writeXYZ(2,Q),writeZ(Q),Vo(1).dim(1),Vo(1).dim(2))), p);
                end
                if doPrint & numel(dir(fullfile(pwd,'ALE','Images',[study '_cFWE05*.png'])))==0
                    se_render_image(Vo(1).fname,0,0)
                    print('-dpng',fullfile(pwd,'ALE','Images',[spm_str_manip(Vo(1).fname,'rt') '.png']))
                    se_ImageCut(fullfile(pwd,'ALE','Images',[spm_str_manip(Vo(1).fname,'rt') '.png']),'X');
                    delete(fullfile(pwd,'ALE','Images',[spm_str_manip(Vo(1).fname,'rt') '.png']))
                end
            else
                Vo          = VZ;
                Vo.fname    = fullfile(pwd,'ALE','Results',[ study '_cFWE05_' c  '_0.nii']);
                Vo = spm_write_vol(Vo,zeros(Vo.dim));
            end
            
        end
    end
end

if doPrint
    fil = dir(fullfile(pwd,'ALE','Results',[study '_cFWE05*.nii']));
    if numel(fil)==1 & numel(dir(fullfile(pwd,'ALE','Images',[study '_cFWE05*.png'])))==0
        if numel(strfind(fil(1).name,'_0.nii'))==0
            se_render_image(fullfile(pwd,'ALE','Results',fil(1).name),0,0)
            print('-dpng',fullfile(pwd,'ALE','Images',[spm_str_manip(fullfile(pwd,'ALE','Results',fil(1).name),'rt') '.png']))
            se_ImageCut(fullfile(pwd,'ALE','Images',[spm_str_manip(fullfile(pwd,'ALE','Results',fil(1).name),'rt') '.png']),'X');
            delete(fullfile(pwd,'ALE','Images',[spm_str_manip(fullfile(pwd,'ALE','Results',fil(1).name),'rt') '.png']))
        else
            A = zeros(100); A(:,:,2) = zeros(100); A(:,:,3) = zeros(100); figure(99), image(A)
            imwrite(A,fullfile(pwd,'ALE','Images',[study '_cFWE05_noSig.png']))
        end
    end
end


fprintf(1,'%s\n',[study ' - done!'])
