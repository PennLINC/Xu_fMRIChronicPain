function se_computeALEtfce(Experiments,study,doPrint)

doPrint = strcmp(doPrint,'+');


        
TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));
prior    = spm_read_vols(TEMPLATE)>.1;


uc  = 0.001;
c = num2str(uc); c = c(3:end);



[~, ~, ~] = mkdir('./ALE/Volumes');
[~, ~, ~] = mkdir('./ALE/NullDistributions');
[~, ~, ~] = mkdir('./ALE/VolumesZ');
[~, ~, ~] = mkdir('./ALE/VolumesTFCE');
[~, ~, ~] = mkdir('./ALE/Results');
[~, ~, ~] = mkdir('./ALE/Images');
[~, ~, ~] = mkdir('./ALE/Images/Foci');
[~, ~, ~] = mkdir('./ALE/Images/ALE');
[~, ~, ~] = mkdir('./ALE/Images/TFCE');
[~, ~, ~] = mkdir('./ALE/Foci');



AnzV     = numel(Experiments);
mB = 1;
for i=1:AnzV
    mB = mB*(1-max(Experiments(i).Kernel(:)));
end

bin = [0:.0001:(1-mB)+.001];
step  = 1/mean(diff(bin));


if exist(fullfile(pwd,'ALE/Foci',[study '.nii']),'file')==0
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
    Vo.fname = fullfile(pwd,'ALE/Foci',[study '.nii']);
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo       = spm_write_vol(Vo,xleer);
    
end
if doPrint & numel(dir(fullfile(pwd,'ALE/Images','Foci',[study '_LTR.png'])))==0
    se_render_image(fullfile(pwd,'ALE/Foci',[study '.nii']),0,0)
    print('-dpng',fullfile(pwd,'ALE/Images','Foci',[study '.png']))
    se_ImageCut(fullfile(pwd,'ALE/Images','Foci',[study '.png']),'X');
    delete(fullfile(pwd,'ALE/Images','Foci',[study '.png']))
end




if exist(fullfile(pwd,'ALE/NullDistributions',[ study '.mat']),'file')==0
    
    fprintf(1,'%s\n',[study ' - computing ALE'])
    xleer      = zeros(TEMPLATE.dim+[30 30 30]);
    leerALE    = ones(TEMPLATE.dim);
    Hx         = zeros(AnzV,numel(bin));
    
    for i=1:AnzV
        data   = xleer; % create a template
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


if exist(fullfile(pwd,'ALE/Volumes', [study '.nii']),'file')==0
    Vo       = TEMPLATE;
    Vo.fname = fullfile(pwd,'ALE/Volumes', [study '.nii']);
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo       = spm_write_vol(Vo,leerALE);
end
if doPrint & numel(dir(fullfile(pwd,'ALE/Images','ALE',[study '_LTR.png'])))==0
    se_render_image(fullfile(pwd,'ALE/Volumes', [study '.nii']),0,0,.5)
    print('-dpng',fullfile(pwd,'ALE/Images','ALE',[study '.png']))
    se_ImageCut(fullfile(pwd,'ALE/Images','ALE',[study '.png']),'X');
    delete(fullfile(pwd,'ALE/Images','ALE',[study '.png']))
end



        
if exist(fullfile(pwd,'ALE/NullDistributions',[ study '.mat']),'file')==0
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


if exist(fullfile(pwd,'ALE/VolumesZ',[study '.nii']),'file')==0
    
    fprintf(1,'%s\n',[study ' - computing p-values'])
    load(['./ALE/NullDistributions/' study '.mat'])
    Vi  = spm_vol(['./ALE/Volumes/' study '.nii']);
    ALE = leerALE;
    
    ALE(isnan(ALE)) = 0;
    p  = cNULL(round(ALE*step)+1);
    Z  = spm_invNcdf(1-p);
    Z(p<eps) = spm_invNcdf(1-eps,0,1)+(ALE(p<eps)*2);
    
    deltaT = max(Z(:))/100;
    tfce0 = tfceMex_pthread(Z,deltaT,0.6,2,0,0);
    tfce0(~prior) = NaN; tfce0(tfce0<0) = 0;
    
    Vo       =  rmfield(spm_vol(['./ALE/Volumes/' study '.nii']),'pinfo');;
    Vo.dt    = [64 1];
    Vo.fname = fullfile(pwd,'ALE/VolumesTFCE',[study '.nii']);
    spm_write_vol(Vo,tfce0);
    
    
    p(~prior) = NaN;
    Z(~prior) = NaN;
    
    Vo       = Vi;
    Vo.dt    = [64 1];
    Vo       = rmfield(Vo,'pinfo');
    Vo.fname = fullfile(pwd,'ALE/VolumesZ',[study '.nii']);
    spm_write_vol(Vo,Z);
    
end



if doPrint & numel(dir(fullfile(pwd,'ALE/Images','TFCE',[study '_LTR.png'])))==0
    se_render_image(fullfile(pwd,'ALE/VolumesTFCE', [study '.nii']),0,0,.5)
    print('-dpng',fullfile(pwd,'ALE/Images','TFCE',[study '.png']))
    se_ImageCut(fullfile(pwd,'ALE/Images','TFCE',[study '.png']),'X');
    delete(fullfile(pwd,'ALE/Images','TFCE',[study '.png']))
end








load(fullfile(pwd,'ALE/NullDistributions',[ study '.mat'])) 

toRepeat = 10000;

try
    load(fullfile(pwd,'ALE/NullDistributions',[study '_clustP.mat']));
    if xi == ceil(toRepeat/500)
        done=1;
    else
        startxi = xi+1;
        done = 0;
    end
catch
    fprintf(1,'%s\n',[study ' - simulating noise']); done = 0; startxi = 1;
end


if done == 0
        NM       = nan(1,toRepeat);
        NN       = nan(1,toRepeat);
        NT       = nan(1,toRepeat);
        
        
        load(fullfile(pwd,'MaskenEtc','permSpace5.mat'))
        totalVoxel = numel(indices);
        AnzXYZ = anzXYZ;
        cnull = cNULL;
            
        xleer       = single(zeros(TEMPLATE.dim+[30 30 30]));
        xleer2      = nan(TEMPLATE.dim+[30 30 30]);
        
        Z           = spm_read_vols(spm_vol(fullfile(pwd,'ALE/VolumesZ',[study '.nii'])));
        deltaT      = max(Z(:))/100;
        
        Kernel = {Experiments.Kernel};
        Peaks     = [Experiments.Peaks];
        clear Experiments
        
        tic
        for xi=startxi:ceil(toRepeat/1000)
            for runde = ((xi-1)*1000)+1:min((xi*1000),toRepeat)
                
                    Vx         = single(nan(AnzV,totalVoxel));
                    for i=1:AnzV
                        nowXYZ = allXYZ(:,ceil(rand(1,Peaks(i))*AnzXYZ));
                        data   = xleer;
                        for ii = 1:Peaks(i)
                            data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30) = ...
                                max(data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30),Kernel{i});
                        end
                        Vx(i,:) = (data(indices));
                    end
                
                data  = 1-prod(1-Vx);

                % Peak ALE threshold
                NM(runde) = max(data(:));

                % Cluster level threshold
                p           = cnull(round(data*step)+1);
                A           = spm_clusters(double(allXYZ(1:3,p<uc)));
                if isempty(A); Q = 0;
                else           Q = hist(A,1:max(A));
                end
                NN(runde) = max(Q);
                
                Z     = xleer2;   Z(indices)   = spm_invNcdf(1-p);
                
                Z(p<eps)    = spm_invNcdf(1-eps,0,1);
                Z = Z(16:end-15,16:end-15,16:end-15);
                
                % TFCE threshold
                tfce      = tfceMex_pthread(Z,deltaT,0.6,2,0,0);
                NT(runde) = max(tfce(:));


            end
            fprintf(1,'%s\n',[int2str(xi*1000) ' finished in ' num2str(toc/60,'%3.1f') ' minutes'])
            save(fullfile(pwd,'ALE/NullDistributions',[ study '_clustP.mat']),'NM','NN','NT','xi')
        end
        
end
    


if numel(dir(fullfile(pwd,'ALE/Results',[study '_cFWE05*.nii']))) ==0
    
    fprintf(1,'%s\n',[study ' - inference and printing'])
    load(fullfile(pwd,'ALE/NullDistributions',[ study '_clustP.mat']))
    
    NM       = sort(NM(~isnan(NM)),'descend');      cutMax   = NM(ceil(numel(NM)*.05));
    NN       = sort(NN(~isnan(NN)),'descend');      cutClust = NN(ceil(numel(NN)*.05));
    NT       = sort(NT(~isnan(NT)),'descend');      cutTFCE  = NT(ceil(numel(NT)*.05));

    VA = spm_vol(['./ALE/Volumes/' study '.nii']);          A = spm_read_vols(VA);          A(prior == 0) = 0;
    VZ = spm_vol(['./ALE/VolumesZ/' study '.nii']);          Z = spm_read_vols(VZ);          Z(prior == 0) = 0;
    VT = spm_vol(['./ALE/VolumesTFCE/' study '.nii']);          T = spm_read_vols(VT);       T(prior == 0) = 0;
    
    Vo          = rmfield(VZ,'pinfo');
    Vo.fname    = fullfile(pwd,'ALE/Results',[ study '_FWE05.nii']);
    if any(A(:)>cutMax)
        Vo = spm_write_vol(Vo,A.*(A>cutMax));
    else
        Vo = spm_write_vol(Vo,zeros(Vo.dim));
    end
               
    
    Vo          = rmfield(VZ,'pinfo');
    Vo.fname    = fullfile(pwd,'ALE/Results',[ study '_TFCE05.nii']);
    if any(T(:)>cutTFCE)
        Vo = spm_write_vol(Vo,T.*(T>=cutTFCE));
    else
        Vo = spm_write_vol(Vo,zeros(Vo.dim));
    end
    
    
    Z(Z < spm_invNcdf(1-uc,0,1)) =   0;
    
    clear XYZ zz; ind = find(Z>0); [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VZ.dim,ind); zz = Z(ind); XYZ = XYZ';
    A = spm_clusters(XYZ);  Q = [];  for i = 1:max(A); if sum(A==i)>= cutClust; Q = [Q find(A==i)]; end; end
    XYZ = XYZ(:,Q); zz = zz(Q);
    
    Vo          = rmfield(VZ,'pinfo');
    Vo.fname    = fullfile(pwd,'ALE/Results',[ study '_cFWE05.nii']);
    if numel(zz)>0
        Vo = spm_write_vol(Vo,accumarray(XYZ',zz,VZ.dim));
    else
        Vo = spm_write_vol(Vo,zeros(Vo.dim));
    end
end



    
    fprintf(1,'%s\n',[study ' - information on best(P)'])
    load(fullfile(pwd,'ALE/NullDistributions',[ study '_clustP.mat']))
    
    NM       = sort(NM(~isnan(NM)),'descend');      cutMax   = NM(ceil(numel(NM)*.05));
    NN       = sort(NN(~isnan(NN)),'descend');      cutClust = NN(ceil(numel(NN)*.05));
    NT       = sort(NT(~isnan(NT)),'descend');      cutTFCE  = NT(ceil(numel(NT)*.05));

    VA = spm_vol(['./ALE/Volumes/' study '.nii']);          A = spm_read_vols(VA);          A(prior == 0) = 0;
    VZ = spm_vol(['./ALE/VolumesZ/' study '.nii']);          Z = spm_read_vols(VZ);          Z(prior == 0) = 0;
    VT = spm_vol(['./ALE/VolumesTFCE/' study '.nii']);          T = spm_read_vols(VT);       T(prior == 0) = 0;
    
    fprintf(1,'%s\n',['Min p-value for vFWE: ' num2str(sum(NM>max(A(:)))/sum(~isnan(NM)),'%4.3f')])
    fprintf(1,'%s\n',['Min p-value for TFCE: ' num2str(sum(NT>max(T(:)))/sum(~isnan(NT)),'%4.3f')])
    
    clear ind XYZ A n
    Z(Z < spm_invNcdf(1-uc,0,1)) = 0; ind = find(Z>0); [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VZ.dim,ind); zz = Z(ind); XYZ = XYZ';
    A = spm_clusters(XYZ); n = hist(A,unique(A));
    fprintf(1,'%s\n',['Min p-value for cFWE: ' num2str(sum(NN>max(n(:)))/sum(~isnan(NN)),'%4.3f')])
    

    



if doPrint
    ext = {'cFWE','FWE','TFCE'};
    for i=1:3
        if numel(dir(fullfile(pwd,'ALE/Images',[study '_' ext{i} '05.png'])))==0
            dat = spm_read_vols(spm_vol(fullfile(pwd,'ALE/Results',[study '_' ext{i} '05.nii'])));
            if sum(dat(:)>0)>1
                se_render_image(fullfile(pwd,'ALE/Results',[study '_' ext{i} '05.nii']),0,0)
                print('-dpng',fullfile(pwd,'ALE/Images',[study '_' ext{i} '05.png']))
                se_ImageCut(fullfile(pwd,'ALE/Images',[study '_' ext{i} '05.png']),'X');
                delete(fullfile(pwd,'ALE/Images',[study '_' ext{i} '05.png']))
            else
                A = zeros(100); A(:,:,2) = zeros(100); A(:,:,3) = zeros(100); figure(99), image(A)
                imwrite(A,fullfile(pwd,'ALE/Images',[study '_' ext{i} '05__empty.png']))
            end
        end
    end
end


fprintf(1,'%s\n\n\n',[study ' - done!'])
