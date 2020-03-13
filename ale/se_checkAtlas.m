function se_checkAtlas(Experiments,study,VOI,doPrint)


doPrint = strcmp(doPrint,'+');

try
    load(fullfile(pwd,'ALE','ROInull',[ study '_' spm_str_manip(VOI.fname,'rt') '.mat']))
catch
    
    TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));
    AnzV     = numel(Experiments);
    
    VOI = VOI(1);
    
        data = spm_read_vols(VOI);
        dat = find(data>0); ids = data(dat); clear XYZ data
        
        in = unique(ids); id = nan(size(ids));
        for xi=1:numel(in)
            id(ids==in(xi)) = xi;
        end
        
        [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VOI.dim,dat);
        XYZ = TEMPLATE.mat \ VOI.mat * [XYZ'; ones(1,size(XYZ,1))];
        XYZ = round([XYZ(1,:)+15; XYZ(2,:)+15; XYZ(3,:)+15]);
        
        for i=1:max(id)
            XYZt = unique([XYZ(1,id==i)',XYZ(2,id==i)',XYZ(3,id==i)'],'rows');
            ind{i} = sub2ind(TEMPLATE.dim+[30 30 30],XYZt(:,1),XYZt(:,2),XYZt(:,3));
        end
        
    xleer      = zeros(TEMPLATE.dim+[30 30 30]);
    leerALE    = ones(TEMPLATE.dim+[30 30 30]);
    
    for i=1:AnzV
        data   = xleer;
        for ii = 1:Experiments(i).Peaks
            data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30) = ...
                max(data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30),Experiments(i).Kernel);
        end
        leerALE = leerALE.*(1-data);
    end
    leerALE         = 1-leerALE;
    
    for i=1:numel(ind)
        benchmarkIntegral(i) = sum(leerALE(ind{i}));
        benchmarkMax(i) = max(leerALE(ind{i}));
    end
    
    
    
    
    fprintf(1,'\n\n%s\n',[study ' - simulating noise for VOI analyses'])
    
    
    useParallel = se_getParallelInfo;
    
    
    toRepeat = 25000;
    xleer    = zeros(TEMPLATE.dim+[30 30 30]);
    
    load(fullfile(pwd,'MaskenEtc','permSpace5.mat'))
    
    Kernel = {Experiments.Kernel};
    Peaks  = [Experiments.Peaks];
    clear Experiments
    
    nowIndices = []; nowVOI = [];
    for i=1:numel(ind)
        nowIndices = [nowIndices; ind{i}];
        nowVOI = [nowVOI; ones(numel(ind{i}),1)*i];
    end
    totalVoxel = numel(nowIndices);
    
    tic
    
    
    NN = nan(numel(ind),toRepeat);
    NM = nan(numel(ind),toRepeat);
    
    AnzXYZ = anzXYZ;
    
    
    for xi=1:ceil(toRepeat/500)
        
        if useParallel == 0
            
            for runde = ((xi-1)*500)+1:min((xi*500),toRepeat)
                Vx         = single(nan(AnzV,totalVoxel));
                for i=1:AnzV
                    nowXYZ = allXYZ(:,ceil(rand(1,Peaks(i))*AnzXYZ));
                    data   = xleer;
                    for ii = 1:Peaks(i)
                        data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30) = ...
                            max(data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30),Kernel{i});
                    end
                    Vx(i,:) = (data(nowIndices));
                end
                
                alldata  = 1-prod(1-Vx);
                tmpN = nan(1,numel(ind));
                tmpM = nan(1,numel(ind));
                for i=1:numel(ind)
                    tmpN(i) = sum(alldata(nowVOI==i));
                    tmpM(i) = max(alldata(nowVOI==i));
                end
                NN(:,runde) = tmpN;
                NM(:,runde) = tmpM;
            end
            
        else
            
            parfor runde = ((xi-1)*500)+1:min((xi*500),toRepeat)
                Vx         = single(nan(AnzV,totalVoxel));
                for i=1:AnzV
                    nowXYZ = allXYZ(:,ceil(rand(1,Peaks(i))*AnzXYZ));
                    data   = xleer;
                    for ii = 1:Peaks(i)
                        data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30) = ...
                            max(data(nowXYZ(1,ii):nowXYZ(1,ii)+30,nowXYZ(2,ii):nowXYZ(2,ii)+30,nowXYZ(3,ii):nowXYZ(3,ii)+30),Kernel{i});
                    end
                    Vx(i,:) = (data(nowIndices));
                end
                
                alldata  = 1-prod(1-Vx);
                tmpN = nan(1,numel(ind));
                tmpM = nan(1,numel(ind));
                for i=1:numel(ind)
                    tmpN(i) = sum(alldata(nowVOI==i));
                    tmpM(i) = max(alldata(nowVOI==i));
                end
                NN(:,runde) = tmpN;
                NM(:,runde) = tmpM;
            end
            
        end
        fprintf(1,'%s\n',[int2str(xi*500) ' finished in ' num2str(toc/60,'%3.1f') ' minutes']);
    end
    
    
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROIresults'));
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROInull'));
    
    save(fullfile(pwd,'ALE','ROInull',[ study '_' spm_str_manip(VOI.fname,'rt') '.mat']),'VOI','benchmarkIntegral','benchmarkMax','NN','NM','in')
end


if doPrint
    for i=1:size(NN,1)
        p(1,i) = sum(NN(i,:)>=benchmarkIntegral(i))/numel(NN(i,:));
        p(2,i) = sum(NM(i,:)>=benchmarkMax(i))/numel(NM(i,:));
    end
    
    FDR(1) = max(spm_uc_FDR(0.05,1,'P',1,sort(p(1,:),'ascend')',[]),0.05/size(p,2));
    FDR(2) = max(spm_uc_FDR(0.05,1,'P',1,sort(p(2,:),'ascend')',[]),0.05/size(p,2));

        
    diary(fullfile(pwd,'ALE','ROIresults',[ study '_' spm_str_manip(VOI.fname,'rt') '.txt'])); diary on
    
    [B I] = sort(p(1,:),'ascend');
    fprintf(1,'%s\n',['Atlas: ' spm_str_manip(VOI.fname)]);
    
    fprintf(1,'%s\n','Analysis of ALE integral (* p<0.05 corrected)');
    for i=1:sum(B<0.05)
        if B(i)<FDR(1)
            fprintf(1,'%s\n',['* p<' num2str(B(i),'%4.3f') ' -> Cluster ' int2str(in(I(i)))]);
        else
            fprintf(1,'%s\n',['p<' num2str(B(i),'%4.3f') ' -> Cluster ' int2str(in(I(i)))]);
        end
    end
    
    
    [B I] = sort(p(2,:),'ascend');
    fprintf(1,'\n\n%s\n','Analysis of ALE peak (* p<0.05 corrected)');
    for i=1:sum(B<0.05)
        if B(i)<FDR(1)
            fprintf(1,'%s\n',['* p<' num2str(B(i),'%4.3f') ' -> Cluster ' int2str(I(i))]);
        else
            fprintf(1,'%s\n',['p<' num2str(B(i),'%4.3f') ' -> Cluster ' int2str(I(i))]);
        end
    end

    diary off
end