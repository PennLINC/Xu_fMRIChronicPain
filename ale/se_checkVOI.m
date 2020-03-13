function se_checkVOI(Experiments,study,VOI,doPrint)


doPrint = strcmp(doPrint,'+');

try
    load(fullfile(pwd,'ALE','ROInull',[ study '.mat']))
catch
    
    TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));
    AnzV     = numel(Experiments);
    
    
    
    for i=1:numel(VOI) 
        dat = spm_read_vols(VOI(i)); % read in the ROI images
        dat = find(dat>0); clear XYZ % get the ones where values are greater than 1
        [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VOI(i).dim,dat);
        XYZ = TEMPLATE.mat \ VOI(i).mat * [XYZ'; ones(1,size(XYZ,1))];
        XYZ = round([XYZ(1,:)+15; XYZ(2,:)+15; XYZ(3,:)+15]);
        ind{i} = sub2ind(TEMPLATE.dim+[30 30 30],XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
    end
    
    xleer      = zeros(TEMPLATE.dim+[30 30 30]);
    leerALE    = ones(TEMPLATE.dim+[30 30 30]);
    
    for i=1:AnzV
        data   = xleer;
        for ii = 1:Experiments(i).Peaks % get experimental peaks and smooth with kernel
            data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30) = ...
                max(data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30),Experiments(i).Kernel);
        end
        leerALE = leerALE.*(1-data); 
    end
    leerALE         = 1-leerALE; % areas with significant convergence
    
    for i=1:numel(VOI)
        benchmarkIntegral(i) = sum(leerALE(ind{i})); 
        benchmarkMax(i) = max(leerALE(ind{i}));
    end
    
    
    
    
    fprintf(1,'\n\n%s\n',[study ' - simulating noise for VOI analyses'])
    
    
    useParallel = se_getParallelInfo;
    
    
    toRepeat = 10000;
    xleer    = zeros(TEMPLATE.dim+[30 30 30]);
    
    tmp = load(fullfile(pwd,'MaskenEtc','permSpace5.mat'));
    allXYZ = tmp.allXYZ;
    anzXYZ = tmp.anzXYZ;
    indice0 = tmp.indice0;
    indices = tmp.indices;
    clear tmp
    
    Kernel = {Experiments.Kernel};
    Peaks  = [Experiments.Peaks];
    clear Experiments
    
    nowIndices = []; nowVOI = [];
    for i=1:numel(VOI)
        nowIndices = [nowIndices; ind{i}];
        nowVOI = [nowVOI; ones(numel(ind{i}),1)*i];
    end
    totalVoxel = numel(nowIndices);
    
    tic
    
    
    alldata = nan(totalVoxel,toRepeat);
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
                
                alldata(:,runde)  = 1-prod(1-Vx);
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
                
                alldata(:,runde)  = 1-prod(1-Vx);
            end
            
        end
        fprintf(1,'%s\n',[int2str(xi*500) ' finished in ' num2str(toc/60,'%3.1f') ' minutes']);
    end
    
    
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROIresults'));
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROInull'));
    
    save(fullfile(pwd,'ALE','ROInull',[ study '.mat']),'VOI','alldata','benchmarkIntegral','benchmarkMax','nowVOI')
end


if doPrint
    for i=1:numel(VOI)
        
        for ii=1:size(alldata,2)
            NN(ii) = sum(alldata(nowVOI==i,ii));
            NM(ii) = max(alldata(nowVOI==i,ii));
        end
        
        figure(99), subplot(1,2,1), cla; hist(NN,50); maxi = get(gca,'YLim'); maxi = maxi(2);
        line([benchmarkIntegral(i) benchmarkIntegral(i)],[0 maxi*.8],'Color','r')
        text('position',[benchmarkIntegral(i) maxi*.8],'string',{'Observed'; 'value'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN)),'descend'); cut = cut(ceil(numel(NN)*.05));
        line([cut cut],[0 maxi*.66],'Color',[0 .4 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.05'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN)),'descend'); cut = cut(ceil(numel(NN)*.01));
        line([cut cut],[0 maxi*.66],'Color',[0 .6 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.01'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN)),'descend'); cut = cut(ceil(numel(NN)*.001));
        line([cut cut],[0 maxi*.66],'Color',[0 .8 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.001'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        title({spm_str_manip(VOI(i).fname,'rt'), ['ALE Integral: p < ' num2str(sum(NN>=benchmarkIntegral(i))/numel(NN),'%5.4f')]},'Interpreter','none','FontSize',16)
        set(gca,'YTickLabel',get(gca,'YTick')/numel(NN)); ylabel(['Percentage of ' int2str(numel(NN)) ' random realizations'],'FontSize',14)
        xlabel('ALE integral in search mask','FontSize',14)
        
        
        figure(99), subplot(1,2,2), cla; hist(NM,25); maxi = get(gca,'YLim'); maxi = maxi(2);
        line([benchmarkMax(i) benchmarkMax(i)],[0 maxi*.8],'Color','r')
        text('position',[benchmarkMax(i) maxi*.8],'string',{'Observed'; 'value'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM)),'descend'); cut = cut(ceil(numel(NM)*.05));
        line([cut cut],[0 maxi*.66],'Color',[0 .4 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.05'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM)),'descend'); cut = cut(ceil(numel(NM)*.01));
        line([cut cut],[0 maxi*.66],'Color',[0 .6 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.01'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM)),'descend'); cut = cut(ceil(numel(NM)*.001));
        line([cut cut],[0 maxi*.66],'Color',[0 .8 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.001'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        title({spm_str_manip(VOI(i).fname,'rt'), ['Maximum ALE: p < ' num2str(sum(NM>=benchmarkMax(i))/numel(NM),'%5.4f')]},'Interpreter','none','FontSize',16)
        set(gca,'YTickLabel',get(gca,'YTick')/numel(NM)); ylabel(['Percentage of ' int2str(numel(NN)) ' random realizations'],'FontSize',14)
        xlabel('Maximum ALE in search mask','FontSize',14)
        
        set(gcf,'Position',get(0,'ScreenSize'))
        hgexport(gcf,fullfile(pwd,'ALE','ROIresults',[study '__in__' spm_str_manip(VOI(i).fname,'rt') '.ps']))
        try; system(['/usr/bin/pstopdf ' fullfile(pwd,'ALE','ROIresults',[study '__in__' spm_str_manip(VOI(i).fname,'rt') '.ps'])]); end
        delete(gcf);
        
    end
    
end