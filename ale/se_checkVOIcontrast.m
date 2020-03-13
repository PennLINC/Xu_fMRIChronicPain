function se_checkVOIcontrast(Experiments1,Experiments2, study,VOI,doPrint)


doPrint = strcmp(doPrint,'+');

try
    load(fullfile(pwd,'ALE','ROInull',[ study '.mat']))
catch
    TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));
    
    Anz1  = numel(Experiments1);
    Anz2  = numel(Experiments2);
    xleer = zeros(TEMPLATE.dim+[30 30 30]);
    
    for i=1:numel(VOI)
        dat = spm_read_vols(VOI(i));
        dat = find(dat>0); clear XYZ
        [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VOI(i).dim,dat);
        XYZ = TEMPLATE.mat \ VOI(i).mat * [XYZ'; ones(1,size(XYZ,1))];
        XYZ = round([XYZ(1,:)+15; XYZ(2,:)+15; XYZ(3,:)+15]);
        ind{i} = sub2ind(TEMPLATE.dim+[30 30 30],XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
        VP{i} = zeros(Anz1+Anz2,numel(ind{i}));
    end
    
    
    
    for i=1:numel(Experiments1)
        data   = xleer;
        for ii = 1:Experiments1(i).Peaks
            data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30) = ...
                max(data(Experiments1(i).XYZ(1,ii):Experiments1(i).XYZ(1,ii)+30,Experiments1(i).XYZ(2,ii):Experiments1(i).XYZ(2,ii)+30,Experiments1(i).XYZ(3,ii):Experiments1(i).XYZ(3,ii)+30),Experiments1(i).Kernel);
        end
        for ii=1:numel(ind)
            VP{ii}(i,:) = 1-data(ind{ii})';
        end
    end
    
    for i=1:numel(Experiments2)
        data   = xleer;
        for ii = 1:Experiments2(i).Peaks
            data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30) = ...
                max(data(Experiments2(i).XYZ(1,ii):Experiments2(i).XYZ(1,ii)+30,Experiments2(i).XYZ(2,ii):Experiments2(i).XYZ(2,ii)+30,Experiments2(i).XYZ(3,ii):Experiments2(i).XYZ(3,ii)+30),Experiments2(i).Kernel);
        end
        for ii=1:numel(ind)
            VP{ii}(i+Anz1,:) = 1-data(ind{ii})';
        end
    end
    
    
    for ii=1:numel(ind)
        voxelALEdiff = (1-prod(VP{ii}(1:Anz1,:))) - (1-prod(VP{ii}(Anz1+1:Anz1+Anz2,:)));
        benchmarkIntegral(ii) = sum(voxelALEdiff);
        benchmarkMax(ii)      = max(voxelALEdiff);
    end
    
    
    fprintf(1,'\n\n%s\n',[study ' - simulating noise for VOI analyses'])
    
    
    toRepeat = 25000;
    
    NN = nan(toRepeat,numel(ind));
    NM = nan(toRepeat,numel(ind));
    
    
    tic
    for xi=1:ceil(toRepeat/5000)
        for repeat = ((xi-1)*5000)+1:min((xi*5000),toRepeat)
            e = randperm(Anz1+Anz2);
            for ii=1:numel(ind)
                voxelALEdiff = (1-prod(VP{ii}(e(1:Anz1),:)))-(1-prod(VP{ii}(e(Anz1+1:end),:)));
                NN(repeat,ii) = sum(voxelALEdiff);
                NM(repeat,ii)      = max(voxelALEdiff);
            end
        end;
        fprintf(1,'%s\n',[int2str(xi*5000) ' finished in ' num2str(toc/60,'%3.1f') ' minutes']);
    end
    
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROIresults'));
    [~, ~, ~] = mkdir(fullfile(pwd,'ALE','ROInull'));
    
    save(fullfile(pwd,'ALE','ROInull',[ study '.mat']),'VOI','NN','NM','benchmarkIntegral','benchmarkMax')
end



if doPrint
    for i=1:numel(VOI)
        
        figure(99), subplot(1,2,1), cla; hist(NN(:,i),50); maxi = get(gca,'YLim'); maxi = maxi(2);
        line([benchmarkIntegral(i) benchmarkIntegral(i)],[0 maxi*.8],'Color','r')
        text('position',[benchmarkIntegral(i) maxi*.8],'string',{'Observed'; 'value'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN(:,i)),i),'descend'); cut = cut(ceil(numel(NN)*.05));
        line([cut cut],[0 maxi*.66],'Color',[0 .4 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.05'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN(:,i)),i),'descend'); cut = cut(ceil(numel(NN)*.01));
        line([cut cut],[0 maxi*.66],'Color',[0 .6 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.01'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NN(~isnan(NN(:,i)),i),'descend'); cut = cut(ceil(numel(NN)*.001));
        line([cut cut],[0 maxi*.66],'Color',[0 .8 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.001'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        title({spm_str_manip(VOI(i).fname,'rt'), ['ALE Integral: p < ' num2str(sum(NN(:,i)>=benchmarkIntegral(i))/size(NN,1),'%5.4f')]},'Interpreter','none','FontSize',16)
        set(gca,'YTickLabel',get(gca,'YTick')/numel(NN)); ylabel(['Percentage of ' int2str(size(NN,1)) ' random realizations'],'FontSize',14)
        xlabel('ALE integral in search mask','FontSize',14)
        
        
        figure(99), subplot(1,2,2), cla; hist(NM(:,i),25); maxi = get(gca,'YLim'); maxi = maxi(2);
        line([benchmarkMax(i) benchmarkMax(i)],[0 maxi*.8],'Color','r')
        text('position',[benchmarkMax(i) maxi*.8],'string',{'Observed'; 'value'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM(:,i)),i),'descend'); cut = cut(ceil(numel(NM)*.05));
        line([cut cut],[0 maxi*.66],'Color',[0 .4 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.05'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM(:,i)),i),'descend'); cut = cut(ceil(numel(NM)*.01));
        line([cut cut],[0 maxi*.66],'Color',[0 .6 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.01'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        cut = sort(NM(~isnan(NM(:,i)),i),'descend'); cut = cut(ceil(numel(NM)*.001));
        line([cut cut],[0 maxi*.66],'Color',[0 .8 0],'LineStyle','--'); text('position',[cut maxi*.66],'string',{'p < 0.001'},'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        title({spm_str_manip(VOI(i).fname,'rt'), ['Maximum ALE: p < ' num2str(sum(NM(:,i)>=benchmarkMax(i))/size(NM,1),'%5.4f')]},'Interpreter','none','FontSize',16)
        set(gca,'YTickLabel',get(gca,'YTick')/numel(NM)); ylabel(['Percentage of ' int2str(size(NN,1)) ' random realizations'],'FontSize',14)
        xlabel('Maximum ALE in search mask','FontSize',14)
        
        set(gcf,'Position',get(0,'ScreenSize'))
        hgexport(gcf,fullfile(pwd,'ALE','ROIresults',[study '__in__' spm_str_manip(VOI(i).fname,'rt') '.ps']))
        try; system(['/usr/bin/pstopdf ' fullfile(pwd,'ALE','ROIresults',[study '__in__' spm_str_manip(VOI(i).fname,'rt') '.ps'])]); end
        delete(gcf);
        
    end
end
