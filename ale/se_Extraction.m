function se_Extraction(Experiments,study,Tasks,s0,Mask)

label = study;

[~, ~, ~] = mkdir('./ALE/Extraction');

fid = fopen(fullfile(pwd,'ALE','Extraction',[label '__in_' spm_str_manip(Mask.fname,'rt') '.txt']),'wt');

fprintf(fid,'\n%s\n\n',['Starting with ' label '!']);

[m x] = unique({Experiments(s0).Author});
fprintf(fid,'%s\n\n\n',[study ': ' int2str(numel(s0)) ' experiments; ' int2str(sum([Experiments(s0(x)).Subjects])) ' unique subjects (average of ' num2str(mean([Experiments(s0).Subjects]),'%3.1f') ' per experiment)']);




Experiments = Experiments(s0);



AnzV     = numel(Experiments);
TEMPLATE = spm_vol('./MaskenEtc/Grey10.nii');
prior    = spm_read_vols(TEMPLATE)>.1;


xleer      = zeros(TEMPLATE.dim+[30 30 30]);

for i=1:AnzV
    data   = xleer;
    for ii = 1:Experiments(i).Peaks
        data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30) = ...
            max(data(Experiments(i).XYZ(1,ii):Experiments(i).XYZ(1,ii)+30,Experiments(i).XYZ(2,ii):Experiments(i).XYZ(2,ii)+30,Experiments(i).XYZ(3,ii):Experiments(i).XYZ(3,ii)+30),Experiments(i).Kernel);
    end
    Experiments(i).dat    = data(16:end-15,16:end-15,16:end-15);
end






dat = spm_read_vols(Mask);


ind = find(dat>0);

if numel(ind)>0
    [X Y Z] = ind2sub(Mask.dim,ind);
    A = spm_clusters([X Y Z]');
    dat = accumarray([X Y Z],A,Mask.dim);
    
    cl  = unique(dat(dat>0));
    
    ALE = spm_read_vols(spm_vol(fullfile(pwd,'ALE','Volumes',[label '.nii'])));
    
    
    for i=1:numel(cl)
        
        ind     = find(dat==cl(i));
        [X Y Z] = ind2sub(Mask.dim,ind);
        XYZ = median((Mask.mat * [X Y Z ones(numel(X),1)]')');
        
        
        fprintf(fid,'\n\n%s\n',['Cluster ' int2str(i) ': ' int2str(numel(ind)) ' voxel [Center: ' int2str(XYZ(1)) '/' int2str(XYZ(2)) '/' int2str(XYZ(3)) ']'  ]);
        
        
        xsum = [];
        
        Ax = [];
        for ii=1:numel(Experiments);
            Ax(ii,:) = Experiments(ii).dat(ind);
        end
        AxF = 1-prod(1-Ax);
        
        for ii=1:numel(Experiments);
            AxR = 1-prod(1-Ax([1:numel(Experiments)]~=ii,:),1);
            
            wig = sum(Experiments(ii).dat(ind));
            xsum(ii,:) = [wig, 100*wig/numel(ind),  100*(1-nanmean(AxR./AxF)),nanmax(100*(1-(AxR./AxF)))];
        end
        
        xsum(:,3) = xsum(:,3)/sum(xsum(:,3))*100;
        
        for ii=1:numel(Experiments);
            stx = repmat(' ',1,max(cellfun('length',{Experiments.Author})));
            stx(1:numel(Experiments(ii).Author)) = Experiments(ii).Author;
            fprintf(fid,'%s\t%7.3f\t%7.3f\t%7.2f\t%7.2f\t%s\n',stx,xsum(ii,:),['(' int2str(Experiments(ii).Subjects) ')']  );
        end
        fprintf(fid,'\n\n');        
            

        
        for ii=1:numel(Tasks);
            use = zeros(1,max([Tasks.Experiments]));  use(Tasks(ii).Experiments) = 1;
            stx = repmat(' ',1,max(cellfun('length',{Tasks.Name})));
            stx(1:numel(Tasks(ii).Name)) = Tasks(ii).Name;
            buf = sum(xsum(use(s0)>0,2)); if isnan(buf), buf=0; end
            fprintf(fid,'%s\t%7.3f\t%7.3f\t%7.2f\n',stx,sum(xsum(use(s0)>0,1)), buf,  sum(xsum(use(s0)>0,3))  );
            
        end
        
    end
    
    
    fprintf(fid,'\n\n%s\n\n\n',['Done with ' label '!']);
    status = fclose(fid);
    
else
    fprintf(fid,'\n\n%s\n\n\n',['No clusters in ' label '!']);
    status = fclose(fid);
end

