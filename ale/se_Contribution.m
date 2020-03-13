function se_Contribution(Experiments,study,Tasks,s0)

label = study;

[~, ~, ~] = mkdir('./ALE/Contribution');



tags = {'TFCE','FWE','cFWE'};


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



for run = 1:3
    
    fid = fopen(fullfile(pwd,'ALE','Contribution',[label '_' tags{run} '.txt']),'wt');
    
    fprintf(fid,'\n%s\n\n',['Starting with ' label '!']);
    
    [m x] = unique({Experiments.Author});
    fprintf(fid,'%s\n\n\n',[study ': ' int2str(numel(s0)) ' experiments; ' int2str(sum([Experiments(x).Subjects])) ' unique subjects (average of ' num2str(mean([Experiments.Subjects]),'%3.1f') ' per experiment)']);
    
    
    
    fil = dir(fullfile(pwd,'ALE','Results',[label '_' tags{run} '*.nii']));
    
    if numel(fil)>0
        
        Vi = spm_vol(fullfile(pwd,'ALE','Results',fil(1).name));
        dat = spm_read_vols(Vi);
        
        
        ind = find(dat>0);
        
        if numel(ind)>0
            [X Y Z] = ind2sub(Vi.dim,ind);
            A = spm_clusters([X Y Z]');
            dat = accumarray([X Y Z],A,Vi.dim); % array with the clusters
            
            cl  = unique(dat(dat>0)); % significant clusters
            
            ALE = spm_read_vols(spm_vol(fullfile(pwd,'ALE','Volumes',[label '.nii'])));
            
            % for the number of clusters that exist
            for i=1:numel(cl)
                % find the index of the significat cluster
                ind     = find(dat==cl(i)); % number of voxels within the cluster
                [X Y Z] = ind2sub(Vi.dim,ind);
                XYZ = median((Vi.mat * [X Y Z ones(numel(X),1)]')',1); % get the center of that cluster
                
                if numel(ind)>5 % only do this for clusters greater than size 5
                    fprintf(fid,'\n\n%s\n',['Cluster ' int2str(i) ': ' int2str(numel(ind)) ' voxel [Center: ' int2str(XYZ(1)) '/' int2str(XYZ(2)) '/' int2str(XYZ(3)) ']'  ]);
                    
                    xsum = [];
                    
                            
                             
                    Ax = [];
                    for ii=1:numel(Experiments); % for each experiment
                        Ax(ii,:) = Experiments(ii).dat(ind); % save, in the array, the probability-value of the voxels for this cluster
                             % and for experiment ii
                    end
                    AxF = 1-prod(1-Ax); % calculates the probabilities by taking the complement of the products of pr-values
                             % for all the areas for which you would not get this voxel
                             % i.e., 1-Ax=1 if this voxel has a pr-value of 0
                             % product of all these values gives you the probability of not getting the areas
                             % 1-prod would give you probability of getting the areas
                    
                    if max(abs([AxF - ALE(ind)']))>eps
                        error('ALE does not match union of MA-maps. Aborting')
                    end
                    
                    for ii=1:numel(Experiments);
                                % AxR is where you leave the experiment in question out
                                % AxR can be broken down to:
                                % 1-Ax([1:numel(Experiments)]~=ii,:) 1 - which is wherever there is not experiment ii
                                % 1-prod(1-Ax([1:numel(Experiments)]~=ii,:),1) is 1 - ^
                                % probability of wherever there was significance in cluster (taking out experiment ii)
                                % paragraph of this from white paper:
                                % "This is done by computing the ratio of the summarized test-values of all voxels of a specific cluster with and
                                % without the experiment in question, thus estimating how much the summarized test-value of this cluster would
                                % decrease when removing the experiment in question."
                        AxR = 1-prod(1-Ax([1:numel(Experiments)]~=ii,:),1);
                        
                        wig = sum(Experiments(ii).dat(ind)); % get the summed probability for the number of voxels
                                % in xsum(ii,:)
                                % wig is the sum of the probabilities for the cluster in experiment ii
                                % 100*wig/numel(ind) scales it by the number of voxels within that cluster
                                % (1-(AxR./AxF)) calculcates the pr of getting the value
                                    % by taking the complement of not seeing significance when the experiment in q is not there
                                    % we then take (1-(AxR./AxF)) and multiply it by 100 to get percentage
                                    % and get the max through all the indices without experiment ii
                                    % essentially gives you the contribution of that experiment
                        xsum(ii,:) = [wig, 100*wig/numel(ind),  100*(1-nanmean(AxR./AxF)),nanmax(100*(1-(AxR./AxF)))];
                                % consider (100*(1-(AxR./AxF)) as how much AxR./AxF differs from the ratio of 1
                                            % if AxR/AxF equals 1, then that would signify no change
                                            % if AxR/AxF equals 0, then that means AxR completely contributes
                    end
                    
                    xsum(:,3) = xsum(:,3)/sum(xsum(:,3))*100; % this calculcates the contribution of each experiment by
                                % (dividing curr contribution) / (total contribution possible)
                    
                    for ii=1:numel(Experiments);
                        if xsum(ii,3)>.1 | xsum(ii,4)>5
                            stx = repmat(' ',1,max(cellfun('length',{Experiments.Author})));
                            stx(1:numel(Experiments(ii).Author)) = Experiments(ii).Author;
                                % stx is the name of the experiment
                            fprintf(fid,'%s\t%7.3f\t%7.3f\t%7.2f\t%7.2f\t%s\n',stx,xsum(ii,:),['(' int2str(Experiments(ii).Subjects) ')']  );
                        end
                        
                    end
                    fprintf(fid,'\n\n');
                    
                    for ii=1:numel(Tasks);
                        % initialize an array of zeros the size of the number of experiments with the tag corresponding to Tasks
                        use = zeros(1,max([Tasks.Experiments]));  use(Tasks(ii).Experiments) = 1; % fill this array with its experiments
                        stx = repmat(' ',1,max(cellfun('length',{Tasks.Name})));
                        stx(1:numel(Tasks(ii).Name)) = Tasks(ii).Name;
                        buf = sum(xsum(use(s0)>0,2)); if isnan(buf), buf=0; end
                                % essentially, the contribution % scales with the number of studies that contributed to the ALE value
                                % i.e., if most studies that contributed were thermal tag, that'd be the tag
                        fprintf(fid,'%s\t%7.3f\t%7.3f\t%7.2f\n',stx,sum(xsum(use(s0)>0,1)), buf,  sum(xsum(use(s0)>0,3))  );
                        
                    end
                    
                end
            end
            
            
            fprintf(fid,'%s\n\n\n',['Done with ' label '!']);
            status = fclose(fid);
            
        else
            fprintf(fid,'%s\n\n\n',['No significant clusters in ' label '!']);
            status = fclose(fid);
        end
        
    else
        fprintf(fid,'%s\n\n\n',['Could not find ' tags{run} ' results for ' label '!']);
        status = fclose(fid);
        
    end
    
    
end
