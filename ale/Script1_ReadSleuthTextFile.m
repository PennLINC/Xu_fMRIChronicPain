clc, clear

allROI = spm_select(Inf,'\.txt$','Select Sleuth Output');


for rois = 1:size(allROI,1)
    
    inputROI = strrep(allROI(rois,:),' ','');
    clear Experiments XYZmm AllTasks AnzV cntExp Tasks XYZ
    
    TEMPLATE = spm_vol('./MaskenEtc/Grey10.nii');
    fid    = fopen(deblank(inputROI));
    study  = spm_str_manip(inputROI,'rt');
    
    Experiments = struct('Author',{});
    cntExp      = 0;
    XYZ         = [];
    
    line = fgetl(fid); introline = 1;
    if numel(deblank(line))==0; line = fgetl(fid);  end
    
    if numel(strfind(lower(line),'reference=mni'))==0;
        introline = 0;
        choice = questdlg('Are coordinates in MNI space?', 'Could not establish space!', ...
            'Yes','No','No');
        switch choice
            case 'No'
                error('Only works for MNI space data!')
                return
        end
    end
    
    
    while feof(fid) == 0
        if introline==1; line = fgetl(fid); else introline = 1; end
        fullline = line;
        if numel(deblank(line))>0
            if numel(strfind(line,'Reference'))==0
                if strcmp(line(1:2),'//')
                    cntExp = cntExp+1;
                    Experiments(cntExp).Author = fullline(3:end);
                    
                    line = fgetl(fid); fullline = line;
                    if numel(deblank(line))==0; line = fgetl(fid); fullline = line; end
                    Experiments(cntExp).Subjects = str2double(line(strfind(lower(line),'=')+1:end));
                    if isnan(Experiments(cntExp).Subjects)
                        line = fgetl(fid); fullline = line;
                        if numel(deblank(line))==0; line = fgetl(fid); fullline = line; end
                        Experiments(cntExp).Subjects = str2double(line(strfind(lower(line),'=')+1:end));
                    end
                    Experiments(cntExp).UncertainTemplates  =  (  5.7/(2*sqrt(2/pi))  * sqrt(8*log(2))) ;   % Assuming 5.7 mm ED between templates
                    Experiments(cntExp).UncertainSubjects   =  (  11.6/(2*sqrt(2/pi))  * sqrt(8*log(2))) / sqrt(Experiments(cntExp).Subjects);   % Assuming 11.6 mm ED between matching points
                    Experiments(cntExp).Smoothing           = sqrt(Experiments(cntExp).UncertainSubjects.^2 + Experiments(cntExp).UncertainTemplates.^2);
                    
                    line = fgetl(fid); fullline = line; cn = 1;
                    if numel(deblank(line))==0 & feof(fid) == 0; line = fgetl(fid); fullline = line; end
                    if numel(deblank(line))==0 & feof(fid) == 0; line = fgetl(fid); fullline = line; end
                    while strcmp(line(1:2),'//')
                        if strcmp(line(3),' ')
                            Experiments(cntExp).Cond{cn} = deblank(lower(fullline(4:end)));
                        else
                            Experiments(cntExp).Cond{cn} = deblank(lower(fullline(3:end)));
                        end
                        line = fgetl(fid); fullline = line;
                    end
                    Experiments(cntExp).XYZmm = str2num(strrep(line,',',' '))';
                    
                else
                    Experiments(cntExp).XYZmm = [Experiments(cntExp).XYZmm str2num(strrep(line,',',' '))'];
                end
            end
        end
    end
    
    fclose(fid);
    
    
    AnzV     = numel(Experiments);
    xleer    = zeros(31,31,31);  xleer(16,16,16) = 1;
    
    fprintf(1,'%s\n',[study ' - preparing ALE'])
    for cntExp=1:AnzV
        if size(Experiments(cntExp).XYZmm,2)>0
            Experiments(cntExp).XYZ   = round(inv(TEMPLATE.mat) * [(Experiments(cntExp).XYZmm(1:3,:)); ones(1,size(Experiments(cntExp).XYZmm,2))]);
            Experiments(cntExp).XYZ(1,Experiments(cntExp).XYZ(1,:)>TEMPLATE.dim(1)) = TEMPLATE.dim(1);
            Experiments(cntExp).XYZ(2,Experiments(cntExp).XYZ(2,:)>TEMPLATE.dim(2)) = TEMPLATE.dim(2);
            Experiments(cntExp).XYZ(3,Experiments(cntExp).XYZ(3,:)>TEMPLATE.dim(3)) = TEMPLATE.dim(3);
            Experiments(cntExp).XYZ(Experiments(cntExp).XYZ<1) = 1;
            Experiments(cntExp).Kernel = single(my_MemSmooth64bit(xleer,Experiments(cntExp).Smoothing,struct('dim',[31 31 31],'mat',TEMPLATE.mat),zeros(31)));
            Experiments(cntExp).Peaks  = size(Experiments(cntExp).XYZ,2);
        else
            Experiments(cntExp).Peaks  = 0;
        end
    end
    
    Experiments = Experiments([Experiments.Peaks]>0);
    Experiments = Experiments(~isnan([Experiments.Subjects]));
    
    
    AllTasks = {}; cnt = 1;
    try
        A = {Experiments.Cond};
        for i=1:numel(A)
            for ii=1:numel(A{i})
                AllTasks{cnt,1} = A{i}{ii}; cnt = cnt+1;
            end
        end
        AllTasks = unique(AllTasks);
        
        
        Tasks = cell2struct(AllTasks,'Name',2);
        spm_progress_bar('Init',numel(Tasks),'Preparing Tasks','Keywords Complete');
        cnt = 1;
        for i = 1:numel(Tasks)
            finder = zeros(size(Experiments));
            for ii = 1:numel(Experiments)
                for iii=1:numel(Experiments(ii).Cond)
                    finder(ii) = finder(ii) + strcmp(Experiments(ii).Cond{iii},Tasks(i).Name)>0;
                end
            end
            if sum(finder>0)>0 && numel(strrep(Tasks(i).Name,' ',''))>0
                Tasks(i).Experiments = find(finder>0);
                Tasks(i).Available   = sum(finder>0);
                Tasks(i).Wer         = {Experiments(Tasks(i).Experiments).Author};          %who was the author
                Tasks(i).WieViele    = sum([Experiments(Tasks(i).Experiments).Subjects]);   %how many subjects
                cnt = cnt+1;
            else
                Tasks(i).Available   = 0;
            end
            
        end
        
        Tasks(end+1).Name        = 'ALL';
        
    catch
        Tasks(1).Name        = 'ALL';
    end
    
    Tasks(end).Experiments   = 1:numel(Experiments);
    Tasks(end).Available     = numel(Experiments);
    Tasks(end).Wer           = {Experiments(Tasks(end).Experiments).Author};
    Tasks(end).WieViele    = sum([Experiments(Tasks(end).Experiments).Subjects]);   %how many subjects
    
    
    
    [B I] = sort([Tasks.Available],'descend');
    Tasks = Tasks(I);
    
    Tasks = Tasks([Tasks.Available]>0);
    
    
    for i = 1:numel(Tasks)
        
        fprintf(1,'%4.0f%s\t%5.0f%s\n',Tasks(i).Available,' Experiments',Tasks(i).WieViele,[' Subjects for ' Tasks(i).Name]);
        
        spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear')
    
    
    fprintf(1,'\n%s\n\n\n\n',[int2str(numel(Experiments)) ' Experiments in total ! '])
    
    save(fullfile(pwd,'DataMatlab',[study '.mat']),'Experiments','Tasks')
    
end