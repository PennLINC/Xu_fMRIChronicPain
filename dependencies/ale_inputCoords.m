% Read in coordinate spreadsheets and creates objects to use for analysis 
% Just a wrapper of Eickhoff's ALE scripts
% I/P: string representing the path to the coordinate .xls spreadsheet
% (see documentation)
% -------------------------------------
% commented by AX, 12/10/2018

function ale_inputCoords(coordXLS) 
dataname = coordXLS; 
[A B]    = xlsread(dataname);

%use header data from the Grey10 template
TEMPLATE = spm_vol(fullfile(pwd,'MaskenEtc','Grey10.nii'));

if size(B,1)>size(A,1)
    A = [nan(size(B,1)-size(A,1),size(A,2)); A];
end

% creates a struct object 'Experiments' that will take in the excel and
% change it to a MATLAB object with structure

Experiments = struct('Author',{});

Exlines = struct('Author',{});

cntExp  = 1; XYZ = [];
cntLine = 1;

%loops through the excel sheet pulling data into array Exlines
for loop = 2:size(A,1)
    
    if ~isempty(B{loop,1}) && sum(~isnan(A(loop,2:4)))==3
        %Read in author, subject & xyz data
        Exlines(cntLine).Author        = B{loop,1};
        Exlines(cntLine).Subjects      = A(loop,1);
        Exlines(cntLine).XYZ           = A(loop,2:4)';
        %Change the things these point to if the column values differ
        %(eg cov currently looks for covariates in columns 5-7)
        Exlines(cntLine).Space         = B{loop,6};
        %look for contrast and conditions tags?
        Exlines(cntLine).Cond = lower({B{loop,6+find(~cellfun('isempty',{B{loop,7:end}}))}});
        
        if cntLine>1
            cntExp = cntExp+1;
            if strcmp(Exlines(cntLine).Author,Exlines(cntLine-1).Author)
                if Exlines(cntLine).Subjects==Exlines(cntLine-1).Subjects
                    if numel(Exlines(cntLine).Cond) == numel(Exlines(cntLine-1).Cond)
                        if numel(cell2mat(Exlines(cntLine).Cond)) == numel(cell2mat(Exlines(cntLine-1).Cond))
                            if ~any(cell2mat(Exlines(cntLine).Cond) ~= cell2mat(Exlines(cntLine-1).Cond))
                                cntExp = cntExp-1;
                                try
                                    if isempty(B{loop-1,1})
                                        cntExp = cntExp+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
        Exlines(cntLine).ExpCount  = cntExp;
        cntLine = cntLine+1;
    end
    
end




Experiments = struct('Author',{});
%array of points corresponding to their experiment number
index       = [Exlines.ExpCount];

xleer    = zeros(31,31,31);  xleer(16,16,16) = 1;

 
for cntExp = 1:Exlines(end).ExpCount
    start = find(index==cntExp,1,'first');
    Experiments(cntExp).Author      = Exlines(start).Author;
    Experiments(cntExp).Subjects    = Exlines(start).Subjects;
    
    if Experiments(cntExp).Subjects == 0
        if any([Exlines(index==cntExp).Subjects])
            Experiments(cntExp).Subjects = mean([Exlines(index==cntExp).Subjects]);
        else
            Experiments(cntExp).Subjects = 1;
        end
    end
    Experiments(cntExp).Space      = Exlines(start).Space;
    Experiments(cntExp).Cond       = deblank(Exlines(start).Cond);
    
    Experiments(cntExp).XYZmm       = [Exlines(index==cntExp).XYZ];
    
    %Applies conversions if necessary, shifting to MNI space
    if isempty(Experiments(cntExp).Space), Experiments(cntExp).Space = 'MNI'; end
    if strcmpi(Experiments(cntExp).Space(1),'T')
        Experiments(cntExp).XYZmm = my_tal2icbm_spm([Experiments(cntExp).XYZmm nan(3,10)]);
        Experiments(cntExp).XYZmm = Experiments(cntExp).XYZmm(:,~isnan(sum(Experiments(cntExp).XYZmm)));
        Experiments(cntExp).XYZmm = round(Experiments(cntExp).XYZmm(1:3,:));
    end
    %Set values for gausian stamp depending on experiment variables(?)
    Experiments(cntExp).UncertainTemplates  =  (  5.7/(2*sqrt(2/pi))  * sqrt(8*log(2))) ;   % Assuming 5.7 mm ED between templates
    Experiments(cntExp).UncertainSubjects   =  (  11.6/(2*sqrt(2/pi))  * sqrt(8*log(2))) / sqrt(Experiments(cntExp).Subjects);   % Assuming 11.6 mm ED between matching points
    
    Experiments(cntExp).Smoothing = sqrt(Experiments(cntExp).UncertainSubjects.^2 + Experiments(cntExp).UncertainTemplates.^2);
    
    Experiments(cntExp).XYZ   = round(inv(TEMPLATE.mat) * [Experiments(cntExp).XYZmm; ones(1,size(Experiments(cntExp).XYZmm,2))]);
    Experiments(cntExp).XYZ(1,Experiments(cntExp).XYZ(1,:)>TEMPLATE.dim(1)) = TEMPLATE.dim(1);
    Experiments(cntExp).XYZ(2,Experiments(cntExp).XYZ(2,:)>TEMPLATE.dim(2)) = TEMPLATE.dim(2);
    Experiments(cntExp).XYZ(3,Experiments(cntExp).XYZ(3,:)>TEMPLATE.dim(3)) = TEMPLATE.dim(3);
    Experiments(cntExp).XYZ(Experiments(cntExp).XYZ<1) = 1;
    Experiments(cntExp).Kernel = single(my_MemSmooth64bit(xleer,Experiments(cntExp).Smoothing,struct('dim',[31 31 31],'mat',TEMPLATE.mat),zeros(31)));
    Experiments(cntExp).Peaks  = size(Experiments(cntExp).XYZ,2);
    
    
    cntExp = cntExp+1;
end


Experiments = Experiments(~isnan([Experiments.Subjects]));


AllTasks = {}; cnt = 1;
A = {Experiments.Cond};
for i=1:numel(A)
    for ii=1:numel(A{i})
        AllTasks{cnt,1} = A{i}{ii}; cnt = cnt+1;
    end
end
AllTasks = unique(AllTasks);

try
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
    Tasks(end).Experiments   = 1:numel(Experiments);
    Tasks(end).Available     = numel(Experiments);
    Tasks(end).Wer           = {'Alle'};
    Tasks(end).WieViele    = sum([Experiments(Tasks(end).Experiments).Subjects]);   %how many subjects
    
catch
    Tasks(1).Name        = 'ALL';
    Tasks(1).Experiments   = 1:numel(Experiments);
    Tasks(1).Available     = numel(Experiments);
    Tasks(1).Wer           = {'Alle'};
    Tasks(1).WieViele    = sum([Experiments(Tasks(end).Experiments).Subjects]);   %how many subjects
    
end



[B I] = sort([Tasks.Available],'descend');
Tasks = Tasks(I);

Tasks = Tasks([Tasks.Available]>0);


for i = 1:numel(Tasks)
    fprintf(1,'%4.0f%s\t%5.0f%s\n',Tasks(i).Available,' Experiments',Tasks(i).WieViele,[' Subjects for ' Tasks(i).Name]);
    spm_progress_bar('Set',i);
end


fprintf(1,'\n%s\n',[int2str(numel(Experiments)) ' Experiments in total ! '])



save(fullfile(pwd,'DataMatlab',[spm_str_manip(dataname,'rt') '.mat']),'Experiments','Tasks')
[spm_str_manip(dataname,'rt') '.mat']
end 