% Wrapper of Eickhoff's ALE scripts for computing ALE
% Takes in the .xlsx file of contrasts to be computed
% See documentation
% -------------------------------
% Commented by AX

function ale_estimateALE(analysisXLSX)

% aROIs = spm_select(Inf,'any','Select analyses speadsheets','',pwd,'.xlsx');
aROIs = analysisXLSX;
aROIs = mat2cell(aROIs,ones(1,size(aROIs,1)),size(aROIs,2));


for roi = 1:numel(aROIs)
    
    clear ROIs
    aROIs{roi} = (deblank(aROIs{roi}));
    
    [roiA roiB] = xlsread(aROIs{roi});
    
    
    for i=1:size(roiB,1)
        if numel(roiB{i,1})>0
            
            
            if strcmpi(roiB{i,1}(1),'M') % Main Effect
                [Experiments, Tasks] = se_readDatafile(roiB{i,3});
                s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                
                
                if numel(s0)>=12
                    fprintf(1,'%s\n',[roiB{i,2} ': ' int2str(numel(s0)) ' experiments; average of ' num2str(mean([Experiments(s0).Subjects]),'%3.1f') ' subjects per experiment'])
                    se_computeALEtfce(Experiments(s0),roiB{i,2},roiB{i,4})
                    se_Contribution(Experiments,roiB{i,2},Tasks,s0);
                    
                    % look for ROIs
                    if any(find(cellfun('length',strfind({roiB{i,3:end}},'$'))>0)) 
                        % read in ROI file
                        A = spm_vol(strrep({roiB{i,2+find(cellfun('length',strfind({roiB{i,3:end}},'$'))>0)}},'$',''));
                        for xi = 1:numel(A); aVOIs(xi) = A{xi}; end
                        % run function for ROI analysis
                        se_checkVOI(Experiments(s0),roiB{i,2},aVOIs,roiB{i,4})
                    end
                    
                    if any(find(cellfun('length',strfind({roiB{i,3:end}},'%'))>0))
                        
                        idx = 2+find(cellfun('length',strfind({roiB{i,3:end}},'%'))>0);
                        for xi=1:numel(idx)
                            aVOIs(xi) = spm_vol(fullfile(pwd,'VOIs',roiB{i,idx(xi)}(2:end)));
                        end
                        se_checkAtlas(Experiments(s0),roiB{i,2},aVOIs,roiB{i,4})
                    end
                else
                    fprintf(1,'%s\n\n\n',[roiB{i,2} ': only ' int2str(numel(s0)) ' experiments - not analyzed'])
                end
                
                
                
            elseif strcmpi(roiB{i,1}(1),'C') % Contrast
                
                if ~strcmpi(roiB{i,3},roiB{i+1,3}) % Different datafiles
                    
                    [Experiments, Tasks] = se_readDatafile(roiB{i,3});
                    s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                    Experiments1 = Experiments(s0); i = i+1;
                    
                    [Experiments, Tasks] = se_readDatafile(roiB{i,3});
                    s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                    Experiments2 = Experiments(s0);
                    
                    
                    if numel(Experiments1)>=12 && numel(Experiments2)>=12
                        
                        if numel(dir(fullfile(pwd,'ALE','Results',[ roiB{i-1,2} '_cFWE05*.nii'])))==0
                            se_computeALEtfce(Experiments1,roiB{i-1,2},roiB{i-1,4})
                        end
                        if numel(dir(fullfile(pwd,'ALE','Results',[ roiB{i,2} '_cFWE05*.nii'])))==0
                            se_computeALEtfce(Experiments2,roiB{i,2},roiB{i,4})
                        end
                        se_computeContrasts(Experiments1, Experiments2, roiB{i-1,2}, roiB{i,2},[1 0 0; 0 1 0],roiB{i-1,4})
                        
                        if any(find(cellfun('length',strfind({roiB{i-1,3:end}},'$'))>0))
                            A = spm_vol(strrep({roiB{i-1,2+find(cellfun('length',strfind({roiB{i-1,3:end}},'$'))>0)}},'$',''));
                            for xi = 1:numel(A);; aVOIs(xi) = A{xi}; end
                            se_checkVOIcontrast(Experiments1, Experiments2,[roiB{i-1,2} '--' roiB{i,2}],aVOIs,roiB{i-1,4})
                        end
                    end
                    
                    
                else % Same datafile
                    
                    [Experiments, Tasks] = se_readDatafile(roiB{i,3});
                    s1 = se_compileStudies(Experiments,Tasks,roiB,i);
                    s2 = se_compileStudies(Experiments,Tasks,roiB,i+1);
                    
                    if numel(s1)>=12 && numel(s2)>=12
                        
                        if numel(dir(fullfile(pwd,'ALE','Results',[ roiB{i,2} '_cFWE05*.nii'])))==0
                            se_computeALEtfce(Experiments(s1),roiB{i,2},roiB{i,4})
                        end
                        if numel(dir(fullfile(pwd,'ALE','Results',[ roiB{i+1,2} '_cFWE05*.nii'])))==0
                            se_computeALEtfce(Experiments(s2),roiB{i+1,2},roiB{i,4})
                        end
                        
                        q1 = ones(1,numel(s1)); q2 = ones(1,numel(s2));
                        for xi=1:numel(s1); for xii=1:numel(s2);; if s1(xi) == s2(xii); q1(xi)  = 0; q2(xii) = 0; end; end; end
                        useIt1 = s1(q1==1); useIt2 = s2(q2==1);
                        
                        i = i+1; se_computeContrasts(Experiments(useIt1), Experiments(useIt2), roiB{i-1,2}, roiB{i,2},[1 0 0; 0 1 0],roiB{i-1,4})
                        
                        if any(find(cellfun('length',strfind({roiB{i-1,3:end}},'$'))>0))
                            A = spm_vol(strrep({roiB{i-1,2+find(cellfun('length',strfind({roiB{i-1,3:end}},'$'))>0)}},'$',''));
                            for xi = 1:numel(A); aVOIs(xi) = A{xi}; end
                            se_checkVOIcontrast(Experiments(useIt1), Experiments(useIt2),[roiB{i-1,2} '--' roiB{i,2}],aVOIs,roiB{i-1,4})
                        end
                        
                    end
                    
                end
                
            elseif strcmpi(roiB{i,1}(1),'I') % Interaktion
                
                [Experiments, Tasks] = se_readDatafile(roiB{i,3}); s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                IAExperiments{1} = Experiments(s0); weight(1) = str2num(roiB{i,1}(2:3)); i = i+1;
                [Experiments, Tasks] = se_readDatafile(roiB{i,3}); s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                IAExperiments{2} = Experiments(s0); weight(2) = str2num(roiB{i,1}(2:3)); i = i+1;
                [Experiments, Tasks] = se_readDatafile(roiB{i,3}); s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                IAExperiments{3} = Experiments(s0); weight(3) = str2num(roiB{i,1}(2:3)); i = i+1;
                [Experiments, Tasks] = se_readDatafile(roiB{i,3}); s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                IAExperiments{4} = Experiments(s0); weight(4) = str2num(roiB{i,1}(2:3));
                
                for ii=1:4;
                    names{ii} = roiB{i-4+ii,2};
                    if numel(dir(fullfile(pwd,'ALE','Results',[ roiB{i-4+ii,2} '_cFWE05*.nii'])))==0; se_computeALEtfce(IAExperiments{ii},roiB{i-4+ii,2},roiB{i-4+ii,4}); end
                end
                
                se_computeInteraction(IAExperiments, names, weight,roiB{i-3,4})
                
                
                
            elseif strcmpi(roiB{i,1}(1),'E') % Extraction
                
                [Experiments, Tasks] = se_readDatafile(roiB{i,3});
                s0 = se_compileStudies(Experiments,Tasks,roiB,i);
                
                if any(find(cellfun('length',strfind({roiB{i,3:end}},'$'))>0))
                    A = spm_vol(strrep({roiB{i,2+find(cellfun('length',strfind({roiB{i,3:end}},'$'))>0)}},'$',''));
                    for xi = 1:numel(A); aVOIs(xi) = A{xi}; end
                    for xi = 1:numel(A); se_Extraction(Experiments,roiB{i,2},Tasks,s0,aVOIs(xi)); end
                    
                    
                end
                
            end
            
            
        else
            fprintf(1,'\n')
        end
    end
    
end
end 
