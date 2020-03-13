function toUse = se_compileStudies(Experiments,Tasks,B,Analysis)

    useIt = zeros(1,numel(Experiments));
    toUse = []; conJ = 0;   warOr = 0;

    for conditions = 5:size(B,2)

        if numel(B{Analysis,conditions})>0
            if B{Analysis,conditions}(1) == '+'
                useIt(Tasks(find(strcmpi({Tasks.Name},B{Analysis,conditions}(2:end)))).Experiments) = ...
                    useIt(Tasks(find(strcmpi({Tasks.Name},B{Analysis,conditions}(2:end)))).Experiments) +1;
                conJ = conJ+1;
                
            elseif B{Analysis,conditions}(1) == '?'
                useIt = double(useIt>warOr);
                conJ = 1;
                warOr = warOr+1;
                
            elseif B{Analysis,conditions}(1) == '-'
                useIt(Tasks(find(strcmpi({Tasks.Name},B{Analysis,conditions}(2:end)))).Experiments) = -Inf;
                
            elseif B{Analysis,conditions}(1) == '#'
                msk = spm_vol(B{Analysis,conditions}(2:end));
                msk.dat   = spm_read_vols(msk);
                
                nochDabei = 1:numel(useIt);
                useIt = double(useIt);

                for ii=1:numel(useIt)

                    XYZt = [];
                    XYZ0 = round(inv(msk.mat) * Experiments(ii).pMap.mat * [Experiments(ii).XYZ; ones(1,size(Experiments(ii).XYZ,2))]);
                    XYZ0 = XYZ0(1:3,:);
                    
                    for dx = -2:2
                        for dy=-2:2
                            for dz=-2:2
                                XYZt = [XYZt XYZ0+repmat([dx dy dz]',1,size(XYZ0,2))];
                            end
                        end
                    end

                    XYZt(XYZt<1) = 1;
                    XYZt(1,XYZt(1,:)>msk.dim(1)) = msk.dim(1);
                    XYZt(2,XYZt(2,:)>msk.dim(2)) = msk.dim(2);
                    XYZt(3,XYZt(3,:)>msk.dim(3)) = msk.dim(3);
                    indices   = sub2ind(msk.dim,XYZt(1,:),XYZt(2,:),XYZt(3,:));

                    if any(msk.dat(indices)>0)
                        useIt(nochDabei(ii)) = useIt(nochDabei(ii))+1;
                    end
                end
                conJ = conJ+1;

            elseif B{Analysis,conditions}(1) == '~'

                msk = spm_vol(B{Analysis,conditions}(2:end));
                nochDabei = find(useIt==conJ);

                useIt = double(useIt);
                
                for ii=1:numel(nochDabei)

                    XYZt = [];
                    for dx = -1:1
                        for dy=-1:1
                            for dz=-1:1
                                XYZt = [XYZt Experiments(nochDabei(ii)).XYZ+repmat([dx dy dz]',1,size(Experiments(nochDabei(ii)).XYZ,2))];
                            end
                        end
                    end

                    if any(spm_get_data(msk,Experiments(nochDabei(ii)).XYZ))
                        useIt(nochDabei(ii)) = -Inf;
                    end

                end
                
                
            end
        end
    end
    toUse = find(useIt==conJ);
