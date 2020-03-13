function useParallel = se_getParallelInfo

try
    try
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            useParallel = 0;
        else
            if poolobj.NumWorkers>0
                useParallel = 1;
            else
                useParallel = 0;
            end;
        end
    catch
        isOpen = matlabpool('size') > 1;
        if isOpen
            useParallel = 1;
        else
            useParallel = 0;
        end
    end
catch
    useParallel = 0;
end
