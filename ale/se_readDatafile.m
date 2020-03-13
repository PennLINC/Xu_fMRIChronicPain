function [Experiments, Tasks] = se_readDatafile(filename)

if strcmpi(filename(end-2:end),'mat')
    try;   load(fullfile(pwd,'DataMatlab',(filename)))
    catch;  error('Could not find data matfiles');
    end
else
    try;   load(fullfile(pwd,'DataMatlab','MetaAnalysis.mat'))
    catch; error('Could not find default data matfiles');
    end
end
