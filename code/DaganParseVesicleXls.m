function [] = DaganParseVesicleXls(dname,fname)

if nargin < 2
    dname = 'C:\Users\assafza\Google Drive\Research\PostDoc\Collaborations\Dagan\rhogapi_rockinh';
    fname = 'vesicles_timepoints_rhogapi_20170815.xlsx';
end

dbFname = [dname filesep fname];

[pathstr,name,ext] = fileparts(dbFname);

firstRow = 2;

metaData = getExperiments(dbFname,firstRow);

metaData.groupsByExperiment = groupByExperiments(metaData);

outFname = [pathstr '\' name '.mat'];

save(outFname,'metaData');
end

% function [exps] = getExperiments(fname,firstRow,lastRow)
function [exps] = getExperiments(fname,firstRow)
exps.experiments = {};
exps.cells = {};
exps.vesicles = {};
exps.fnames = {};
exps.tStart = nan;
exps.tEnd = nan;

[numbers strings misc] = xlsread(fname);
lastRow = size(strings,1);
exps.N = lastRow - firstRow + 1;

[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:H%d',firstRow,lastRow));  

exps.tStart = numsAll(firstRow-1 : lastRow-1,1);
exps.tEnd = numsAll(firstRow-1 : lastRow-1,2);

for line = firstRow : lastRow 
    i = line - firstRow + 1;
    curExp = txtAll{i,1};
    curCell = txtAll{i,3};
    curVes = txtAll{i,4};
        
    exps.experiments{i} = curExp;
    exps.cells{i} = curCell;
    exps.vesicles{i} = curVes;
    exps.fnames{i} = [curCell curVes];
end
end


% groupsByTreatments - cell array tha holds pairs <treatmentStr,[list of
% indices]>
function [out] = groupByExperiments(exps)
out = {};
for r = 1 : exps.N
    t = 1;
    nGroups = length(out);
    while t <= nGroups
        if strcmp(exps.experiments(r),out{t}.experiment) == 1
            out{t}.inds = [out{t}.inds,r];            
            break;
        else
            t = t + 1;
        end
    end
    if t > nGroups
        out{t}.experiment = exps.experiments(r);
        out{t}.inds = [r];
    end
end
end
