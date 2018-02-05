function [] = DaganAnalyseExperiments(baseDname,tifStacksDnames,params)

nexp = length(tifStacksDnames);

% analysis
for iexp = 1 : nexp
    DaganAnalyseExperiment([baseDname filesep tifStacksDnames{iexp} filesep],params);
end
end