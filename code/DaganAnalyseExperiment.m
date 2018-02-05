function [] = DaganAnalyseExperiment(tifStackDname, params)

close all;

for icell = 1 : 100
    for iv = 1 : 100
        fname = ['c' num2str(icell) 'v' num2str(iv) '.tif'];
        if exist([tifStackDname fname],'file')
            DaganRegister(tifStackDname, fname,params);
        end
    end
end
end