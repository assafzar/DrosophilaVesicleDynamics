%% Assess background and record raw intensities (as average of nPixwlMaxWinSize = 5 maximal values in a window)
function [] = DaganMetaBackgroundAnalysis(dname,fname)

nBackgroundTimePoints = 2;
nPixwlMaxWinSize = 5;

load([dname fname]);

nExps = length(metaData.groupsByExperiment);

for iExp = 1 : nExps
    expName = metaData.groupsByExperiment{iExp}.experiment{1};
    expDname = [dname expName];
    indsVes = metaData.groupsByExperiment{iExp}.inds;
    nVes = length(indsVes);
    
    dataTime = [];
    backgroundCh1 = [];
    backgroundCh2 = [];
            
    for iVes = 1 : nVes
        curVesInd = indsVes(iVes);
        
        curStime = metaData.tStart(curVesInd);
        
        assert(~isnan(curStime));
        
        % Get vesicle data (processed from kymograph)
        vesDname = [expDname filesep metaData.fnames{curVesInd}];
        vesFname = [vesDname filesep metaData.fnames{curVesInd} '_out.mat'];
        if ~exist(vesFname,'file')
            warning('%s does not exists',vesFname);
            continue;
        end
        load(vesFname);                      
        
        % Define intensities based on mean of nBackgroundTimePoints window pixels
        % at each time point
        dataCh1TimeRaw = nan(1,size(dataCh1,2));
        dataCh2TimeRaw = nan(1,size(dataCh2,2));
        for irow = 1 : size(dataCh1,1)
            dataCh1Sorted = sort(dataCh1(irow,:));
            dataCh1TimeRaw(irow) = mean(dataCh1Sorted(end-nPixwlMaxWinSize+1:end));
        end
        
        for irow = 1 : size(dataCh2,1)
            dataCh2Sorted = sort(dataCh2(irow,:));
            dataCh2TimeRaw(irow) = mean(dataCh2Sorted(end-nPixwlMaxWinSize+1:end));
        end
                
        curDataTime = curStime : curStime + length(dataCh1TimeRaw) - 1;

        % Update the time and background data        
        dataTime = [dataTime, curStime:(curStime+nBackgroundTimePoints-1)]; 
        backgroundCh1 = [backgroundCh1 dataCh1TimeRaw(1:nBackgroundTimePoints)];
        backgroundCh2 = [backgroundCh2 dataCh2TimeRaw(1:nBackgroundTimePoints)];

        assert(length(dataTime) == length(backgroundCh1));
        
        save([vesDname filesep metaData.fnames{curVesInd} '_vesRaw.mat'],'curDataTime','dataCh1TimeRaw','dataCh2TimeRaw','vesSize');
    end
    
    h = figure;
    hold on;
    plot(dataTime,backgroundCh1,'or','MarkerSize',8);
    plot(dataTime,backgroundCh2,'og','MarkerSize',8);
    xlabel('Frame #');
    ylabel('Background intensity');
    hold off;
    saveas(h,[expDname filesep expName '_backgroundAnalysis.tif']);
    
    
    meanBackgroundCh1 = mean(backgroundCh1); stdBackgroundCh1 = std(backgroundCh1);
    meanBackgroundCh2 = mean(backgroundCh2); stdBackgroundCh2 = std(backgroundCh2);
    [rBackgroundCh1,pvalBackgroundCh1] = corr(dataTime',backgroundCh1');
    [rBackgroundCh2,pvalBackgroundCh2] = corr(dataTime',backgroundCh2');
    N = length(dataTime);
    fprintf('\n\n **** background stats **** \n\n');
    fprintf(sprintf('%s (N = %d)\n',expName,N));
    fprintf(sprintf('Bleaching channel 1: pval = %g, rho = %g \n',pvalBackgroundCh1,rBackgroundCh1));
    fprintf(sprintf('Bleaching channel 2: pval = %g, rho = %g \n',pvalBackgroundCh2,rBackgroundCh2));    
    
    % Background model
    backgroundModelCh1 = nan; backgroundModelCh2 = nan;
    if pvalBackgroundCh1 < 0.05
        backgroundModelCh1 = polyfit(dataTime,backgroundCh1,2); 
    end
    if pvalBackgroundCh2 < 0.05
        backgroundModelCh2 = polyfit(dataTime,backgroundCh2,2); 
    end
    
    % Estimate background 
    if ~isnan(backgroundModelCh1)
        backgroundCh1Residual = backgroundCh1 - polyval(backgroundModelCh1,dataTime);
    else
        backgroundCh1Residual = backgroundCh1 - mean(backgroundCh1);
    end
    
    backgroundCh1ResidualStd = std(backgroundCh1Residual);
    
    if ~isnan(backgroundModelCh2)
        backgroundCh2Residual = backgroundCh2 - polyval(backgroundModelCh2,dataTime);
    else
        backgroundCh2Residual = backgroundCh2 - mean(backgroundCh2);
    end
    
    backgroundCh2ResidualStd = std(backgroundCh2Residual);
    
    save([expDname filesep expName '_backgroundAnalysis.mat'],...
        'expName','dataTime','backgroundCh1','backgroundCh2','N',...
        'rBackgroundCh1','pvalBackgroundCh1','rBackgroundCh2','pvalBackgroundCh2',...
        'meanBackgroundCh1','stdBackgroundCh1','meanBackgroundCh2','stdBackgroundCh2',...
        'backgroundModelCh1','backgroundModelCh2','backgroundCh1ResidualStd','backgroundCh2ResidualStd'...
        );
end

end