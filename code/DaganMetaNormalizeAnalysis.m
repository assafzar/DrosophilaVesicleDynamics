function [] = DaganMetaNormalizeAnalysis(dname,fname,params)

load([dname fname]);

if ~exist([dname filesep 'out'],'dir')
    mkdir([dname filesep 'out']);
end

nExps = length(metaData.groupsByExperiment);

smoothKernel = gausswin(params.smoothWinSize,params.smoothSigma);
smoothKernelNorm = smoothKernel / sum(smoothKernel); % Normalize.

allExpsVesData = cell(1,nExps);
allExpsStr = cell(1,nExps);
for iExp = 1 : nExps
    timePerFrame = params.timePerFrame(iExp);
    
    expName = metaData.groupsByExperiment{iExp}.experiment{1};
    expDname = [dname expName];
    indsVes = metaData.groupsByExperiment{iExp}.inds;
    nVes = length(indsVes);
    
    
    %% Linear model for background, if needed...
    load([expDname filesep expName '_backgroundAnalysis.mat']);
    %     'expName','dataTime','backgroundCh1','backgroundCh2','N',...
    %         'rBackgroundCh1','pvalBackgroundCh1','rBackgroundCh2','pvalBackgroundCh2',...
    %         'meanBackgroundCh1','stdBackgroundCh1','meanBackgroundCh2','stdBackgroundCh2',...
    %         'backgroundModelCh1','backgroundModelCh2','backgroundCh1ResidualStd','backgroundCh2ResidualStd'...
    
    
    
    expVesData = cell(1,nVes);
    
    for iVes = 1 : nVes
        curVesInd = indsVes(iVes);               
        % Get vesicle data (processed from kymograph)
        vesDname = [expDname filesep metaData.fnames{curVesInd}];
        vesRawFname = [vesDname filesep metaData.fnames{curVesInd} '_vesRaw.mat'];
        if ~exist(vesRawFname,'file')
            warning('%s does not exists',vesRawFname);
            continue;
        end        
        
        %% Raw vesicle data in both channels & time
        load(vesRawFname);%'curDataTime','dataCh1TimeRaw','dataCh2TimeRaw'
        
        %% Smooth data
        dataCh1TimeSmooth = conv(dataCh1TimeRaw, smoothKernelNorm);
        dataCh2TimeSmooth = conv(dataCh2TimeRaw, smoothKernelNorm);
        
        dataCh1TimeSmooth = dataCh1TimeSmooth(ceil(params.smoothWinSize/2):(end-floor(params.smoothWinSize/2)));
        dataCh2TimeSmooth = dataCh2TimeSmooth(ceil(params.smoothWinSize/2):(end-floor(params.smoothWinSize/2)));

        %% Background correction
        if ~isnan(backgroundModelCh1)
            curBackgroundCh1 = polyval(backgroundModelCh1,curDataTime);
        else
            curBackgroundCh1 = ones(1,length(curDataTime)) .* meanBackgroundCh1;
        end
        dataCh1TimeNorm = (dataCh1TimeSmooth - curBackgroundCh1) ./ backgroundCh1ResidualStd;
        
        
        if ~isnan(backgroundModelCh2)
            curBackgroundCh2 = polyval(backgroundModelCh2,curDataTime);
        else
            curBackgroundCh2 = ones(1,length(curDataTime)) .* meanBackgroundCh2;
        end
        dataCh2TimeNorm = (dataCh2TimeSmooth - curBackgroundCh2) ./ backgroundCh2ResidualStd;
        
        %% Derivative
        derivCh1 = gradient(dataCh1TimeNorm);
        derivCh2 = gradient(dataCh2TimeNorm);
        
        deriv2Ch1 = gradient(derivCh1);
        deriv2Ch2 = gradient(derivCh2);

        %% find maximum location
        [maxValCh1,maxIndCh1] = max(dataCh1TimeNorm);
        [maxValCh2,maxIndCh2] = max(dataCh2TimeNorm);

        %% osillation
        
        %% prepare output for vesicle
        vesData.dataTime = curDataTime;
        
        vesData.dataCh1TimeRaw = dataCh1TimeRaw;
        vesData.dataCh1TimeSmooth = dataCh1TimeSmooth;
        vesData.dataCh1TimeNorm = dataCh1TimeNorm;
        vesData.maxValCh1 = maxValCh1;
        vesData.maxIndCh1 = maxIndCh1;
        vesData.derivCh1 = derivCh1;
        vesData.deriv2Ch1 = deriv2Ch1;

        vesData.dataCh2TimeRaw = dataCh2TimeRaw;
        vesData.dataCh2TimeSmooth = dataCh2TimeSmooth;
        vesData.dataCh2TimeNorm = dataCh2TimeNorm;
        vesData.maxValCh2 = maxValCh2;
        vesData.maxIndCh2 = maxIndCh2;
        vesData.derivCh2 = derivCh2;
        vesData.deriv2Ch2 = deriv2Ch2;
        
        vesData.peaks = DaganFindOscillations(vesData.dataCh1TimeNorm,params,timePerFrame);
        
        vesData.size = vesSize;      
        
        save([vesDname filesep metaData.fnames{curVesInd} '_vesNorm.mat'],'vesData');

        expVesData{iVes} = vesData;

        %% Debug
        % raw signal channel #1
        h = figure;
        hold on;
        plot(curDataTime,dataCh1TimeRaw,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Raw intensity');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_rawCh1.tif']);

        % raw signal channel #2
        h = figure;
        hold on;
        plot(curDataTime,dataCh2TimeRaw,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Raw intensity');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_rawCh2.tif']);

        % normalize signal channel #1
        h = figure;
        hold on;
        plot(curDataTime,dataCh1TimeNorm,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Normalized intensity');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_normCh1.tif']);

        % normalize signal channel #2
        h = figure;
        hold on;
        plot(curDataTime,dataCh2TimeNorm,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Normalized intensity');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_normCh2.tif']);

        close all;
        
        % derivative signal channel #1
        h = figure;
        hold on;
        plot(curDataTime,derivCh1,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Derivative');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_derivCh1.tif']);

        % derivative signal channel #2
        h = figure;
        hold on;
        plot(curDataTime,derivCh2,'ok','MarkerSize',8);        
        xlabel('Frame #');
        ylabel('Derivative');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_derivCh2.tif']);        

        % 2nd derivative signal channel #1
        h = figure;
        hold on;
        plot(curDataTime,deriv2Ch1,'ok','MarkerSize',8);
        xlabel('Frame #');
        ylabel('2nd derivative');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_deriv2Ch1.tif']);
        
        % 2nd derivative signal channel #2
        h = figure;
        hold on;
        plot(curDataTime,deriv2Ch2,'ok','MarkerSize',8);
        xlabel('Frame #');
        ylabel('2nd derivative');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_deriv2Ch2.tif']);

        % Peaks        
        h = figure;
        hold on;
        plot(curDataTime,dataCh1TimeNorm,'ok','MarkerSize',8);        
        plot(curDataTime(vesData.peaks.locs),dataCh1TimeNorm(vesData.peaks.locs),'xr','MarkerSize',15,'LineWidth',2);        
        xlabel('Frame #');
        ylabel('Normalized intensity');
        hold off;
        saveas(h,[vesDname filesep metaData.fnames{curVesInd} '_peaksCh1.tif']);
        
        close all;
    end
    allExpsVesData{iExp} = expVesData;
    allExpsStr{iExp} = metaData.groupsByExperiment{iExp}.experiment{1};
end
save([dname filesep 'out' filesep 'allExpsVesData.mat'],'allExpsVesData','allExpsStr');
end