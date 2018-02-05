function [] = DaganMetaAnalysis_SustainedVesicles(dname,timePerFrame)

close all;
clc;


load([dname filesep 'out' filesep 'allExpsVesData.mat']);%'allExpsVesData','allExpsStr'

ntime = 1000;
maxTime = -inf;

nExps = length(allExpsVesData);

nVesSustained = zeros(1,nExps);
nVes = zeros(1,nExps);

%% get experiments data in matrices
for iExp = 1 : nExps
    expName = allExpsStr{iExp};        
    curExpVesData = allExpsVesData{iExp};
            
    nVes(iExp) = length(curExpVesData);
    
    isSustained()
    
    masterCh1AlignedCh1 = nan(nVes,ntime*2+1);
    masterCh1AlignedCh2 = nan(nVes,ntime*2+1);
    masterCh2AlignedCh1 = nan(nVes,ntime*2+1);
    masterCh2AlignedCh2 = nan(nVes,ntime*2+1);
    
    vesSizes = nan(1,nVes);
    
    vesOscillations = nan(1,nVes);
    vesOscillationsNPeaks = nan(1,nVes);
    vesOscillationsFreq = cell(1,nVes); % each is a list of frequencies
    
    vesMaxSignalCh1 = nan(1,nVes);
    vesMaxSignalCh2 = nan(1,nVes);    
    
    for iVes = 1 : nVes
        curVesData = curExpVesData{iVes};
        
        ch1Data = curVesData.dataCh1TimeNorm;
        ch2Data = curVesData.dataCh2TimeNorm;
        
        ch1MaxInd = curVesData.maxIndCh1;
        ch2MaxInd = curVesData.maxIndCh2;
        
        curVesTime = length(ch1Data);
        curVesMaxTime = max([curVesTime - ch1MaxInd+1,ch1MaxInd,curVesTime - ch2MaxInd+1,ch2MaxInd]);
        if curVesMaxTime > maxTime
            maxTime = curVesMaxTime;
        end
        
        indsTimeCh1 = ((ntime+1)-ch1MaxInd+1) : ((ntime+1)+(curVesTime-ch1MaxInd));
        masterCh1AlignedCh1(iVes,indsTimeCh1) = ch1Data;
        masterCh1AlignedCh2(iVes,indsTimeCh1) = ch2Data;
        
        indsTimeCh2 = ((ntime+1)-ch2MaxInd+1) : ((ntime+1)+(curVesTime-ch2MaxInd));
        masterCh2AlignedCh1(iVes,indsTimeCh2) = ch1Data;
        masterCh2AlignedCh2(iVes,indsTimeCh2) = ch2Data;
        
        vesSizes(iVes) = curVesData.size;
        vesMaxSignalCh1(iVes) = curVesData.maxValCh1;
        vesMaxSignalCh2(iVes) = curVesData.maxValCh2;                        
        
        vesOscillationsNPeaks(iVes) = curVesData.peaks.n;
        if curVesData.peaks.n > 0
            vesOscillations(iVes) = curVesData.peaks.score;            
            vesOscillationsFreq{iVes} = curVesData.peaks.freqs;
            if ~isempty(curVesData.peaks.freqs)
                allOscillationFrequencies = [allOscillationFrequencies curVesData.peaks.freqs];                
            else
                curVesData.peaks.freqs = nan;                
            end
            
            if curVesData.peaks.n > 1
                isOscillating = true;
            else
                isOscillating = false;
            end
            
            allVesOscillationsNPeaks = [allVesOscillationsNPeaks curVesData.peaks.n];
            allVescilesIsOscillating = [allVescilesIsOscillating isOscillating];
            allVescilesTime = [allVescilesTime length(curVesData.dataTime)*timePerFrame(iExp)];
        end
    end
    expsData{iExp}.masterCh1AlignedCh1 = masterCh1AlignedCh1;
    expsData{iExp}.masterCh1AlignedCh2 = masterCh1AlignedCh2;
    expsData{iExp}.masterCh2AlignedCh1 = masterCh2AlignedCh1;
    expsData{iExp}.masterCh2AlignedCh2 = masterCh2AlignedCh2;
    expsData{iExp}.sizes = vesSizes;
    expsData{iExp}.oscillations = vesOscillations;
    expsData{iExp}.oscillationsFreq = vesOscillationsFreq;
    expsData{iExp}.oscillationsNPeaks = vesOscillationsNPeaks;
    expsData{iExp}.vesMaxSignalCh1 = vesMaxSignalCh1;
    expsData{iExp}.vesMaxSignalCh2 = vesMaxSignalCh2;        
end

indsData = ((ntime+1)-maxTime+1):((ntime+1)+maxTime-1);
time = -maxTime+1 : maxTime-1;

%% Oscillation frequencies
h = figure; hold on;
bins = (6:4:38)*0.35;
counts = hist(allOscillationFrequencies(~isnan(allOscillationFrequencies)),bins);
bar(bins,counts);
saveas(h,[dname filesep 'oscillationFreq_1.tif']);
hold off;

h = figure; hold on;
bins = (6:8:38)*0.35;
counts = hist(allOscillationFrequencies(~isnan(allOscillationFrequencies)),bins);
bar(bins,counts);
saveas(h,[dname filesep 'oscillationFreq_2.tif']);
hold off;

timeTH = prctile(allOscillationFrequencies,90); % 15.5
inds = allVescilesTime > (timeTH+3);
nAll = sum(inds);
nOscillating = sum(allVescilesIsOscillating(inds));

fprintf(sprintf('time TH = %.1f, oscillating = %d/%d, median vesicle time = %.1f\n',timeTH,nOscillating,nAll,median(allVescilesTime(inds))));

h = figure; 
hold on;
plot(allVescilesTime,allVesOscillationsNPeaks,'ok','MarkerSize',12,'LineWidth',2);
hold off;
saveas(h,[dname filesep 'timeVsNOscillations.tif']);

[r,p] = corr(allVescilesTime',allVesOscillationsNPeaks');

allOscillationFrequenciesNoNans = allOscillationFrequencies(~isnan(allOscillationFrequencies));
save([dname filesep 'oscillationFrequencies.mat'],'allOscillationFrequenciesNoNans');

% h = figure; hold on;
% counts = hist(allOscillationFrequencies,6:4:38);
% bar(6:4:38,counts./sum(counts));
% saveas(h,[dname filesep 'oscillationFreq.tif']);
% hold off;

%% Plot
caxisVals = [-1,10];
for iExp = 1 : nExps
    expName = allExpsStr{iExp};    
    curExpData = expsData{iExp};
    masterCh1AlignedCh1 = curExpData.masterCh1AlignedCh1(:,indsData); %masterCh1AlignedCh1(masterCh1AlignedCh1 < caxisVals(1)) = nan;
    masterCh1AlignedCh2 = curExpData.masterCh1AlignedCh2(:,indsData); %masterCh1AlignedCh2(masterCh1AlignedCh2 < caxisVals(1)) = nan;
    masterCh2AlignedCh1 = curExpData.masterCh2AlignedCh1(:,indsData); %masterCh2AlignedCh1(masterCh2AlignedCh1 < caxisVals(1)) = nan;
    masterCh2AlignedCh2 = curExpData.masterCh2AlignedCh2(:,indsData); %masterCh2AlignedCh2(masterCh2AlignedCh2 < caxisVals(1)) = nan;
    
    vesSizes = curExpData.sizes;
    
    oscillations = curExpData.oscillations;
    oscillationsNPeaks = curExpData.oscillationsNPeaks;    
    oscillationsFreq = [expsData{iExp}.oscillationsFreq{:}];
    
    [counts,centers] = hist(oscillationsFreq);
    h = figure; hold on; bar(centers,counts); xlabel('Time');ylabel('Count');hold off; saveas(h,[dname expName filesep 'oscillationFreq.tif']);
    
    h = figure; hold on; imagescnan(masterCh1AlignedCh1); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh1AlignedCh1.tif']);
    h = figure; hold on; imagescnan(masterCh1AlignedCh2); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh1AlignedCh2.tif']);
    h = figure; hold on; imagescnan(masterCh2AlignedCh1); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh2AlignedCh1.tif']);
    h = figure; hold on; imagescnan(masterCh2AlignedCh2); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh2AlignedCh2.tif']);
    close all;

    % sorted
    [~,indsSortCh1] = sort(max(masterCh1AlignedCh1,[],2));
    [~,indsSortCh2] = sort(max(masterCh2AlignedCh2,[],2));
    
    [~,indsSortSize] = sort(vesSizes);
    [sortedOscillations,indsSortOscillaions] = sort(oscillations);

    h = figure; hold on; imagescnan(masterCh1AlignedCh1(indsSortCh1,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh1AlignedCh1_sort.tif']);
    h = figure; hold on; imagescnan(masterCh1AlignedCh2(indsSortCh1,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh1AlignedCh2_sort.tif']);
    h = figure; hold on; imagescnan(masterCh2AlignedCh1(indsSortCh2,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh2AlignedCh1_sort.tif']);
    h = figure; hold on; imagescnan(masterCh2AlignedCh2(indsSortCh2,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'masterCh2AlignedCh2_sort.tif']);
    
    h = figure; hold on; imagescnan(masterCh2AlignedCh1(indsSortSize,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'sortBySizeCh1.tif']);
    h = figure; hold on; imagescnan(masterCh2AlignedCh2(indsSortSize,:)); caxis(caxisVals); colorbar; hold off; saveas(h,[dname expName filesep 'sortBySizeCh2.tif']);
    
    h = figure; hold on; plot(1:length(sortedOscillations),sortedOscillations,'ok','MarkerSize',8); hold off; saveas(h,[dname expName filesep 'sortByOscillationsCh1.tif']);
    
    [r1,p1] = corr(curExpData.sizes',curExpData.vesMaxSignalCh1');
    h = figure; hold on; plot(curExpData.sizes,curExpData.vesMaxSignalCh1,'ok','MarkerSize',8); 
    title(sprintf('Corr size vs. ch1 intensity: r = %.2f, p = %.5f',r1,p1));
    xlabel('Size'); ylabel('Signal Ch1'); hold off; saveas(h,[dname expName filesep 'sizeVsSignalCh1.tif']);
    [r2,p2] = corr(curExpData.sizes',curExpData.vesMaxSignalCh2');
    h = figure; hold on; plot(curExpData.sizes,curExpData.vesMaxSignalCh2,'ok','MarkerSize',8); 
    title(sprintf('Corr size vs. ch2 intensity: r = %.2f, p = %.5f',r2,p2));
    xlabel('Size'); ylabel('Signal Ch2'); hold off; saveas(h,[dname expName filesep 'sizeVsSignalCh2.tif']);
    close all;

    [rho,pval] = corr(max(masterCh1AlignedCh1,[],2),max(masterCh2AlignedCh2,[],2));
    h = figure; hold on; plot(max(masterCh1AlignedCh1,[],2),max(masterCh2AlignedCh2,[],2),'ok','MarkerSize',8); xlabel('Ch1'); ylabel('Ch2'); 
    title(sprintf('corr = %g, pval = %g',rho,pval)); hold off; saveas(h,[dname expName filesep 'maxCh1VsCh2.tif']);

    plotAligned(time,masterCh1AlignedCh1,[dname expName filesep 'masterCh1AlignedCh1_overlay.tif']);
    plotAligned(time,masterCh1AlignedCh2,[dname expName filesep 'masterCh1AlignedCh2_overlay.tif']);
    plotAligned(time,masterCh2AlignedCh1,[dname expName filesep 'masterCh2AlignedCh1_overlay.tif']);
    plotAligned(time,masterCh2AlignedCh2,[dname expName filesep 'masterCh2AlignedCh2_overlay.tif']);
    
    save([dname expName filesep 'alignedMetaData.mat'],'time','masterCh1AlignedCh1','masterCh1AlignedCh2','masterCh2AlignedCh1','masterCh2AlignedCh2');
end
end

function [] = plotAligned(time,data,outFname)
nVes = size(data,1);
h = figure;
hold on;
for iVes = 1 : nVes
    plot(time,data(iVes,:),'k--');
end
xlabel('Frame #');
ylabel('Normalized intensity');
hold off;
saveas(h,outFname);
end






