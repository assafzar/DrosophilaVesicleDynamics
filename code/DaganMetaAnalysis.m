function [] = DaganMetaAnalysis(dname,timePerFrame)

close all;
clc;

% single vesicle analysis parameters
singleVesParams.deltaZscore = 2; % # stds above that vesicle background
singleVesParams.nBackgroundFrames = 2; % number of frames to consider as background

%
% if nargin < 1
%     %     dname = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\';
%     dname = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\arpinh\';
% end

load([dname filesep 'out' filesep 'allExpsVesData.mat']);%'allExpsVesData','allExpsStr'

ntime = 1000;
maxTime = -inf;

nExps = length(allExpsVesData);


%% get experiments data in matrices
expsData = cell(1,nExps);

nVesSustained = zeros(1,nExps);
nVes = zeros(1,nExps);

allOscillationFrequencies = [];
allVescilesIsOscillating = [];
allVescilesTime = [];
allVesOscillationsNPeaks = [];

allVesCycleTimeActin = [];
allVesSustainedActin = [];

allVesCycleTimeRho = [];
allVesSustainedRho = [];

allVesOnsetTimeLag = [];
allVesDeclineTimeLag = [];
allVesPeakTimeLag = [];

allVesPeakZScoreActin = [];
allVesPeakZScoreRho = [];

allVesTEndActin = [];
allVesTEndRho = [];

for iExp = 1 : nExps
    expName = allExpsStr{iExp};    
    expDname = [dname expName];
    
    sustainedDname = [expDname filesep 'sustained'];
    if ~exist(sustainedDname,'dir')
        mkdir(sustainedDname);
    end
    
    curExpVesData = allExpsVesData{iExp};
    
    nVes(iExp) = length(curExpVesData);
    
    masterCh1AlignedCh1 = nan(nVes(iExp),ntime*2+1);
    masterCh1AlignedCh2 = nan(nVes(iExp),ntime*2+1);
    masterCh2AlignedCh1 = nan(nVes(iExp),ntime*2+1);
    masterCh2AlignedCh2 = nan(nVes(iExp),ntime*2+1);
    
    vesSizes = nan(1,nVes(iExp));
    
    vesOscillations = nan(1,nVes(iExp));
    vesOscillationsNPeaks = nan(1,nVes(iExp));
    vesOscillationsFreq = cell(1,nVes(iExp)); % each is a list of frequencies
    
    vesMaxSignalCh1 = nan(1,nVes(iExp));
    vesMaxSignalCh2 = nan(1,nVes(iExp));    
    
    for iVes = 1 : nVes(iExp)
        curVesData = curExpVesData{iVes};
        
        ch1Data = curVesData.dataCh1TimeNorm;
        ch2Data = curVesData.dataCh2TimeNorm;
        
        ch1MaxInd = curVesData.maxIndCh1;
        ch2MaxInd = curVesData.maxIndCh2;
        
        %% sustained        
        % Actin
        [isSustainedVesicleActin,sustainedFramesActin,tendActin] = isSustained(ch1Data,ch1MaxInd,timePerFrame(iExp),sustainedDname,singleVesParams.deltaZscore,singleVesParams.nBackgroundFrames);
        if(isSustainedVesicleActin)
            nVesSustained(iExp) = nVesSustained(iExp) + 1;
        end
        if ~isnan(sustainedFramesActin)            
            allVesCycleTimeActin = [allVesCycleTimeActin sustainedFramesActin*timePerFrame(iExp)];            
            allVesSustainedActin = [allVesSustainedActin isSustainedVesicleActin];                        
        end
        
        % Rho
        [isSustainedVesicleRho,sustainedFramesRho,tendRho] = isSustained(ch2Data,ch2MaxInd,timePerFrame(iExp),sustainedDname,singleVesParams.deltaZscore,singleVesParams.nBackgroundFrames);
        if(isSustainedVesicleRho)
            nVesSustained(iExp) = nVesSustained(iExp) + 1;
        end
        if ~isnan(sustainedFramesRho)            
            allVesCycleTimeRho = [allVesCycleTimeRho sustainedFramesRho*timePerFrame(iExp)];            
            allVesSustainedRho = [allVesSustainedRho isSustainedVesicleRho];                        
        end
        %%
        
        if ~isnan(sustainedFramesActin) && ~isnan(sustainedFramesRho)
            lagTimeOnset = getOnsetLagTime(ch1Data,ch2Data,ch1MaxInd,ch2MaxInd,singleVesParams.deltaZscore,singleVesParams.nBackgroundFrames);
            lagTimeDecline = getDeclineLagTime(ch1Data,ch2Data,ch1MaxInd,ch2MaxInd,singleVesParams.deltaZscore+2);
            lagTimePeak = ch1MaxInd - ch2MaxInd;
            if ~isnan(lagTimeOnset)
                allVesOnsetTimeLag = [allVesOnsetTimeLag lagTimeOnset*timePerFrame(iExp)];
            end
            if ~isnan(lagTimeDecline)
                allVesDeclineTimeLag = [allVesDeclineTimeLag lagTimeDecline*timePerFrame(iExp)];
            end
            allVesPeakTimeLag = [allVesPeakTimeLag lagTimePeak*timePerFrame(iExp)];
            
            allVesPeakZScoreActin = [allVesPeakZScoreActin ch1Data(ch1MaxInd)];
            allVesPeakZScoreRho = [allVesPeakZScoreRho ch2Data(ch2MaxInd)];
            
            allVesTEndActin = [allVesTEndActin tendActin*timePerFrame(iExp)];
            allVesTEndRho = [allVesTEndRho tendRho*timePerFrame(iExp)];
        end
        
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
    
    
    %% Sustained
    fprintf(sprintf('\n %s\n',expName));
    fprintf(sprintf('Sustatrueined vesicles: %d/%d\n',nVesSustained(iExp),nVes(iExp)));            
end

fprintf(sprintf('\nExperiment cycle time stas:\n Actin (n = %d)): mean(%.1f), std(%.1f)\n Rho (n = %d)): mean(%.1f), std(%.1f)\n',...
    length(allVesCycleTimeActin),mean(allVesCycleTimeActin),std(allVesCycleTimeActin),...
    length(allVesCycleTimeRho),mean(allVesCycleTimeRho),std(allVesCycleTimeRho)));    

pOnsetTimeLag = signrank(allVesOnsetTimeLag);
fprintf(sprintf('\nOnset lag time stas (n = %d)): mean(%.3f), std(%.3f), p-value(%.5f)\n',length(allVesOnsetTimeLag),mean(allVesOnsetTimeLag),std(allVesOnsetTimeLag),pOnsetTimeLag));    
pDeclineTimeLag = signrank(allVesDeclineTimeLag);
fprintf(sprintf('Decline lag time stas (n = %d)): mean(%.3f), std(%.3f), p-value(%.5f)\n',length(allVesDeclineTimeLag),mean(allVesDeclineTimeLag),std(allVesDeclineTimeLag),pDeclineTimeLag));    
pPeakTimeLag = signrank(allVesPeakTimeLag);
fprintf(sprintf('Peak lag time stas (n = %d)): mean(%.3f), std(%.3f), p-value(%.5f)\n',length(allVesPeakTimeLag),mean(allVesPeakTimeLag),std(allVesPeakTimeLag),pPeakTimeLag));    

pvalPeakZScore = signrank(allVesPeakZScoreActin-allVesPeakZScoreRho);
fprintf(sprintf('\nPeak z-score stas Actin (n = %d)): mean(%.1f), std(%.1f)\n',length(allVesPeakZScoreActin),mean(allVesPeakZScoreActin),std(allVesPeakZScoreActin)));    
fprintf(sprintf('Peak z-score stas Rho (n = %d)): mean(%.1f), std(%.1f)\n',length(allVesPeakZScoreRho),mean(allVesPeakZScoreRho),std(allVesPeakZScoreRho)));    
fprintf(sprintf('Peak z-score Actin vs. Rho p-value: p-value(%.5f)\n',pvalPeakZScore));    

allVesTEndTimeLag = allVesTEndActin - allVesTEndRho;
pTEndTimeLag = signrank(allVesTEndTimeLag);
fprintf(sprintf('\n End cycle lag time stats (n = %d)): mean(%.1f), std(%.1f), p-value(%.5f)\n',...
    length(allVesTEndTimeLag),mean(allVesTEndTimeLag),std(allVesTEndTimeLag),pTEndTimeLag));    

fprintf(sprintf('\n\n video length (pooled): mean(%.1f), std(%.1f)\n',mean(allVescilesTime),std(allVescilesTime)));

save([dname filesep 'out' filesep 'allVesCycleTime.mat'],...
    'allVesCycleTimeActin','allVesSustainedActin',...
    'allVesCycleTimeRho','allVesSustainedRho',...
    'allVesOnsetTimeLag','allVesDeclineTimeLag','allVesPeakTimeLag',...
    'pOnsetTimeLag','pDeclineTimeLag','pPeakTimeLag',...
    'allVesPeakZScoreActin','allVesPeakZScoreRho','pvalPeakZScore',...
    'allVesTEndTimeLag','pTEndTimeLag');

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

%% 

function [isSustainedVesicleActin,nFrames,tend] = isSustained(data,maxInd,timePerFrame,sustainedDname,deltaZscore,nBackgroundFrames)

assert(maxInd > nBackgroundFrames);
backgroundZScore = mean(data(1:nBackgroundFrames));

signalAboveBackground = (data > (backgroundZScore + deltaZscore));

if(sum(signalAboveBackground) == 0)
    warning('no signal above background of %.1f + %.1f',backgroundZScore,deltaZscore);
    isSustainedVesicleActin = false;
    nFrames = nan;
    tend = nan;
    return;
end


tstrat = nan;
tend = nan;
for i = (maxInd-1) : -1 : 1
    if ~signalAboveBackground(i)
        tstart = i + 1;
        break;
    end
end
assert(~isnan(tstart));

for i = (maxInd+1) : length(signalAboveBackground)
    if ~signalAboveBackground(i)
        tend = i - 1;
        break;
    end
end
if isnan(tend)
    tend = length(signalAboveBackground);
end
% tstart = find(signalAboveBackground,1,'first');

nFrames = tend - tstart + 1;
isSustainedVesicleActin = sum(signalAboveBackground(tstart:end)) == (length(signalAboveBackground)-tstart+1);

% h = figure; hold on; plot(timePerFrame*(1:length(data)),data,'or'); plot(timePerFrame*(tstart:tend),data(tstart:tend),'*b');
% for i = 1:100
%     outFname = [sustainedDname filesep num2str(i) 'a.jpg'];
%     if ~exist(outFname,'file')
%         saveas(h,outFname);
%         break;
%     end
% end
% close all;
end

%%
function lagTime = getOnsetLagTime(ch1Data,ch2Data,ch1MaxInd,ch2MaxInd,deltaZscore,nBackgroundFrames)

onsetCh1 = getOnsetFrame(ch1Data,ch1MaxInd,deltaZscore,nBackgroundFrames);
onsetCh2 = getOnsetFrame(ch2Data,ch2MaxInd,deltaZscore,nBackgroundFrames);

if isnan(onsetCh1) || isnan(onsetCh2)
    lagTime = nan;
else
    lagTime = onsetCh1 - onsetCh2;
end
end


function onsetT = getOnsetFrame(data,maxInd,deltaZscore,nBackgroundFrames)
backgroundZScore = mean(data(1:nBackgroundFrames));

signalAboveBackground = (data > (backgroundZScore + deltaZscore));

if(sum(signalAboveBackground) == 0)
    onsetT = nan;
    return;
end


onsetT = nan;
for i = (maxInd-1) : -1 : 1
    if ~signalAboveBackground(i)
        onsetT = i + 1;
        break;
    end
end
assert(~isnan(onsetT));
end

%%
function lagTime = getDeclineLagTime(ch1Data,ch2Data,ch1MaxInd,ch2MaxInd,deltaZscore)

declineCh1 = getDeclineFrame(ch1Data,ch1MaxInd,deltaZscore);
declineCh2 = getDeclineFrame(ch2Data,ch2MaxInd,deltaZscore);

if isnan(declineCh1) || isnan(declineCh2)
    lagTime = nan;
else
    lagTime = declineCh1 - declineCh2;
end
end


function declineT = getDeclineFrame(data,maxInd,deltaZscore)

peakZScore = data(maxInd);

signalBelowPeak = (data < (peakZScore - deltaZscore));

if(sum(signalBelowPeak(maxInd:end)) == 0)
    declineT = nan;
    return;
end

declineT = nan;
for i = maxInd : length(signalBelowPeak)
    if signalBelowPeak(i)
        declineT = i;
        break;
    end
end
assert(~isnan(declineT));
end




