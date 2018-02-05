function [] = DaganDynamicsMainSynch()

close all;

%% data
% dname = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\20170301\data_rock_inh\';
% actinFname = 'rock_inh_actin.csv';
% rhoFanme = 'rock_inh_rho.csv';
dname = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\20170301\selectByPeaks\';
actinFname = 'select_data_actin.csv';
rhoFanme = 'select_data_rho.csv';

%% parameters
peakIntervalSearchRadius = 5;
signalFracTH = 0.2;
timeLag = -5:5;

%%
actinData = csvread([dname actinFname]);
rhoData = csvread([dname rhoFanme]);

[ntimeActin,nActin] = size(actinData);
[ntimeRho,nRho] = size(rhoData);

actinRhoPeaks = detectPeaks(actinData,rhoData,peakIntervalSearchRadius);
rhoActinPeaks = detectPeaks(rhoData,actinData,peakIntervalSearchRadius);

actinRhoIntervals = getIntervals(actinData,rhoData,actinRhoPeaks,signalFracTH);
rhoActinIntervals = getIntervals(rhoData,actinData,rhoActinPeaks,signalFracTH);

% todo: plot intervals

actinRhoTemporalCorr = getTemporalCorr(actinData,rhoData,actinRhoIntervals,timeLag,'Actin->Rho');
rhoActinTemporalCorr = getTemporalCorr(rhoData,actinData,rhoActinIntervals,timeLag,'Rho->Actin');

end

%%
function peaks = detectPeaks(master,slave,peakIntervalSearchRadius)
[ntimeMaster,nobsMaster] = size(master);
[ntimeSlave,nobsSlave] = size(slave);

assert(ntimeMaster == ntimeSlave && nobsMaster == nobsSlave);

ntime = ntimeMaster;
nobs = nobsMaster;

peaks.masterVal = cell(1,nobs);
peaks.masterInd = cell(1,nobs);
peaks.slaveVal = cell(1,nobs);
peaks.slaveInd = cell(1,nobs);

for io = 1 : nobs
    curMaster = master(:,io);
    curSlave = slave(:,io);
    indMasterMax = find(curMaster == 1);
    
    % exclude these lines
    maxVal = curMaster(indMasterMax);    
    
    % find max slave in the interval
    indsSlaveInterval = max(indMasterMax-peakIntervalSearchRadius,1): min(indMasterMax+peakIntervalSearchRadius,ntime);
    tmpSlave = curSlave; tmpSlave(~indsSlaveInterval) = -inf;
    [slaveMaxVal,indSlaveMax] = max(tmpSlave);
    
    peaks.masterVal{io} = 1;
    peaks.masterInd{io} = indMasterMax;
    peaks.slaveVal{io} = slaveMaxVal;
    peaks.slaveInd{io} = indSlaveMax;
end

end

%%
function interval = getIntervals(master,slave,peaks,signalFracTH)
interval.master = getInterval(master,peaks.masterInd,peaks.masterVal,signalFracTH);
interval.slave = getInterval(slave,peaks.slaveInd,peaks.slaveVal,signalFracTH);
end

%%
function interval = getInterval(data,maxInds,maxVals,signalFracTH)
[ntime,nobs] = size(data);

interval.start = cell(1,nobs);
interval.end = cell(1,nobs);

for io = 1 : nobs
    curData = data(:,io);
    maxVal = maxVals{io};
    maxInd = maxInds{io};
    TH = signalFracTH * maxVal;
    
    %% begining of interval
    if maxInd <= 3
        if TH >= curData(2)
                interval.start{io} = 2;
            else
                interval.start{io} = 1;
        end
    end
    
    for ipast = maxInd-1 : -1 : 3 % 2
        if TH >= curData(ipast)
            % (2 time points below threshold)
            if TH >= curData(ipast-1) && TH >= curData(ipast-2)
                interval.start{io} = ipast;
                break;
            end
        end
        if ipast == 3
            if TH >= curData(2)
                interval.start{io} = 2;
            else
                interval.start{io} = 1;
            end
        end
    end
    
    %% end of interval
    assert(maxInd+1 <= ntime-2);
    
    for ifuture = maxInd+1 : ntime-2 % 1
        if TH > curData(ifuture)
            if TH >= curData(ifuture+1) && TH >= curData(ifuture+2)
                interval.end{io} = ifuture;
                break;
            end
        end
        if ifuture == ntime-2
            if TH >= curData(ntime-1)
                interval.end{io} = ntime-1;
            else
                interval.end{io} = ntime;
            end
        end
    end    
end
end


function temporalCorr = getTemporalCorr(master,slave,intervals,timeLag,masterSlaveStr)
[ntime,nobs] = size(master);
nTimeLag = length(timeLag);

fprintf('\n\n');
fprintf('%s:\n',masterSlaveStr);

temporalCorr = cell(1,nobs);

for io = 1 : nobs
    curMasterData = master(:,io);
    curSlaveData = slave(:,io);  
        
    masterStartInd = intervals.master.start{io};
    masterEndInd = intervals.master.end{io};
    
    masterInterval = masterStartInd:masterEndInd;
    
    temporalCorr{io}.corr = nan(1,nTimeLag);
    temporalCorr{io}.maxCorr = nan;
    temporalCorr{io}.timeLag = nan;
    
    for iTimeLag = timeLag
        iCorr = iTimeLag - timeLag(1) + 1;
        
        if (masterStartInd + iTimeLag < 1) || (masterEndInd - iTimeLag > ntime)
            continue;
        end
                
        slaveInterval = (masterStartInd + iTimeLag):(masterEndInd + iTimeLag);
        
        temporalCorr{io}.corr(iCorr) = corr(curMasterData(masterInterval),curSlaveData(slaveInterval));
    end
        
    if ~isnan(temporalCorr{io}.corr(ceil(nTimeLag/2)+1)) && ~isnan(temporalCorr{io}.corr(ceil(nTimeLag/2)-1)) %% making sure that we can go both directions!
        [temporalCorr{io}.maxCorr,indMax] = max(temporalCorr{io}.corr);
        temporalCorr{io}.timeLag = timeLag(indMax);
        fprintf('obs %d: time lag = %.2f, corr = %.2f\n',io,temporalCorr{io}.timeLag,temporalCorr{io}.maxCorr);
    else
        fprintf('obs %d: not enough data\n',io);  
    end
end

end


% assert(ntimeActin == ntimeRho && nActin == nRho);
%
% thresholds = 0:0.05:1;
%
% outActin = analyzeDaganDynamics(dname,actinFname,thresholds);
% outRho = analyzeDaganDynamics(dname,rhoFanme,thresholds);
%
% figure; imagesc(outRho.up(1:20,:)); caxis([0,1]); colorbar; title('RhoA up');
% figure; imagesc(outActin.up(1:20,:)); caxis([0,1]); colorbar; title('Actin up');
% figure; imagesc(outRho.down(1:20,:)); caxis([0,1]); colorbar; title('RhoA down');
% figure; imagesc(outActin.down(1:20,:)); caxis([0,1]); colorbar; title('Actin down');
%
% figure;
% hold on;
% plot(thresholds(1:end),outRho.upMedian,'--g','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outRho.downMedian,'--r','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outActin.upMedian,'--c','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outActin.downMedian,'--m','MarkerSize',8,'LineWidth',2);
% xlabel('Relative intensity'); ylabel('Median recruitment time');
% legend({'RhoA up','RhoA down','Actin up','Actin down'});
% hold off;
%
% figure;
% hold on;
% plot(thresholds(1:end-1),outRho.upMedian(1:end-1)-outRho.upMedian(2:end),'--g','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outRho.downMedian(1:end-1)-outRho.downMedian(2:end),'--r','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outActin.upMedian(1:end-1)-outActin.upMedian(2:end),'--c','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outActin.downMedian(1:end-1)-outActin.downMedian(2:end),'--m','MarkerSize',8,'LineWidth',2);
% xlabel('Relative intensity'); ylabel('{\Delta}Median recruitment time');
% legend({'RhoA up','RhoA down','Actin up','Actin down'});
% hold off;
%
%
%
% end