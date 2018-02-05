
function peaks = DaganFindOscillations(data,params,timePerFrame)
[peaks.pks,peaks.locs]  = findpeaks(data,'MinPeakHeight',params.MinPeakHeight,'MinPeakDistance',params.MinPeakDistance);
peaks = filterMinPeakProminence(data,peaks,params.MinPeakProminence);
peaks.n = length(peaks.locs);
peaks.score = (peaks.n - 1)./ length(data);
if peaks.n > 1
    peaks.freqs = (peaks.locs(2:end) - peaks.locs(1:end-1)) * timePerFrame;
else
    peaks.freqs = nan;
end
end


%% Wrong implementation - MinPeakProminence means drop off on both sides before encountering a larger value
function peaks = filterMinPeakProminence(data,peaks,MinPeakProminence)

newLocs = [];
nPeaks = length(peaks.locs);

if nPeaks == 1
    return;
end

for iloc = 1 : nPeaks
    curPeakVal = peaks.pks(iloc);
    curPeakLoc = peaks.locs(iloc);
    includePeak = true;
            
    if iloc == 1
        if curPeakVal < peaks.pks(iloc+1) && ...
                curPeakVal < min(data(peaks.locs(iloc):peaks.locs(iloc+1))) + MinPeakProminence
            includePeak = false;            
        end
    else if iloc == nPeaks
            if curPeakVal < peaks.pks(iloc-1) && ...
                    curPeakVal < min(data(peaks.locs(iloc-1):peaks.locs(iloc))) + MinPeakProminence
                includePeak = false;                
            end
        else
            if (curPeakVal < peaks.pks(iloc-1) || curPeakVal < peaks.pks(iloc+1)) && ...
                    ((curPeakVal < min(data(peaks.locs(iloc):peaks.locs(iloc+1))) + MinPeakProminence) || ...
                    (curPeakVal < min(data(peaks.locs(iloc-1):peaks.locs(iloc))) + MinPeakProminence))
                includePeak = false;                
            end
        end
    end
    if includePeak
        newLocs = [newLocs curPeakLoc];
    end
end

peaks.locs = newLocs;
peaks.pks = data(peaks.locs);

end

% 
%    % back 
%     for ilocBack = iloc : -1 : 1
%         if peaks.locs(ilocBack)  < curPeakLoc - MinPeakProminence
%             break;
%         end
%         if peaks.pks(ilocBack) > curPeakVal
%             includePeak = false;
%             break;
%         end
%     end
%     % forth
%     for ilocForth = iloc : nPeaks
%         if peaks.locs(ilocForth)  > curPeakLoc + MinPeakProminence
%             break;
%         end
%         if peaks.pks(ilocForth) > curPeakVal
%             includePeak = false;
%             break;
%         end
%     end