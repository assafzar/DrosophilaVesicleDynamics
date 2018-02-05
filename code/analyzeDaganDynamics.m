function out = analyzeDaganDynamics(dname,fanme,thresholds)

nthresholds = length(thresholds);

data = csvread([dname fanme]);

[ntime,n] = size(data);

out.up = zeros(ntime,nthresholds);
out.down = zeros(ntime,nthresholds);

out.upEMD = nan(1,nthresholds-1);
out.downEMD = nan(1,nthresholds-1);

out.upEMDRef = nan(1,nthresholds);
out.downEMDRef = nan(1,nthresholds);

out.upMedian = nan(1,nthresholds);
out.downMedian = nan(1,nthresholds);



for ith = 1 : nthresholds 
    allUp = [];
    allDown = [];
    curTH = thresholds(ith);
    for io = 1 : n
        obs = data(:,io);
        indMax = find(obs == 1);
        maxVal = obs(indMax);
        TH = maxVal * curTH;
        
        %% going up 
        for iup = indMax-1 : -1 : 3 % 2
            if TH >= obs(iup)
                % (2 time points below threshold)
                if TH >= obs(iup-1) && TH >= obs(iup-2)
                    out.up(indMax-iup,ith) = out.up(indMax-iup,ith) + 1;
                    allUp = [allUp, indMax-iup];
                    break;
                end
            end
            if iup == 3
                out.up(end,ith) = out.up(end,ith) + 1;
                allUp = [allUp, inf];
                break;
            end
        end        
        
        %% going down
        for idown = indMax+1 : ntime-2 % 1
            if TH > obs(idown)
                if TH >= obs(idown+1) && TH >= obs(idown+2) 
                    out.down(idown-indMax,ith) = out.down(idown-indMax,ith) + 1;
                    allDown = [allDown, idown-indMax];
                    break;
                end
            end
            if idown == ntime-2
                out.down(end,ith) = out.down(end,ith) + 1;
                allDown = [allDown, inf];
                break;
            end
        end
    end
    
    %     assert(n == length(allUp));
    %     assert(n == length(allDown));
    
    out.up(:,ith) = out.up(:,ith) ./ sum(out.up(:,ith)); 
    out.down(:,ith) = out.down(:,ith) ./ sum(out.down(:,ith));
    out.upMedian(ith) = median(allUp);
    out.downMedian(ith) = median(allDown);
end

%% EMD
for ith = 1 : nthresholds-1  
    out.upEMD(ith) = EmdDagan(out.up(:,ith),out.up(:,ith+1));
    out.downEMD(ith) = EmdDagan(out.down(:,ith),out.down(:,ith+1));
end

%% EMD reference
for ith = 1 : nthresholds 
    referenceDistribution = zeros(size(out.up,1),1);
    referenceDistribution(1) = 1;
    out.upEMDRef(ith) = EmdDagan(out.up(:,ith),referenceDistribution);
    out.downEMDRef(ith) = EmdDagan(out.down(:,ith),referenceDistribution);
end


end

function similarity= EmdDagan(x,y)
similarity = sum(abs(cumsum(x)-cumsum(y)));
end