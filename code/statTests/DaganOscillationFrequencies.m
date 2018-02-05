function [] = DaganOscillationFrequencies()

close all;
clc;

%% experiments
dname = 'C:\Users\assafza\Google Drive\Research\PostDoc\Collaborations\Dagan\';


rockDname = [dname 'rockinh\'];
rockArpDname = [dname  'rockinh_arpinh\'];

rockOscFreq = load([rockDname 'oscillationFrequencies.mat']);
rockArpOscFreq = load([rockArpDname 'oscillationFrequencies.mat']);

rockOscFreq = rockOscFreq.allOscillationFrequenciesNoNans;
rockArpOscFreq = rockArpOscFreq.allOscillationFrequenciesNoNans;

rockOscInds = find(~isnan(rockOscFreq));
rockArpOscInds = find(~isnan(rockArpOscFreq));

pval = ranksum(rockOscFreq(rockOscInds),rockArpOscFreq(rockArpOscInds));

fold = mean(rockArpOscFreq(rockArpOscInds))/mean(rockOscFreq(rockOscInds));

fprintf(sprintf('pval (n = %d,%d) = %.4f, fold = %.2f\n',length(rockOscInds),length(rockArpOscInds),pval,fold));

end