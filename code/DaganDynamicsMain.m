function [] = DaganDynamicsMain()

close all;

dname = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\20170301\data_rock_inh\';
actinFname = 'rock_inh_actin.csv';
rhoFanme = 'rock_inh_rho.csv';

thresholds = 0:0.05:1;

outActin = analyzeDaganDynamics(dname,actinFname,thresholds);
outRho = analyzeDaganDynamics(dname,rhoFanme,thresholds);

figure; imagesc(outRho.up(1:20,:)); caxis([0,1]); colorbar; title('RhoA up');
figure; imagesc(outActin.up(1:20,:)); caxis([0,1]); colorbar; title('Actin up');
figure; imagesc(outRho.down(1:20,:)); caxis([0,1]); colorbar; title('RhoA down');
figure; imagesc(outActin.down(1:20,:)); caxis([0,1]); colorbar; title('Actin down');
% 
% figure; imagesc(outRho.up); caxis([0,0.1]); colorbar; title('RhoA up');
% figure; imagesc(outActin.up); caxis([0,0.1]); colorbar; title('Actin up');
% figure; imagesc(outRho.down); caxis([0,0.1]); colorbar; title('RhoA down');
% figure; imagesc(outActin.down); caxis([0,0.1]); colorbar; title('Actin down');
% 
% figure; 
% hold on; 
% plot(thresholds(1:end),outRho.upEMDRef,'--g','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outRho.downEMDRef,'--r','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outActin.upEMDRef,'--c','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end),outActin.downEMDRef,'--m','MarkerSize',8,'LineWidth',2);
% xlabel('Relative intensity'); ylabel('Recruitment time (EMD to ref)');
% legend({'RhoA up','RhoA down','Actin up','Actin down'});
% hold off;
% 
% figure; 
% hold on; 
% plot(thresholds(1:end-1),outRho.upEMD,'--g','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outRho.downEMD,'--r','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outActin.upEMD,'--c','MarkerSize',8,'LineWidth',2);
% plot(thresholds(1:end-1),outActin.downEMD,'--m','MarkerSize',8,'LineWidth',2);
% xlabel('Relative intensity'); ylabel('{\Delta}Recruitment time (EMD)');
% legend({'RhoA up','RhoA down','Actin up','Actin down'});
% hold off;

figure; 
hold on; 
plot(thresholds(1:end),outRho.upMedian,'--g','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end),outRho.downMedian,'--r','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end),outActin.upMedian,'--c','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end),outActin.downMedian,'--m','MarkerSize',8,'LineWidth',2);
xlabel('Relative intensity'); ylabel('Median recruitment time');
legend({'RhoA up','RhoA down','Actin up','Actin down'});
hold off;

figure; 
hold on; 
plot(thresholds(1:end-1),outRho.upMedian(1:end-1)-outRho.upMedian(2:end),'--g','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end-1),outRho.downMedian(1:end-1)-outRho.downMedian(2:end),'--r','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end-1),outActin.upMedian(1:end-1)-outActin.upMedian(2:end),'--c','MarkerSize',8,'LineWidth',2);
plot(thresholds(1:end-1),outActin.downMedian(1:end-1)-outActin.downMedian(2:end),'--m','MarkerSize',8,'LineWidth',2);
xlabel('Relative intensity'); ylabel('{\Delta}Median recruitment time');
legend({'RhoA up','RhoA down','Actin up','Actin down'});
hold off;



end