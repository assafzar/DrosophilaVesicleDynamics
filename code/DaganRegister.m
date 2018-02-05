function [] = DaganRegister(tifStackDname, fname,params)

close all;

inFname = [tifStackDname fname];

info = imfinfo(inFname);
nTime = numel(info);

[tAnnotated,roiAnnotated] = findAnnotatedFrame(inFname,info,params.dilateR); 

vesSize = size(roiAnnotated,1)*size(roiAnnotated,2);

if tAnnotated < 3
    warning([inFname ': tAnnotate = ' num2str(tAnnotated)]);
end

[translations,ROIS,refLines] = calcTranslationsRefFrame(inFname,info,params,nTime,roiAnnotated,tAnnotated);
dxs = translations.dx;
dys = translations.dy;
scores = translations.score;

% tFuture = tAnnotated : nTime;
% tPast = tAnnotated : -1 : 1;
% [pastTranslations,pastROI] = calcTranslations(inFname,info,params,tPast,roiAnnotated);
% [futureTranslations,futureROI] = calcTranslations(inFname,info,params,tFuture,roiAnnotated);
% 
% ROIS = [pastROI(end:-1:1),futureROI(2:end)];
% 
% [dys,dxs,scores] = integrateTrans(pastTranslations,futureTranslations);

[imgs,imgsTrans] = applyTranslations(inFname,info,dys,dxs,params);

outdname = [tifStackDname fname(1:end-4)];
if ~exist(outdname,'dir')
    mkdir(outdname);
end

if exist([outdname filesep fname(1:end-4) '_out.mat'],'file') && ~params.always
    return;
end

outROIsFname = [outdname filesep fname(1:end-4) '_roisOut.tif'];
outCropFname = [outdname filesep fname(1:end-4) '_croppedOut.tif'];
% outImgsFname = [tifStackDname fname(1:end-4) '_imgsOut.tif'];
% outImgsTransROIsFname = [tifStackDname fname(1:end-4) '_imgsTransRoisOut.tif'];

kymographCh1Fname = [outdname filesep fname(1:end-4) '_kymographCh1.tif'];
kymographCh2Fname = [outdname filesep fname(1:end-4) '_kymographCh2.tif'];

outCorrFname = [outdname filesep fname(1:end-4) '_corr.tif'];

if exist(outROIsFname,'file')
    delete(outROIsFname);
    delete(outCropFname);
    delete(outCorrFname);
    delete(kymographCh1Fname);
    delete(kymographCh2Fname);
end

dataCh1 = [];
dataCh2 = [];
for i = 1 : nTime
    [Iroi,Icrop,refLineData,refLineDataCrop] = getComposite(imgs{i},ROIS{i},refLines{i},params);
    if i == 1
        imwrite(Iroi,outROIsFname);
        imwrite(Icrop,outCropFname);
    else
        imwrite(Iroi,outROIsFname, 'writemode', 'append');
        imwrite(Icrop,outCropFname,'writemode', 'append');
    end
    
    dataCh1 = [dataCh1;refLineData.dataCh1];
    dataCh2 = [dataCh2;refLineData.dataCh2];
    
    % Buggy code :-( Something does not work with the translation??
    %     imwrite(imgsTrans{i},outImgsFname, 'writemode', 'append'); % works for 1st iteration??
    %     ITransRoi = getComposite(imgsTrans{i},pastROI{end});
    %     imwrite(ITransRoi,outImgsTransROIsFname, 'writemode', 'append');
    %     disp([num2str(i) ': (dy,dx) = (' num2str(dys(i)) ',' num2str(dxs(i)) '), score = ' num2str(scores(i))]);
end


h = figure;
hold on;
plot(1:length(scores),scores,'ok','MarkerSize',8);
xlabel('Frame');
ylabel('Correlation');
hold off;
saveas(h,outCorrFname);

% Todo: add caxis based on the video! (move this to meta analysis)
plotKymograph(dataCh1,kymographCh1Fname);
plotKymograph(dataCh2,kymographCh2Fname);

outFname = [outdname filesep fname(1:end-4) '_out.mat'];
save(outFname,'dataCh1','dataCh2','vesSize','tAnnotated');

end

%%
function [translations,bbs,refLines] = calcTranslationsRefFrame(inFname,info,params,ntime, roiAnnotated,tAnnotated)
stats = regionprops(roiAnnotated,'BoundingBox');
BB.xstart = ceil(stats.BoundingBox(1));
BB.ystart = ceil(stats.BoundingBox(2));
BB.xend = BB.xstart + stats.BoundingBox(3) - 1;
BB.yend = BB.ystart + stats.BoundingBox(4) - 1;

assert(BB.ystart < BB.yend);

refLineLength = round((BB.xend - BB.xstart)/2);
refLine.xstart = BB.xstart + round(refLineLength/2);
refLine.xend = BB.xend - round(refLineLength/2);
refLine.y = BB.yend;

translations.dx = nan(1,ntime);
translations.dy = nan(1,ntime);
translations.score = nan(1,ntime);

bbs = cell(1,ntime);
refLines = cell(1,ntime);

I0_orig = imread(inFname, tAnnotated, 'Info', info);
I0 = imfilter(I0_orig,params.gaussFilter);% See if this is helpful
I0Ref = I0(:,:,params.ch4Registration);
I0Ref = enhanceRef(I0Ref);

for i = 1 : ntime
    I1_orig = imread(inFname, i, 'Info', info);
    I1 = imfilter(I1_orig,params.gaussFilter); % See if this is helpful
    I1Ref = I1(:,:,params.ch4Registration);
    I1Ref = enhanceRef(I1Ref);
    [translations.dy(i),translations.dx(i),translations.score(i),newBB,newRefLine]...
        = calcTranslation(I0Ref,I1Ref,BB,refLine,params.searchR); 
    
    bbs{i} = newBB;
    refLines{i} = newRefLine;
end
end

%%
function [translations,bbs] = calcTranslations(inFname,info,params,ttime, roiAnnotated)
ntime = length(ttime);

translations.dx = nan(1,ntime);
translations.dy = nan(1,ntime);
translations.score = nan(1,ntime);
bbs = cell(1,ntime);

I0 = imread(inFname, ttime(1), 'Info', info);
stats = regionprops(roiAnnotated,'BoundingBox');
curBB.xstart = ceil(stats.BoundingBox(1));
curBB.ystart = ceil(stats.BoundingBox(2));
curBB.xend = curBB.xstart + stats.BoundingBox(3) - 1;
curBB.yend = curBB.ystart + stats.BoundingBox(4) - 1;
% curROI = false(size(roiAnnotated));
% curROI(BB(2):BB(2)+BB(4)-1,BB(1):BB(1)+BB(2)-1) = true;

translations.dx(1) = 0;
translations.dy(1) = 0;
translations.score(1) = 1;
bbs{1} = curBB;

for i = 2 : length(ttime)
    curT = ttime(i);
    I1 = imread(inFname, curT, 'Info', info);
    [translations.dy(i),translations.dx(i),translations.score(i),newBB]...
        = calcTranslation(I0(:,:,params.ch4Registration),I1(:,:,params.ch4Registration),curBB,params.searchR);    
    curBB = newBB;
    bbs{i} = curBB;
end
end

%%
function [dy,dx,score,newBB,newRefLine] = calcTranslation(im0,im1,BB,refLine,searchR)

im0 = double(im0);
im1 = double(im1);

itr = -searchR:searchR;
nItr = length(itr);

bb0 = im0(BB.ystart:BB.yend,BB.xstart:BB.xend);
bb0 = bb0 ./ sum(bb0(:));

allScores = nan(nItr);
for idy = itr
    for idx = itr
        try %% out of range...
            bb1 = im1((BB.ystart:BB.yend)+idy,(BB.xstart:BB.xend)+idx);
        catch ee
            continue;
        end
        bb1 = bb1 ./ sum(bb1(:));
        corr = sqrt(bb0 .* bb1); % Bhattacharyya coefficient
        allScores(idy+searchR+1,idx+searchR+1) = sum(corr(:));
    end
end

[maxVal,maxInd] = max(allScores(:));
[idy,idx] = ind2sub(size(allScores),maxInd);

dx = idx - 1 - searchR;
dy = idy - 1 - searchR;
score = maxVal;

newBB.xstart = BB.xstart + dx;
newBB.ystart = BB.ystart + dy;
newBB.xend = BB.xend + dx;
newBB.yend = BB.yend + dy;

newRefLine.xstart = refLine.xstart + dx;
newRefLine.xend = refLine.xend + dx;
newRefLine.y = refLine.y + dy;

end

%%
function [dys,dxs,scores] = integrateTrans(pastTranslations,futureTranslations)
dxPast = cumsum(pastTranslations.dx);
dyPast = cumsum(pastTranslations.dy);
dxFuture = cumsum(futureTranslations.dx);
dyFuture = cumsum(futureTranslations.dy);

dxs = [dxPast(end:-1:1),dxFuture(2:end)];
dys = [dyPast(end:-1:1),dyFuture(2:end)];
scores = [pastTranslations.score futureTranslations.score(2:end)];
end

%%
function [imgs,imgsTrans] = applyTranslations(inFname,info,dys,dxs,params)
ntime = length(dys);
imgs = cell(1,ntime);
imgsTrans = cell(1,ntime);
for i = 1 : ntime
    I_orig = imread(inFname, i, 'Info', info);
    I = imfilter(I_orig,params.gaussFilter);
    ITrans = imtranslate(I,[dys(i),dxs(i),0]);
    imgs{i} = I;
    imgsTrans{i} = ITrans;
end
end

%%
function [Iroi,Icrop,refLineData,refLineDataCrop] = getComposite(im,bb,refLine,params)

refLineData = optimizeLineRef(im,refLine,params);

sizeY = size(im,1);
sizeX = size(im,2);

Iroi = im;
tmp = zeros(sizeY,sizeX);
tmp(max(1,bb.ystart):min(sizeY,bb.yend),max(1,bb.xstart):min(sizeX,bb.xend)) = 50;
Iroi(:,:,3) = tmp;
Iroi(refLineData.y,refLineData.xstart:refLineData.xend,3) = 100;

refLineDataCrop.y = refLineData.y - bb.ystart + 1;
refLineDataCrop.xstart = refLineData.xstart - bb.xstart + 1;
refLineDataCrop.xend = refLineData.xend - bb.xstart + 1;

Icrop = im(max(1,bb.ystart):min(sizeY,bb.yend),max(1,bb.xstart):min(sizeX,bb.xend),:);
Icrop(refLineDataCrop.y,refLineDataCrop.xstart:refLineDataCrop.xend,3) = 150;
end


%%
function Iout = enhanceRef(Iin)
Iout = Iin;
% maxVal = max(Iin(:));
% p85 = prctile(Iin(:),85);
% Iout = Iin;
% Iout(Iout > p85) = maxVal;
end

%% search for reference line
function refLineData = optimizeLineRef(im,refLine,params)
maxIntensity = -inf;
maxYPos = -inf;
for iy = refLine.y : -1 : (refLine.y - params.dilateBB)
    curIntensity = median(im(iy,refLine.xstart:refLine.xend,params.ch4Registration));
    if curIntensity > maxIntensity
        maxIntensity = curIntensity;
       maxYPos = iy;
    end
end
refLineData.y = maxYPos;
refLineData.xstart = refLine.xstart; 
refLineData.xend = refLine.xend;
refLineData.dataCh1 = im(refLineData.y,refLine.xstart:refLine.xend,1);
refLineData.dataCh2 = im(refLineData.y,refLine.xstart:refLine.xend,2);
end

%%

% imwrite(I0, outFname);
% accDx = 0;
% accDy = 0;
% 
% % inds = 2 : n;
% % inds = 3 : 2 : n;
% inds = 4 : 2 : n;
% % for i = 2 : n
% % for i = 3 : 2 : n
% for i = 1 : length(inds)
%     curInd = inds(i);
%     I1 = imread(inFname, curInd, 'Info', info);        
%     translations{i+1} = calcTranslation(I0(:,:,params.ch4Registration),I1(:,:,params.ch4Registration),params.searchR);     
%     I1Trans = imtranslate(I1,-[translations{i+1}.dy, translations{i+1}.dx]); % old implementation
%     %     I1Trans = imtranslate(I1,[translations{i+1}.dx, translations{i+1}.dy]); % old implementation    
%     imwrite(I1Trans, outFname, 'writemode', 'append');    
%     disp([num2str(i+1) ': (dy,dx) = (' num2str(translations{i+1}.dy) ',' num2str(translations{i+1}.dx) ')']);
%     I0 = I1;
% end
% end

%%
function [tAnnotated,roiAnnotated] = findAnnotatedFrame(inFname,info,dilateR)
for it = 1 : numel(info)
    Icur = imread(inFname, it, 'Info', info);
    Iblue = Icur(:,:,3); 
    if sum(Iblue(:)) == 0
        continue;
    else
        tAnnotated = it;
        roiAnnotated = logical(Iblue);
        roiAnnotated = imfill(roiAnnotated,'holes');
        roiAnnotated = imdilate(roiAnnotated,strel('square',dilateR));
        return;
    end    
end
error('Could not find ROI in frame');    
end



% Todo: add caxis based on the video! (move this to meta analysis)
function [] = plotKymograph(dataCh,kymographFname)
h = figure;
hold on;
imagesc(dataCh);
xlabel('Space');
ylabel('Time');
colorbar;
hold off;
saveas(h,kymographFname);
end