function [] = tifStackZeroToNan(tifStackFnameIn)

if nargin < 1
    tifStackFnameIn = 'C:\Users\CellBio1\Google Drive\Research\PostDoc\Collaborations\Dagan\zero2nan\c8v2_1.tif';
end
    
tifStackFnameOut = [tifStackFnameIn(1:end-4) '_nan.tif'];

close all;

info = imfinfo(tifStackFnameIn);
nTime = numel(info);

for i = 1 : nTime    
    Icur = imread(tifStackFnameIn, i, 'Info', info);
    Icur(Icur == 0) = nan;
    
    if i == 1
        imwrite(Icur,tifStackFnameOut);        
    else
        imwrite(Icur,tifStackFnameOut, 'writemode', 'append');        
    end    
end
end