%% DaganMain
%   Quantifying molecular dynamics at non-contracting vesicles in the
%   Drosophila Salivry gland
%
%   Please cite the following manuscript if using this code: 
%       Formation and disassembly of a contractile actomyosin network mediates content release from large secretory vesicles
%       https://doi.org/10.1101/216044 (accepted to JCB)
%
%   Input: 
%       dname - main directory path
%       baseDname - condition (e.g., rockinh - ROCK inhibition)
%       tifStacksDnames - cell array of independent experiments / replicates (directories)
%       params.timePerFrame - time resolution for each experiment ()
%       dbFname - excel file name 
%       params - additional parameters for the analysis
%
%       Assumptions: 
%           Each vesicle has a code civj, cell i and vesicle j
%           Each vesicle is represented by (1) data file, 2 data channels,
%           the 3rd channel for defining the reference frame (which defines the bounding box for registration)
%
%       Processing:
%           DaganAnalyseExperiments - analyze vesicles (registration and intensity measurements)
%           DaganParseVesicleXls -  parse xls file (parsing and partitioning to groups)
%           DaganMetaBackgroundAnalysis - calculate background model (time points before vesicle formation, to normalize per experiment + bleaching correction)
%           DaganMetaNormalizeAnalysis - normalized vesicle data based on
%               background model, find oscillations (via DaganFindOscillations)
%           DaganMetaAnalysis - accumulate meta data per experiment (including sustained vesicles analysis)
%
%   Assaf Zaritsky, 2017 - for Dagan Segal's project
function [] = DaganMain()

close all;
clc;

%% experiments
dname = 'C:\Users\assafza\Google Drive\Research\PostDoc\Collaborations\Dagan\';

baseDname = [dname 'testData\'];
tifStacksDnames = {'1rockinh','2rockinh'};
params.timePerFrame = [0.469,0.348];
dbFname = 'metaData.xlsx';

% baseDname = [dname 'rockinh\'];
% tifStacksDnames = {'1rockinh','2rockinh','7rockinh','8rockinh'};
% params.timePerFrame = [0.469,0.348,0.339,0.348];
% dbFname = 'vesicles_timepoints_20170622.xlsx';

% baseDname = [dname 'rockinh_arpinh\'];
% tifStacksDnames = {'3rockinh_arpinh','4rockinh_arpinh','5rockinh_arpinh'};
% % tifStacksDnames = {'4rockinh_arpinh'};
% dbFname = 'vesicles timepoints_arpinh_20170626.xlsx'; % 23
% params.timePerFrame = [0.352,0.352,0.352];

% baseDname = [dname 'rhogapi_rockinh\'];
% tifStacksDnames = {'1rhogapi_rockinh','6rhogapi_rockinh','8rhogapi_rockinh'};
% params.timePerFrame = [0.355,0.335,0.335];
% % dbFname = 'vesicles_timepoints_rhogapi_20170815.xlsx';
% dbFname = 'vesicles_timepoints_rhogapi_20170918_corrected.xlsx';

% baseDname = [dname 'Rho_hightimeres\'];
% % tifStacksDnames = {'2RhoGlaR_short','3RhoGlaR_short','4RhoGlaR_short'};
% tifStacksDnames = {'3RhoGlaR_short','4RhoGlaR_short'};
% params.timePerFrame = [2/60.0,2/60.0,2/60.0];
% % dbFname = 'rhoGlaRHighTimeres.xlsx';
% dbFname = 'rhoGlaRHighTimeres26092017_partial.xlsx';

% baseDname = [dname 'RhoGAP_hightimeres\'];
% tifStacksDnames = {'1_3rhoGAPGFPshort','1rhoGAPGFPshort','3rhoGAPGFPshort'};
% params.timePerFrame = [4.55/60.0,4.55/60.0,4.55/60.0];
% % dbFname = 'rhogapGFPlaRubyhightimeres.xlsx';
% dbFname = 'rhogapGFPlaRubyhightimeres26092017.xlsx';



%% parameters

% experiments parameetrs
params.ch4Registration = 1; % 1 - Rho, 2 - Actin

% vesicle registration and qunatification
params.always = true;
params.searchR = 10;
params.dilateR = 3;
params.dilateBB = 5;
hsize = [5,5]; % default: [3,3]
sigma = 0.7; % default: 0.5
params.gaussFilter = fspecial('gaussian',hsize,sigma); % fspecial('gaussian', hsize, sigma)

% for 1D smoothing of vesicle signal over time
params.smoothWinSize = 5;
params.smoothSigma = 2.5;

% oscillations
params.MinPeakHeight = 3;% 2; % peaks must be 3 std above background
params.MinPeakProminence = 1.5; % 1; % drop off on both sides before encountering a larger value
params.MinPeakDistance = 7; % Minimum peak separation

%% flags
flag.doAnalyzeExp = true;
flag.doParseVesicleXls = true;
flag.doMetaBackgroundAnalysis = true;
flag.doMetaNormalizeAnalysis = true;
flag.doMetaAnalysis = true;

%% Execute
if flag.doAnalyzeExp
    DaganAnalyseExperiments(baseDname,tifStacksDnames,params); % local analysis and identify all missing ROIs
end

if flag.doParseVesicleXls
    DaganParseVesicleXls(baseDname,dbFname); % parse xls file
end

dbFname = [dbFname(1:end-5) '.mat'];

if flag.doMetaBackgroundAnalysis
    DaganMetaBackgroundAnalysis(baseDname,dbFname); % calculate background statistics
end

if flag.doMetaNormalizeAnalysis
    DaganMetaNormalizeAnalysis(baseDname,dbFname,params); % normalized data, find oscillations
end

if flag.doMetaAnalysis
    DaganMetaAnalysis(baseDname,params.timePerFrame); % meta analysis, sustained vesicles
end

end