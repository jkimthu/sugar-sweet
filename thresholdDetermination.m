% threshold determination - for isolating CFP and YFP populations


% Goal: manually determine threshold intensity value, above which to define
%       a cell as positive for a fluorophore.

%       because each experiment has different absolute intensities, due to
%       different channels, PDMS thickness, cleanliness, etc., this
%       determination needs to be checked for each experiment.



% Input data is the QC-trimmed .mat file from particle tracking of raw timelapse
% videos with three channels: phase, YFP and CFP.



% Strategy:

% (pre-code, manual)
%
%                i. make a list of YFP+ and CFP+ cells
%               ii. confirm that i can visually detect only one fluorophore in each cell
%              iii. do this for one xy position at 4h and 8h into experiment


% (in following code)
%
%                A. initialize data and input manually selected cell IDs
%                B. calculate mean and stdev of signal intensities from:
%                          i. YFP+ cells from 4h image
%                         ii. CFP+ cells from 4h image
%                        iii. YFP+ cells from 8h image
%                         iv. CFP+ cells from 8h image
%                C. plot and save bar plots of summary statistics


% last updated: jen, 2019 Feb 11
% commit: first commit, threshold determination for 2018-11-23

% OK let's go!

%% A. initialize data and input manually selected cell IDs

clc
clear

% 0. initialize data
xy = 2;
date = '2018-11-23';
cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date,'_xy02'))
load(strcat('glycogen-',date,'-earlyEdits-jiggle-0p5.mat'),'D5');


% 0 . initialize image frequency
selected_tpt = linspace(1,241,17);


% 0. initialize visually scored (+) and (-) fluor signal particle IDs

% from image 9 in xy02 reduced frequency dataset
yfp_plus_9 = [25,28,27,7,18,13,3,1,24,31,42,21,31,9,2,40,35,45,28,33,29,26,32]; % yfp+ cells 
yfp_neg_9 = [36,20,15,4,14,23,16,39,37,30,6,10,19,22,34,5]; % yfp- cells
yfp_maybes_9 = [46,11,8,47,44,12,43,41]; % ambiguous yfp signal (very light signal)

cfp_plus_9 = 4;


% from image 17 in xy02 reduced frequency dataset
yfp_plus_17 = [24,25,27,18,13,7,3,1,31,21,17,9,40,35,28,33,29,26,32];
yfp_neg_17 = [39,38,20,15,4,14,23,16,37,30,6,2,10,19,22,5];
yfp_maybes_17 = [11,8,12,41];

cfp_plus_17 = 4;


% double-check data entry

% 1. should not have co-occuring particles between sets of same fluorophore
check1 = intersect(yfp_plus_17,yfp_neg_17);
check2 = intersect(yfp_plus_17,yfp_maybes_17);
check3 = intersect(yfp_neg_17,yfp_maybes_17);

% 2. should not have co-occuring particles between fluorophores
check4 = intersect(yfp_plus_17,cfp_plus_17);
check5 = intersect(yfp_maybes_17,cfp_plus_17);

% 3. SHOULD have co-occurrence between yfp- and cfp+, represeting 100% of cfp+
check6 = intersect(yfp_neg_17,cfp_plus_17);
check7 = cfp_plus_17 == check6;

% 4. report: if checks fail, throw error message
check_summary = check1 + check2 + check3 + check4 + check5 + check6;
if isempty(check_summary) == 0
    error('checkpoint failure! double-check fluorescent subsets')
elseif check7 == 0
    error('checkpoint alert! not all CFP+ cells are YFP-')
else
    disp('-  quality report: fluorescent data checkpoints passed!')
    clear check1 check2 check3 check4 check5 check6 check7 check_summary
end



% 1. compile experiment data matrix
dt_min = 2;
xy_start = xy;
xy_end = xy;
xyData = buildDM_glycogen(D5, xy_start, xy_end,dt_min);
clear xy_start xy_end



% 2. gather fluorescence intensities for all particle sets

%% B. (i) YFP-based particle sets, image 9

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp cfp_neg plus neg maybes idx_plus idx_neg idx_maybes
clc

% 1. isolate image9 data
image_interest = 9;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
indeces = zeros(length(particleID),1);
concat_cfp = [yfp_plus_9'; yfp_neg_9'; yfp_maybes_9'];
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 4. resort indeces into cell sets
plus = length(yfp_plus_9);
neg = length(yfp_neg_9);

idx_plus = indeces(1:plus);
idx_neg = indeces(plus+1:plus+neg);
idx_maybes = indeces(plus+neg+1:end);


% 5. isolate intensity based on cell set
yInt{1} = int_yfp(idx_plus);
yInt{2} = int_yfp(idx_neg);
yInt{3} = int_yfp(idx_maybes);

cInt{1} = int_cfp(idx_plus);
cInt{2} = int_cfp(idx_neg);
cInt{3} = int_cfp(idx_maybes);


% 6. calculate mean and standard devation of cell intensity
meanYFP = cellfun(@mean,yInt);
meanCFP = cellfun(@mean,cInt);
stdYFP = cellfun(@std,yInt);
stdCFP = cellfun(@std,cInt);


% 7. plot
figure(1)
bar(meanYFP)
hold on
errorbar(1:3,meanYFP,stdYFP,'.')
ylim([95 130])
ylabel('Mean YFP Intensity')
title('YFP signal in in YFP-based cell sets, xy02 reduced freq image 9')
%xticklabels({'YFP+','YFP-','YFP maybes'}) % when updating matlab versions

% CFP signal intensity from YFP-based sets
figure(2)
bar(meanCFP)
hold on
errorbar(1:3,meanCFP,stdCFP,'.')
ylim([95 130])
ylabel('Mean CFP Intensity')
title('CFP signal in YFP-based cell sets, in xy02 reduced freq image 9')

%% B. (ii) CFP-based particle sets, image 9

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp cfp_neg plus neg maybes idx_plus idx_neg idx_maybes
clc

% 1. isolate image data
image_interest = 9;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. determine CFP- cells based on CFP+
particleIDs = dm_image(:,12);    % col 12 = track number, which for single xy dm is also ID
cfp_neg = particleIDs;
cfp_neg(cfp_plus_9) = [];


% 4. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
indeces = zeros(length(particleID),1);
concat_cfp = [cfp_plus_9'; cfp_neg];
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 5. resort indeces into cell sets
plus = length(cfp_plus_9);
neg = length(cfp_neg);

idx_plus = indeces(1:plus);
idx_neg = indeces(plus+1:plus+neg);


% 6. isolate intensity based on cell set
yInt{1} = int_yfp(idx_plus);
yInt{2} = int_yfp(idx_neg);

cInt{1} = int_cfp(idx_plus);
cInt{2} = int_cfp(idx_neg);


% 7. calculate mean and standard deviation of cell intensity
meanYFP = cellfun(@mean,yInt);
meanCFP = cellfun(@mean,cInt);
stdYFP = cellfun(@std,yInt);
stdCFP = cellfun(@std,cInt);


% 8. plot
figure(3)
bar(meanYFP)
hold on
errorbar(1:2,meanYFP,stdYFP,'.')
ylim([95 130])
ylabel('Mean YFP Intensity')
title('YFP signal in cfp-based cells, xy02 reduced freq image 9')

figure(4)
bar(meanCFP)
hold on
errorbar(1:2,meanCFP,stdCFP,'.')
ylim([95 130])
ylabel('Mean CFP Intensity')
title('CFP signal in cfp-based cells, xy02 reduced freq image 9')

%% B. (iii) YFP-based particle sets, image 17

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp cfp_neg plus neg maybes idx_plus idx_neg idx_maybes
clc


% 1. isolate image data
image_interest = 17;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
indeces = zeros(length(particleID),1);
concat_cfp = [yfp_plus_17'; yfp_neg_17'; yfp_maybes_17'];
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 4. resort indeces into cell sets
plus = length(yfp_plus_17);
neg = length(yfp_neg_17);

idx_plus = indeces(1:plus);
idx_neg = indeces(plus+1:plus+neg);
idx_maybes = indeces(plus+neg+1:end);


% 5. isolate intensity based on cell set
yInt{1} = int_yfp(idx_plus);
yInt{2} = int_yfp(idx_neg);
yInt{3} = int_yfp(idx_maybes);

cInt{1} = int_cfp(idx_plus);
cInt{2} = int_cfp(idx_neg);
cInt{3} = int_cfp(idx_maybes);


% 6. calculate mean and standard devation of cell intensity
meanYFP = cellfun(@mean,yInt);
meanCFP = cellfun(@mean,cInt);
stdYFP = cellfun(@std,yInt);
stdCFP = cellfun(@std,cInt);


% 7. plot
figure(5) 
bar(meanYFP)
hold on
errorbar(1:3,meanYFP,stdYFP,'.')
ylim([95 130])
ylabel('Mean YFP Intensity')
title('YFP signal in YFP-based cell sets, xy02 reduced freq image 17')
%xticklabels({'YFP+','YFP-','YFP maybes'}) % when updating matlab versions

figure(6)
bar(meanCFP)
hold on
errorbar(1:3,meanCFP,stdCFP,'.')
ylim([95 130])
ylabel('Mean CFP Intensity')
title('CFP signal in YFP-based cell sets, xy02 reduced freq image 17')

%% B. (iv) CFP-based particle sets, image 17

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp cfp_neg plus neg maybes idx_plus idx_neg idx_maybes
clc


% 1. isolate image data
image_interest = 17;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. determine CFP- cells based on CFP+
particleIDs = dm_image(:,12);    % col 12 = track number, which for single xy dm is also ID
cfp_neg = particleIDs;
cfp_neg(cfp_plus_17) = [];


% 4. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
indeces = zeros(length(particleID),1);

concat_cfp = [cfp_plus_17'; cfp_neg];
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 5. resort indeces into cell sets
plus = length(cfp_plus_17);
neg = length(cfp_neg);

idx_plus = indeces(1:plus);
idx_neg = indeces(plus+1:plus+neg);



% 6. isolate intensity based on cell set
yInt{1} = int_yfp(idx_plus);
yInt{2} = int_yfp(idx_neg);

cInt{1} = int_cfp(idx_plus);
cInt{2} = int_cfp(idx_neg);


% 7. calculate mean and standard deviation of cell intensity
meanYFP = cellfun(@mean,yInt);
meanCFP = cellfun(@mean,cInt);
stdYFP = cellfun(@std,yInt);
stdCFP = cellfun(@std,cInt);


% 8. plot
figure(7)
bar(meanYFP)
hold on
errorbar(1:2,meanYFP,stdYFP,'.')
ylim([95 130])
ylabel('Mean YFP Intensity')
title('YFP signal in cfp-based cells, xy02 reduced freq image 17')

figure(8)
bar(meanCFP)
hold on
errorbar(1:2,meanCFP,stdCFP,'.')
ylim([95 130])
ylabel('Mean CFP Intensity')
title('CFP signal in cfp-based cells, xy02 reduced freq image 17')

%% C. save figures

for f = 1:8
    figure(f)
    plotName = strcat('bar-figure',num2str(f));
    saveas(gcf,plotName,'epsc')
    close(gcf)
end


