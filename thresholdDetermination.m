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
% commit: threshold determination for 2019-02-06

% OK let's go!

%% A. initialize data and input manually selected cell IDs

clc
clear

% 0. initialize data
xy = 1;
date = '2019-02-06';
cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date))
load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');


% 0. initialize visually scored (+) and (-) fluor signal particle IDs

% from image 80 in xy01 (dt = 3 min)
yfp_plus_4h = [11,265,308,285,66,88,163,267,148,76,284,264,330,265,147,294,332,18]; % yfp+ cells 
yfp_maybes_4h = [179,313,329,302,445,135,47]; % ambiguous yfp signal (very light signal)
cfp_plus_4h = [196,24,242,20,250,91,60,65,258,191,224,96,259,77,129,94,236,153,74,73,455,116,234,218]; % cfp+

yfp_neg_4h = cfp_plus_4h;
cfp_neg_4h = yfp_plus_4h; % cfp-

% yfp_neg_4h = [202,225,92,128,117,398,205,105,39,6,16,114]; % yfp- cells


% from image 161 in xy01 
yfp_plus_8h = [11,265,308,285,88,38,490,88,38,490,267,148,284,330,147,294,332,18,287,62,150];
yfp_maybes_8h = [163,277,172,349,72,50,30,179];
cfp_plus_8h = [196,24,242,20,91,60,65,258,191,224,96,259,77,129,94,236,153,74,73,218,116,234,145,204,122,71];

yfp_neg_8h = cfp_plus_8h;
cfp_neg_8h = yfp_plus_8h;

% yfp_neg_8h = [60,65,91,202,90,225,205,364,71,143,224,168];
% cfp_neg_8h = [116,128,205,50,105,195,247,165,216,293,232,14];


% double-check data entry

% 1. should not have co-occuring particles between sets of same fluorophore
check1 = intersect(yfp_plus_8h,yfp_neg_8h);
check2 = intersect(yfp_plus_8h,yfp_maybes_8h);
check3 = intersect(yfp_neg_8h,yfp_maybes_8h);

% 2. should not have co-occuring particles between fluorophores
check4 = intersect(yfp_plus_8h,cfp_plus_8h);
check5 = intersect(yfp_maybes_8h,cfp_plus_8h);


% 3. report: if checks fail, throw error message
check_summary = check1 + check2 + check3 + check4 + check5;
if isempty(check_summary) == 1
    disp('-  quality report: fluorescent data checkpoints passed!')
    clear check1 check2 check3 check4 check5 check_summary
end



% 1. compile experiment data matrix
dt_min = 3;
xy_start = xy;
xy_end = xy;
xyData = buildDM_glycogen(D5, xy_start, xy_end,dt_min);
clear xy_start xy_end


%% B. (i) YFP-based particle sets, image 80

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp plus neg maybes idx_plus idx_neg idx_maybes
clc

% 1. isolate image80 data
image_interest = 80;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
concat_yfp = [yfp_plus_4h'; yfp_neg_4h'; yfp_maybes_4h'];
indeces = zeros(length(concat_yfp),1);
for i = 1:length(indeces)
    idx = find(particleID == concat_yfp(i));
    indeces(i) = idx;
end
clear i idx


% 4. resort indeces into cell sets
plus = length(yfp_plus_4h);
neg = length(yfp_neg_4h);

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
title('YFP signal in in YFP-based cell sets, xy01 at 4h')
%xticklabels({'YFP+','YFP-','YFP maybes'}) % when updating matlab versions

% CFP signal intensity from YFP-based sets
figure(2)
bar(meanCFP)
hold on
errorbar(1:3,meanCFP,stdCFP,'.')
ylim([95 165])
ylabel('Mean CFP Intensity')
title('CFP signal in YFP-based cell sets, xy01 at 4h')

%% B. (ii) CFP-based particle sets, image 80

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp plus neg maybes idx_plus idx_neg idx_maybes
clc

% 1. isolate image data
image_interest = 80;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
concat_cfp = [cfp_plus_4h'; cfp_neg_4h'];
indeces = zeros(length(concat_cfp),1);
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 5. resort indeces into cell sets
plus = length(cfp_plus_4h);
neg = length(cfp_neg_4h);

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
title('YFP signal in cfp-based cells, xy01 at 4h')

figure(4)
bar(meanCFP)
hold on
errorbar(1:2,meanCFP,stdCFP,'.')
ylim([95 165])
ylabel('Mean CFP Intensity')
title('CFP signal in cfp-based cells, xy01 at 4h')

%% B. (iii) YFP-based particle sets, image 161

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp plus neg maybes idx_plus idx_neg idx_maybes
clc


% 1. isolate image data
image_interest = 161;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity


% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
concat_yfp = [yfp_plus_8h'; yfp_neg_8h'; yfp_maybes_8h'];
indeces = zeros(length(concat_yfp),1);
for i = 1:length(indeces)
    idx = find(particleID == concat_yfp(i));
    indeces(i) = idx;
end
clear i idx


% 4. resort indeces into cell sets
plus = length(yfp_plus_8h);
neg = length(yfp_neg_8h);

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
title('YFP signal in YFP-based cell sets, xy01 at 8h')
%xticklabels({'YFP+','YFP-','YFP maybes'}) % when updating matlab versions

figure(6)
bar(meanCFP)
hold on
errorbar(1:3,meanCFP,stdCFP,'.')
ylim([95 165])
ylabel('Mean CFP Intensity')
title('CFP signal in YFP-based cell sets, xy01 at 8h')

%% B. (iv) CFP-based particle sets, 8h

clear meanYFP meanCFP stdYFP stdCFP yInt cInt image_interest frame dm_image
clear int_yfp int_cfp plus neg maybes idx_plus idx_neg idx_maybes
clc


% 1. isolate image data
image_interest = 161;
frame = xyData(:,9);        % col 9 = frame number in buildDM_glycogen
dm_image = xyData(frame == image_interest,:);


% 2. for all particles in yfp+, gather yfp and cfp intensities
int_yfp = dm_image(:,14);    % col 14 = mean yfp intensity
int_cfp = dm_image(:,13);    % col 13 = mean cfp intensity



% 3. identify indeces of each particle
particleID = dm_image(:,15);  % col 15 = particle ID from tracking
concat_cfp = [cfp_plus_8h'; cfp_neg_8h'];
indeces = zeros(length(concat_cfp),1);
for i = 1:length(indeces)
    idx = find(particleID == concat_cfp(i));
    indeces(i) = idx;
end
clear i idx


% 5. resort indeces into cell sets
plus = length(cfp_plus_8h);
neg = length(cfp_neg_8h);

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
title('YFP signal in cfp-based cells, xy01 at 8h')

figure(8)
bar(meanCFP)
hold on
errorbar(1:2,meanCFP,stdCFP,'.')
ylim([95 165])
ylabel('Mean CFP Intensity')
title('CFP signal in cfp-based cells, xy01 at 8h')

%% C. save figures

for f = 1:8
    figure(f)
    plotName = strcat('bar-figure',num2str(f));
    saveas(gcf,plotName,'epsc')
    close(gcf)
end


