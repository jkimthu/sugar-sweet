%% quantifyGrowth

% Goal: quantify and visualize growth rates of WT and mutants for glycogen
%       utilization under a pulsling nutrient environment.
%
%       input data is data structure produced from glycogen_analysis.m

%       this script contains methods of measuring growth rate:
%
%           section 1. instantaneous rate of volumetric doubling
%           section 2. total area occupied, normalized per cell and initial value
%           section 3. dA/dt, where A is normalized by the number of cells
%           - for specifics, see strategy at the start of each section


% last update: Jen, 2019 Feb 22
% commit: analysis for 2019-02-19 experiment, steady control

% ok let's go!

%% section 1. instantaneous rate of volumetric doubling


% Strategy:
%                0. initialize experiment parameters and data
%                1. assemble data matrix
%                2. isolate volume (Va), drop, and track number for growth rate calculations
%                3. calculate growth rate
%                4. truncate data to non-erroneous timestamps (e.g. bubbles) 
%                5. bin growth rate into time bins based on timestamp
%                6. isolate selected specific growth rate and remove nans from data analysis
%                7. isolate YFP and CFP intensities
%                8. convert intensities to (+) or (-) fluorophore
%                       - for details on threshold determination,
%                         see scrap_pile section "YFP and CFP thresholds"
%                9. test threshold, throw error if any data points ID as both
%               10. separate growth rates by fluorophore
%               11. bin growth rate by time
%               12. calculate stats for each bin, including mean growth rate
%               13. plot growth rate over time



clear
clc

% 0. initialize data
%xy = 2;
xy_start = 25;
xy_end = 29;
dt_min = 3;
%dt_min = 30; % reduced frequency dataset

date = '2019-02-19';

%cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date))
load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');


% 0. define growth rates of interest
specificGrowthRate = 'log2';


% 0. define time binning parameters
specificBinning = 60; % in minutes
binsPerHour = 60/specificBinning;


% 0. define fluorescence intensity threshold
threshold = 103.4;


% 1. assemble data matrix
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
clear xy_start xy_end


% 2. isolate volume (Va), drop, and track number for growth rate calculations
volumes = glycogen_data(:,5);        % col 4 = calculated va_vals (cubic um)
isDrop = glycogen_data(:,3);         % col 2 = isDrop, 1 marks a birth event
trackNum = glycogen_data(:,12);      % col 12 = track number (not ID from particle tracking)


% 3. calculate growth rate
dt_sec = dt_min * 60;
growthRates = calculateGrowthRate_glycogen(volumes,isDrop,trackNum,dt_sec);
clear isDrop volumes trackNum dt_min


% 4. truncate data to non-erroneous timestamps (e.g. bubbles) 
maxTime = 8; % in hours
frame = glycogen_data(:,9);      % col 9 = frame in image sequence
timeInSeconds = frame * dt_sec;  % frame = is consequetive images in analysis
timeInHours = timeInSeconds/3600;

if maxTime > 0
    glycogenData_bubbleTrimmed = glycogen_data(timeInHours <= maxTime,:);
    growthRates_bubbleTrimmed = growthRates(timeInHours <= maxTime,:);
else
    glycogenData_bubbleTrimmed = glycogen_data;
    growthRates_bubbleTrimmed = growthRates;
end
clear maxTime frame timeInSeconds timeInHours


% 5. bin growth rate into time bins based on timestamp
frame = glycogenData_bubbleTrimmed(:,9);      % col 9 = frame in image sequence
timeInSeconds = frame * dt_sec;    
timeInHours = timeInSeconds/3600;
bins = ceil(timeInHours*binsPerHour);
clear timeInSeconds frame


% 6. isolate selected specific growth rate and remove nans from data analysis
specificColumn = 3;
growthRate_log2 = growthRates_bubbleTrimmed(:,specificColumn);

growthRt_noNaNs = growthRate_log2(~isnan(growthRate_log2),:);
bins_noNaNs = bins(~isnan(growthRate_log2),:);
glycogenData_noNaNs = glycogenData_bubbleTrimmed(~isnan(growthRate_log2),:);


% 7. isolate YFP and CFP intensities
cfp = glycogenData_noNaNs(:,13);         % col 13 = mean CFP intensity
yfp = glycogenData_noNaNs(:,14);         % col 14 = mean YFP intensity


% 8. convert intensities to (+) or (-) fluorophore
isCFP = cfp > threshold;
isYFP = yfp > threshold;



% 9. test threshold, throw out any tracks that at any point ID as both
isBoth = isCFP+isYFP;
glycogenData_noDoublePos = glycogenData_noNaNs;
growthRt_noDoublePos = growthRt_noNaNs;
bins_noDoublePos = bins_noNaNs;

if sum(isBoth == 2) > 0
    
    disp('removing cells positive for both fluorophores')
    
    % identify trackNum of double-positives (errors)
    idx_errors = find(isBoth == 2);
    trackNum = glycogenData_noNaNs(:,12);   % col 12 = trackNum
    ids_errors = trackNum(idx_errors);
    unique_errors = unique(ids_errors);
    unique_errors_sorted = sort(unique_errors(:,1),'descend');
    
    
    % trim based on trackNum, in reverse order to maintain indeces
    for i = 1:length(unique_errors_sorted)
        
        currentTrack = unique_errors_sorted(i);
        growthRt_noDoublePos(glycogenData_noDoublePos(:,12) == currentTrack) = [];
        bins_noDoublePos(glycogenData_noDoublePos(:,12) == currentTrack) = [];
        glycogenData_noDoublePos(glycogenData_noDoublePos(:,12) == currentTrack,:) = [];
        
    end
        
end


% 10. finalize data in prep for sorting by fluorophore
growthRt_final = growthRt_noDoublePos;
glycogenData_final = glycogenData_noDoublePos;
bins_final = bins_noDoublePos;
isCFP_final = glycogenData_final(:,13) > threshold;   % col 13 = mean CFP intensity
isYFP_final = glycogenData_final(:,14) > threshold;   % col 14 = mean YFP intensity;



% 11. separate growth rates by fluorophore

% growthRt_yfp = growthRt_final(isYFP_final == 1);
% growthRt_cfp = growthRt_final(isCFP_final == 1);

% above method only takes growth rates if the cell is actively labeled.
% however, all the growth rates should be taken into account! instead, pull
% out IDs that are identified with YFP or CFP and use these to 
trackNum = glycogenData_final(:,12); % col 12 = track num
IDs_yfp = unique(trackNum(isYFP_final == 1));
IDs_cfp = unique(trackNum(isCFP_final == 1));

IDs_labeled = [IDs_yfp; IDs_cfp];

growthRt_yfp = [];
growthRt_cfp = [];
bins_yfp = [];
bins_cfp = [];
for id = 1:length(IDs_labeled)
    
    currentID = IDs_labeled(id);
    currentGR = growthRt_final(trackNum == currentID);
    currentBins = bins_final(trackNum == currentID);
    
    if ismember(currentID,IDs_yfp) == 1 % if ID belongs to YFP+
        growthRt_yfp = [growthRt_yfp; currentGR];
        bins_yfp = [bins_yfp; currentBins];
    else % else
        growthRt_cfp = [growthRt_cfp; currentGR];
        bins_cfp = [bins_cfp; currentBins];
    end
    
end
clear currentID currentGR currentBins


% 12. bin growth rate by time
binned_yfp = accumarray(bins_yfp,growthRt_yfp,[],@(x) {x});
binned_cfp = accumarray(bins_cfp,growthRt_cfp,[],@(x) {x});


% 13. calculate mean, standard dev, counts, and standard error
y_bin_means = cellfun(@mean,binned_yfp);
y_bin_stds = cellfun(@std,binned_yfp);
y_bin_counts = cellfun(@length,binned_yfp);
y_bin_sems = y_bin_stds./sqrt(y_bin_counts);

c_bin_means = cellfun(@mean,binned_cfp);
c_bin_stds = cellfun(@std,binned_cfp);
c_bin_counts = cellfun(@length,binned_cfp);
c_bin_sems = c_bin_stds./sqrt(c_bin_counts);



% 14. plot growth rate over time
palette = {'DodgerBlue','GoldenRod'};

yfp_color = rgb(palette(2));
cfp_color = rgb(palette(1));
xmark = 'o';

figure(1)
errorbar((1:length(y_bin_means))/binsPerHour,y_bin_means,y_bin_sems,'Color',yfp_color)
hold on
errorbar((1:length(c_bin_means))/binsPerHour,c_bin_means,c_bin_sems,'Color',cfp_color)
hold on
grid on
axis([0,8.5,-0.05,0.35])
xlabel('Time (hr)')
ylabel('Growth rate (1/hr)')
title(strcat(date,': (',specificGrowthRate,')'))
legend('YFP WT', 'CFP mutant')


%% section 2. Total Area Occupied, normalized per cell and initial value

% Strategy:
%                0. initialize experiment parameters and data
%                1. assemble data matrix
%                2. truncate data to non-erroneous timestamps (e.g. avoid bubbles) 
%                3. bin growth rate into time bins based on timestamp
%                4. isolate YFP and CFP intensities and area
%                5. convert intensities to (+) or (-) fluorophore
%                       - for details on threshold determination,
%                         see scrap_pile section "YFP and CFP thresholds"
%                6. test threshold, throw error if any data points ID as both
%                7. separate area by fluorophore
%                8. bin area by time
%                9. calculate stats for each bin, including summed area
%               10. normalize sums by cell counts
%               11. normalize by initial sum
%               12. plot growth rate over time


clear
clc

% 0. initialize data
xy_start = 1;
xy_end = 24;
dt_min = 3;
%dt_min = 30; % reduced frequency dataset

date = '2019-02-19';
%cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date))
load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');


% 0. define time binning parameters
specificBinning = 30; % in minutes
binsPerHour = 60/specificBinning;


% 0. define fluorescence intensity threshold
threshold = 103.4;


% 1. assemble data matrix
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
clear xy_start xy_end


% 2. truncate data to non-erroneous timestamps (e.g. bubbles) 
maxTime = 8; % in hours
frame = glycogen_data(:,9);      % col 9 = frame in image sequence
dt_sec = dt_min * 60;
timeInSeconds = frame * dt_sec;  % frame = is consequetive images in analysis
timeInHours = timeInSeconds/3600;

if maxTime > 0
    glycogenData_bubbleTrimmed = glycogen_data(timeInHours <= maxTime,:);
else
    glycogenData_bubbleTrimmed = glycogen_data;
end
clear maxTime frame timeInSeconds timeInHours



% 3. bin growth rate into time bins based on timestamp
frame = glycogenData_bubbleTrimmed(:,9);      % col 9 = frame in image sequence
timeInSeconds = frame * dt_sec;
timeInHours = timeInSeconds/3600;
bins = ceil(timeInHours*binsPerHour);
clear timeInSeconds frame



% 4. isolate YFP and CFP intensities and area
cfp = glycogenData_bubbleTrimmed(:,13);         % col 13 = mean CFP intensity
yfp = glycogenData_bubbleTrimmed(:,14);         % col 14 = mean YFP intensity
%area = glycogenData_bubbleTrimmed(:,6);         % col 6 = measured surface area



% 5. convert intensities to (+) or (-) fluorophore
isCFP = cfp > threshold;
isYFP = yfp > threshold;


% 6. test threshold, throw out any tracks that at any point ID as both
isBoth = isCFP+isYFP;
glycogenData_noDoublePos = glycogenData_bubbleTrimmed;
bins_noDoublePos = bins;

if sum(isBoth == 2) > 0
    
    disp('removing cells positive for both fluorophores')
    
    % identify trackNum of double-positives (errors)
    idx_errors = find(isBoth == 2);
    trackNum = glycogenData_bubbleTrimmed(:,12);   % col 12 = trackNum
    ids_errors = trackNum(idx_errors);
    unique_errors = unique(ids_errors);
    unique_errors_sorted = sort(unique_errors(:,1),'descend');
    
    
    % trim based on trackNum, in reverse order to maintain indeces
    for i = 1:length(unique_errors_sorted)
        
        currentTrack = unique_errors_sorted(i);
        glycogenData_noDoublePos(glycogenData_noDoublePos(:,12) == currentTrack,:) = [];
        
    end
        
end


% 7. finalize data in prep for sorting by fluorophore
glycogenData_final = glycogenData_noDoublePos;
bins_final = bins_noDoublePos;
isCFP_final = glycogenData_final(:,13) > threshold;   % col 13 = mean CFP intensity
isYFP_final = glycogenData_final(:,14) > threshold;   % col 14 = mean YFP intensity;
area_final = glycogenData_final(:,6);                 % col 6 = measured surface area


% 8. separate area by fluorophore
trackNum = glycogenData_final(:,12); % col 12 = track num
IDs_yfp = unique(trackNum(isYFP_final == 1));
IDs_cfp = unique(trackNum(isCFP_final == 1));

IDs_labeled = [IDs_yfp; IDs_cfp];

area_yfp = [];
area_cfp = [];
time_yfp = [];
time_cfp = [];
for id = 1:length(IDs_labeled)
    
    currentID = IDs_labeled(id);
    currentSA = area_final(trackNum == currentID);
    currentTime = bins_final(trackNum == currentID);
    
    if ismember(currentID,IDs_yfp) == 1 % if ID belongs to YFP+
        area_yfp = [area_yfp; currentSA];
        time_yfp = [time_yfp; currentTime];
    else % else
        area_cfp = [area_cfp; currentSA];
        time_cfp = [time_cfp; currentTime];
    end
    
end
clear currentID currentSA currentTime


% 9. bin area by time
binned_yfp = accumarray(time_yfp,area_yfp,[],@(x) {x});
binned_cfp = accumarray(time_cfp,area_cfp,[],@(x) {x});


% 10. calculate mean, standard dev, counts, and standard error
y_bin_means = cellfun(@mean,binned_yfp);
y_bin_sums = cellfun(@sum,binned_yfp);
y_bin_stds = cellfun(@std,binned_yfp);
y_bin_counts = cellfun(@length,binned_yfp);
y_bin_sems = y_bin_stds./sqrt(y_bin_counts);

c_bin_means = cellfun(@mean,binned_cfp);
c_bin_sums = cellfun(@sum,binned_cfp);
c_bin_stds = cellfun(@std,binned_cfp);
c_bin_counts = cellfun(@length,binned_cfp);
c_bin_sems = c_bin_stds./sqrt(c_bin_counts);


% 11. normalize sums by counts
y_norm = y_bin_sums./y_bin_counts;
c_norm = c_bin_sums./c_bin_counts;


% 12. normalize by initial sum
y_norm2 = y_norm/y_norm(1);
c_norm2 = c_norm/c_norm(1);


% 13. plot growth rate over time
palette = {'DodgerBlue','GoldenRod'};

yfp_color = rgb(palette(2));
cfp_color = rgb(palette(1));
xmark = 'o';

figure(1)
plot((1:length(y_norm2))/binsPerHour,y_norm2,'Color',yfp_color)
hold on
plot((1:length(c_norm2))/binsPerHour,c_norm2,'Color',cfp_color)
hold on
grid on
axis([0,8.5,.9,1.2])
xlabel('Time (hr)')
ylabel('Area per cell / initial')
title(strcat(date,': growth in total area'))
legend('YFP WT', 'CFP mutant')


%% section 3. dA/dt, where A is normalized by the number of cells


% Strategy:
%                0. initialize experiment parameters and data
%                1. assemble data matrix
%                2. truncate data to non-erroneous timestamps (e.g. avoid bubbles) 
%                3. bin growth rate into time bins based on timestamp
%                4. isolate YFP and CFP intensities and area
%                5. convert intensities to (+) or (-) fluorophore
%                       - for details on threshold determination,
%                         see scrap_pile section "YFP and CFP thresholds"
%                6. test threshold, throw error if any data points ID as both
%                7. separate area by fluorophore
%                8. bin area by time
%                9. calculate stats for each bin, including summed area
%               10. normalize sums by cell counts
%               11. calculate change in area per cell per timestep
%               12. plot growth rate over time


clear
clc

% 0. initialize data
%xy = 2;
xy_start = 1;
xy_end = 16;
dt_min = 2;
%dt_min = 30; % reduced frequency dataset

date = '2018-11-23';
%cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date,'_xy02'))
%load(strcat('glycogen-',date,'-earlyEdits-jiggle-0p5.mat'),'D5');
cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date))
load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');


% 0. define time binning parameters
specificBinning = 2; % in minutes
binsPerHour = 60/specificBinning;


% 0. define fluorescence intensity threshold
threshold = 103.4;


% 1. assemble data matrix
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
clear xy_start xy_end


% 2. truncate data to non-erroneous timestamps (e.g. bubbles) 
maxTime = 8; % in hours
frame = glycogen_data(:,9);      % col 9 = frame in image sequence
dt_sec = dt_min * 60;
timeInSeconds = frame * dt_sec;  % frame = is consequetive images in analysis
timeInHours = timeInSeconds/3600;

if maxTime > 0
    glycogenData_bubbleTrimmed = glycogen_data(timeInHours <= maxTime,:);
else
    glycogenData_bubbleTrimmed = glycogen_data;
end
clear maxTime frame timeInSeconds timeInHours



% 3. bin growth rate into time bins based on timestamp
frame = glycogenData_bubbleTrimmed(:,9);      % col 9 = frame in image sequence
timeInSeconds = frame * dt_sec;
timeInHours = timeInSeconds/3600;
bins = ceil(timeInHours*binsPerHour);
clear timeInSeconds frame



% 4. isolate YFP and CFP intensities and area
cfp = glycogenData_bubbleTrimmed(:,13);         % col 13 = mean CFP intensity
yfp = glycogenData_bubbleTrimmed(:,14);         % col 14 = mean YFP intensity
area = glycogenData_bubbleTrimmed(:,6);         % col 6 = measured surface area



% 5. convert intensities to (+) or (-) fluorophore
isCFP = cfp > threshold;
isYFP = yfp > threshold;



% 6. test threshold, throw error if any data points ID as both
isBoth = isCFP+isYFP;
if sum(isBoth == 2) > 0
    error('threshold fail! some cells positive for both fluorophores')
end

% this step doesn't matter as long as threshold doesn't allow double-positives
glycogenData_final = glycogenData_bubbleTrimmed(isBoth < 2,:);
bins_final = bins(isBoth < 2);
isYFP_final = isYFP(isBoth < 2);
isCFP_final = isCFP(isBoth < 2);
area_final = area(isBoth < 2);


% 7. separate area by fluorophore
trackNum = glycogenData_final(:,12); % col 12 = track num
IDs_yfp = unique(trackNum(isYFP_final == 1));
IDs_cfp = unique(trackNum(isCFP_final == 1));

IDs_labeled = [IDs_yfp; IDs_cfp];

area_yfp = [];
area_cfp = [];
time_yfp = [];
time_cfp = [];
for id = 1:length(IDs_labeled)
    
    currentID = IDs_labeled(id);
    currentSA = area_final(trackNum == currentID);
    currentTime = bins_final(trackNum == currentID);
    
    if ismember(currentID,IDs_yfp) == 1 % if ID belongs to YFP+
        area_yfp = [area_yfp; currentSA];
        time_yfp = [time_yfp; currentTime];
    else % else
        area_cfp = [area_cfp; currentSA];
        time_cfp = [time_cfp; currentTime];
    end
    
end
clear currentID currentSA currentTime


% 8. bin area by time
binned_yfp = accumarray(time_yfp,area_yfp,[],@(x) {x});
binned_cfp = accumarray(time_cfp,area_cfp,[],@(x) {x});


% 9. calculate mean, standard dev, counts, and standard error
y_bin_means = cellfun(@mean,binned_yfp);
y_bin_sums = cellfun(@sum,binned_yfp);
y_bin_stds = cellfun(@std,binned_yfp);
y_bin_counts = cellfun(@length,binned_yfp);
y_bin_sems = y_bin_stds./sqrt(y_bin_counts);

c_bin_means = cellfun(@mean,binned_cfp);
c_bin_sums = cellfun(@sum,binned_cfp);
c_bin_stds = cellfun(@std,binned_cfp);
c_bin_counts = cellfun(@length,binned_cfp);
c_bin_sems = c_bin_stds./sqrt(c_bin_counts);


% 10. normalize sums by counts
y_norm = y_bin_sums./y_bin_counts;
c_norm = c_bin_sums./c_bin_counts;


% 11. calculate change in area per cell per timestep
y_dAdt = diff(y_norm)./dt_min;
c_dAdt = diff(c_norm)./dt_min;


% 12. plot growth rate over time
palette = {'DodgerBlue','GoldenRod'};

yfp_color = rgb(palette(2));
cfp_color = rgb(palette(1));
xmark = 'o';

figure(1)
plot((1:length(y_dAdt))/binsPerHour,y_dAdt,'Color',yfp_color)
hold on
plot((1:length(c_dAdt))/binsPerHour,c_dAdt,'Color',cfp_color)
hold on
grid on
%axis([0,8.5,.9,1.2])
xlabel('Time (hr)')
ylabel('dA/dt')
title(strcat(date,': change in area per cell over time'))
legend('YFP WT', 'CFP mutant')
