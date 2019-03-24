%% quantifyGrowth

% Goal: quantify and visualize growth rates of WT and mutants for glycogen
%       utilization under a pulsling nutrient environment.
%
%       input data is data structure produced from glycogen_analysis.m

%       this plots instantaneous growth rate in a few ways:
%
%           section 1. instantaneous rate of volumetric doubling over full experimental time



% last update: Jen, 2019 Mar 10
% commit: first analysis for 2019-03-09 experiment, steady

% ok let's go!

%% section 1. instantaneous rate of volumetric doubling


% Strategy:
%                0. initialize experiment parameters and data
%                1. assemble data matrix
%                2. isolate volume (Va), drop, and track number for growth rate calculations
%                3. calculate growth rate
%                4. truncate data to non-erroneous timestamps (e.g. bubbles) 
%                5. bin time based on timestamp
%                6. isolate selected specific growth rate and remove nans from data analysis
%                7. isolate YFP and CFP intensities
%                8. convert intensities to (+) or (-) fluorophore
%                       - for details on threshold determination,
%                         see scrap_pile section "YFP and CFP thresholds"
%                9. test threshold, throw error if any data points ID as both
%               10. prep data for fluorophore sorting 
%               11. separate growth rates by fluorophore
%               12. bin growth rate by time
%               13. calculate stats for each bin, including mean growth rate
%               14. plot growth rate over time



clear
clc

% 0. initialize data

xy_start = 11;
xy_end = 20;
dt_min = 3;

date = '2019-03-08';


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
volumes = glycogen_data(:,5);        % col 5 = calculated va_vals (cubic um)
isDrop = glycogen_data(:,3);         % col 3 = isDrop, 1 marks a birth event
trackNum = glycogen_data(:,12);      % col 12 = track number (not ID from particle tracking)


% 3. calculate growth rate
dt_sec = dt_min * 60;
growthRates = calculateGrowthRate_glycogen(volumes,isDrop,trackNum,dt_sec);
clear isDrop volumes trackNum dt_min


% 4. truncate data to non-erroneous timestamps (e.g. bubbles) 

maxTime = 5; % in hours

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


% 5. bin time based on timestamp
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


% 10. prep data for fluorophore sorting 
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
% out IDs that are identified with YFP or CFP and use these to isolate
% growth rates.

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

axis([0,8.5,-0.05,2])

xlabel('Time (hr)')
ylabel('Growth rate (1/hr)')
title(strcat(date,': (',specificGrowthRate,')'))
%legend('YFP WT', 'CFP mutant')
legend('YFP mutant', 'CFP wt')


