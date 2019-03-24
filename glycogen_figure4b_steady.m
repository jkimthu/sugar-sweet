%% glycogen figure4b - STEADY

% Goal: produce final manuscript figures from microfluidic visualizations
%
%       input:  both steady replicates with YFP wt and CFP mutant

%       output: two plots (individual and replicate stats) of growth rate vs time
%               mean is the average growth rate between experiments
%               errorbars are standard deviation between replicate means


% Note: major difference from "pulsing" version of this script is the time
%       cut-off. steady data also goes out to 8h but after 5h the cells are
%       so elongated that any tracked particle after this point is error.


% last update: Jen, 2019 Mar 24
% commit: first commit, final figure for manuscript

% ok let's go!

%% %% A. initialize

clc
clear

% 0. initialize experiments to include in analysis
dates = {'2019-02-19','2019-02-25'};
counter = 0;


% 0. initialize meta data
xy_start = 25;
xy_ends = [29; 30];
dt_min = 3;


% 0. define growth rates of interest
specificGrowthRate = 'log2';


% 0. define time binning parameters
specificBinning = 60; % in minutes
binsPerHour = 60/specificBinning;


% 0. define fluorescence intensity threshold
threshold = 103.4;


% 0. initialize color designations
color_yfp = rgb('GoldenRod'); 
color_cfp = rgb('DodgerBlue'); 


%% Part 1. curves for individual replicates


% 1. loop through each experiment of interest
for rep = 1:length(dates)
    
    
    % 2. initialize experiment meta data
    date = dates{rep};
    disp(strcat(date, ': analyze!'))
    
    
    % 3. load measured data
    filename = strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat');
    load(filename,'D5')
    
    
    % 4. build data matrix
    xy_end = xy_ends(rep);
    repData = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
    clear xy_end
    
    
    % 5. isolate volume (Va), drop, and track number for growth rate calculations
    volumes = repData(:,5);        % col 5 = calculated va_vals (cubic um)
    isDrop = repData(:,3);         % col 3 = isDrop, 1 marks a birth event
    trackNum = repData(:,12);      % col 12 = track number (not ID from particle tracking)
    
    
    % 6. calculate growth rate
    dt_sec = dt_min * 60;
    growthRates = calculateGrowthRate_glycogen(volumes,isDrop,trackNum,dt_sec);
    clear isDrop volumes trackNum
    
    
    % 7. truncate data to non-erroneous timestamps (e.g. bubbles)
    maxTime = 8; % in hours
    frame = repData(:,9);            % col 9 = frame in image sequence
    timeInSeconds = frame * dt_sec;  % frame = is consequetive images in analysis
    timeInHours = timeInSeconds/3600;
    
    if maxTime > 0
        glycogenData_maxTrimmed = repData(timeInHours <= maxTime,:);
        growthRates_maxTrimmed = growthRates(timeInHours <= maxTime,:);
    else
        glycogenData_maxTrimmed = repData;
        growthRates_maxTrimmed = growthRates;
    end
    clear maxTime frame timeInSeconds timeInHours growthRates

    
    
    % 8. bin time based on timestamp
    frame = glycogenData_maxTrimmed(:,9);      % col 9 = frame in image sequence
    timeInSeconds = frame * dt_sec;
    timeInHours = timeInSeconds/3600;
    bins = ceil(timeInHours*binsPerHour);
    clear timeInSeconds frame

    
    
    % 9. isolate selected specific growth rate and remove nans from data analysis
    specificColumn = 3;
    growthRate_log2 = growthRates_maxTrimmed(:,specificColumn);
    
    growthRt_noNaNs = growthRate_log2(~isnan(growthRate_log2),:);
    bins_noNaNs = bins(~isnan(growthRate_log2),:);
    glycogenData_noNaNs = glycogenData_maxTrimmed(~isnan(growthRate_log2),:);
    clear glycogenData_maxTrimmed growthRates_maxTrimmed growthRate_log2 bins
    
     
    
    % 10. isolate YFP and CFP intensities
    cfp = glycogenData_noNaNs(:,13);         % col 13 = mean CFP intensity
    yfp = glycogenData_noNaNs(:,14);         % col 14 = mean YFP intensity
    
    
    
    % 11. convert intensities to (+) or (-) fluorophore
    isCFP = cfp > threshold;
    isYFP = yfp > threshold;
    clear cfp yfp
    
    
    
    % 12. test threshold, throw out any tracks that at any point ID as both
    isBoth = isCFP + isYFP;
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
    clear isCFP isYFP currentTrack i glycogenData_noNaNs bins_noNaNs isBoth
    clear unique_errors_sorted unique_errors ids_errors idx_errors
    
    
    
    % 13. prep data for fluorophore sorting
    growthRt_final = growthRt_noDoublePos;
    glycogenData_final = glycogenData_noDoublePos;
    bins_final = bins_noDoublePos;
    isCFP_final = glycogenData_final(:,13) > threshold;   % col 13 = mean CFP intensity
    isYFP_final = glycogenData_final(:,14) > threshold;   % col 14 = mean YFP intensity;
    clear glycogenData_noDoublePos growthRt_noDoublePos bins_noDoublePos
    
    
    
    % 14. separate growth rates by fluorophore
    trackNum = glycogenData_final(:,12);            % col 12 = track num
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
    clear currentID currentGR currentBins trackNum
    clear isCFP_final isYFP_final IDs_labeled IDs_yfp IDs_cfp id

    
    
    % 15. bin growth rate by time and calculate mean
    binned_yfps(:,rep) = accumarray(bins_yfp,growthRt_yfp,[],@(x) {x});
    binned_cfps(:,rep) = accumarray(bins_cfp,growthRt_cfp,[],@(x) {x});
    
    binned_yfp_mean{rep} = accumarray(bins_yfp,growthRt_yfp,[],@mean);
    binned_cfp_mean{rep} = accumarray(bins_cfp,growthRt_cfp,[],@mean);
    clear growthRt_cfp growthRt_yfp bins_yfp bins_cfp
    
  
    
    % 16. plot response in growth rate for all timescales over time
    cutoff = 5; % hr
    timeVector = (1:length(binned_yfp_mean{rep}))/binsPerHour;
    times = timeVector(timeVector <= cutoff);

    figure(1)   
    plot(times,binned_yfp_mean{rep}(1:length(times)),'Color',color_yfp,'LineWidth',1,'Marker','.')
    hold on
    plot(times,binned_cfp_mean{rep}(1:length(times)),'Color',color_cfp,'LineWidth',1,'Marker','.')
    title(strcat('growth in steady nutrient: YFP wt vs CFP delta-glgP'))
    xlabel('time (hr)')
    ylabel(strcat('growth rate: (', specificGrowthRate,')'))
    axis([0.9,5.1,-0.03,2])
    legend('YFP WT', 'CFP mutant')

     
end
clear rep date D5 growthRt_noNaNs


%% Part 2. find mean and standard deviation between replicates


% 1. compile data from each replicate into a single matrix (row = replicate)
for rep = 1:length(dates)
    
    col = rep;
    
    rep_means_yfp = binned_yfp_mean{col}(1:length(times));
    compiled_yfps(rep,:) = rep_means_yfp;
    
    rep_means_cfp = binned_cfp_mean{col}(1:length(times));
    compiled_cfps(rep,:) = rep_means_cfp;

    
end
clear rep rep_means_cfp rep_means_yfp col


% 2. calculate mean and stdev across replicates (rows)
compiled_yfp_means = mean(compiled_yfps);
compiled_yfp_stds = std(compiled_yfps);

compiled_cfp_means = mean(compiled_cfps);
compiled_cfp_stds = std(compiled_cfps);



% 3. plot mean
figure(2)
plot(times,compiled_yfp_means,'Color',color_yfp,'LineWidth',1,'Marker','.')
hold on
plot(times,compiled_cfp_means,'Color',color_cfp,'LineWidth',1,'Marker','.')
title(strcat('growth in steady nutrient: YFP wt vs CFP delta-glgP'))
ylabel('log 2 growth rate (1/hr)')
xlabel('time (hr)')
axis([0.9,5.1,-0.03,2])


% 4. plot shaded errorbars (st dev)
figure(2)
hold on
ss = shadedErrorBar(times,compiled_yfps,{@nanmean,@nanstd},'lineprops',{'Color',color_yfp},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
ss.patch.FaceColor = color_yfp;
axis([0.9,5.1,-0.03,2])

hold on
ss = shadedErrorBar(times,compiled_cfps,{@nanmean,@nanstd},'lineprops',{'Color',color_cfp},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
ss.patch.FaceColor = color_cfp;
axis([0.9,5.1,-0.03,2])


%% Part 3. save data and figures

% save mean signal of each replicate over time
% in compiled data matrices,
%           each row is a replicate (refer to variable: dates)
%           column is a timepoint (refer to variable:times)

save('glycogen_mu_vs_time_steady.mat','compiled_yfps','compiled_cfps','times','dates')

figure(1)
plotName = strcat('glycogen-mu-v-time-steady-fig1-individualReplicates');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(2)
plotName = strcat('glycogen-mu-v-time-steady-fig2-replicateStats');
saveas(gcf,plotName,'epsc')
close(gcf)


