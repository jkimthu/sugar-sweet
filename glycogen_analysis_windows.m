% Glycogen analysis for Windows

% Goal: determine growth rates of WT and glycogen breakdown mutants under a
%       pusling MPG environment.
%
%       input data is a timelapse with three channels: phase, YFP and CFP
%       stratgey:
%                0. initialize experiment parameters
%                1. particle identification, characterization and tracking
%                2. quality control: clean dataset prior to growth rate calculations
%                3. fun part! compute population growth rates between YFP and CFP populations
%                4. 
%

% BIIIIG thank yous to Cherry Gao for helping me get started!

% last update: Jen, 2019 January 15
% commit: add eccentricity and orientation to 2018-11-23 analysis


% ok let's go!

%% 0. initialize experiment parameters

clear
clc

%  initialize conversion from pixels to um based on camera and magnification
%  scope5 Andor COSMOS = 6.5um pixels / 60x magnification
ConversionFactor = 6.5/60; % units = um/pixel


% initialize experiment folder
experiment = '2018-11-23';


% % open folder for experiment of interest
% imageFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',experiment);
% cd(imageFolder); % move to the folder containing images 


% initialize image names
imageName = strcat('glycogen-combo1-',experiment);


% initialize total number of movies and timepoints
xy_final = 16; % total num of xy positions in analysis
num_tpt = 241; % total num of timepoints in analysis


%% 1. particle identification, characterization and tracking
%
%  Strategy:
%
%         a. identification of particles (cells) in phase
%
%         b. determine particle properties, including size and
%              fluorescent label 
 %      
%                    1) compute the YFP fluorescence intensity of each cell
%                    2) choose a threshold for determining if a cell is on or off in YFP
%                    3) output the cell ID
%                    4) do the same for CFP, and make sure that the cell ID numbers for
%                       YFP-positive cells & the ID numbers for CFP-positive cells are
%                       mutually exclusive
%
%         c. track particles through time based on xy coordinate

%% ONE. TEST ANALYSIS
%  visualize image analysis parameters on single image
xy = 16;
tpt = 50;

% initialize xy name
if xy < 10
    xyName = strcat('xy0',num2str(xy));
elseif xy < 100
    xyName = strcat('xy',num2str(xy));
end

% initialize timepoint name
if tpt < 10
    tptName = strcat('t00',num2str(tpt));
elseif tpt < 100
    tptName = strcat('t0',num2str(tpt));
elseif tpt < 1000
    tptName = strcat('t',num2str(tpt));
end


% define image name of each channel from current tpt, in current xy
ph_name = strcat(imageName,tptName,xyName,'c1.tif');
c_name = strcat(imageName,tptName,xyName,'c2.tif');
y_name = strcat(imageName,tptName,xyName,'c3.tif');


% load image into matlab with function "imread"
ph_image = imread(ph_name);
c_image = imread(c_name);
y_image = imread(y_name);
clear ph_name c_name y_name



%figure;
subplot(131); imshow(imadjust(ph_image)); title('phase');         % "imadjust" automatically adjusts brightness of 16-bit image
subplot(132); imshow(imadjust(c_image)); title('cfp');
subplot(133); imshow(imadjust(y_image)); title('yfp');

%%
% edge detection of cells in phase image
bw = edge(ph_image,'sobel');


% clean up edge detection
se = strel('disk',3); % structuring element; disk-shaped with radius of 3 pixels
bw_dil = imdilate(bw, se); % dilate
bw_fill = imfill(bw_dil, 'holes'); % fill
bw_final = imerode(imerode(bw_fill,se),se); % erode twice to smooth; this is your final mask


% visualize effect of each step of edge detection and cleanup
figure;
ax(1) = subplot(231); imshowpair(imadjust(ph_image),bw_final); title('final mask-phase overlay');
ax(2) = subplot(232); imshow(bw); title('step 1: edge detection');
ax(3) = subplot(233); imshow(bw_dil); title('step 2: dilate edges');
ax(4) = subplot(234); imshow(bw_fill); title('step 3: fill');
ax(5) = subplot(235); imshow(bw_final); title('step 4: erode or smooth');
linkaxes(ax(:),'xy'); % link subplots so that you can zoom-in on all subplots simultaneously

%%
% segment cells
cc = bwconncomp(bw_final);


% gather properties for each identified particle
stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');


% visualize distribution of cell areas
area = extractfield(stats,'Area')';
area_conv = area.*ConversionFactor^2;
figure; histogram(area_conv)
xlabel('Area (sq microns)')
ylabel('Count')



%%
% 2. determine which cells are YFP or CFP labeled
% overlay mask with YFP fluorescence image by dot-product
y_mask_overlay = bw_final .* double(y_image); % convert your uint16 image to double class
c_mask_overlay = bw_final .* double(c_image);

% compute the fluorescence intensity value for each cell
for cell = 1:cc.NumObjects
    
    pixel_id_of_cell = []; % initialize
    pixel_id_of_cell = cc.PixelIdxList{cell}; % pixel index where the cell of interest is
    
    y_fluo_int_cells(cell) = mean(y_mask_overlay(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
    c_fluo_int_cells(cell) = mean(c_mask_overlay(pixel_id_of_cell));
    
end


% look at the distribution of YFP intensities of cells to determine a threshold
figure; histogram(y_fluo_int_cells,50); title('yfp'); xlabel('mean intensity')
figure; histogram(c_fluo_int_cells,50); title('cfp'); xlabel('mean intensity')

%%
% from the histogram above, it seems like a threshold of around 100 might work
threshold_yfp = 101;
threshold_cfp = 103;

% let's see how good the threshold is!
% in this for-loop, you will
% 1) erase the masks of YFP-negative cells
% 2) output the cell ID numbers of YFP-positive cells

bw_final_y = bw_final; % initialize
bw_final_c = bw_final; % initialize

for cell = 1:cc.NumObjects
    
    % yfp
    if y_fluo_int_cells(cell) < threshold_yfp
        bw_final_y(cc.PixelIdxList{cell}) = 0;   % erase the YFP-negative cells
        stats(cell).yfp = 0;
    else
        stats(cell).yfp = 1;
    end
    
    % cfp
    if c_fluo_int_cells(cell) < threshold_cfp
        bw_final_c(cc.PixelIdxList{cell}) = 0;    % erase the CFP-negative cells
        stats(cell).cfp = 0;
    else
        stats(cell).cfp = 1;
    end
    
end

% let's see if the threshold worked
figure; imshowpair(bw_final_y,imadjust(y_image)); title('after thresholding, remaining YFP masks')
figure; imshowpair(bw_final_c,imadjust(c_image)); title('after thresholding, remaining CFP masks')

%%
% add YFP or CFP designation to data structure
%           0 = not above threshold,
%           1 = YES label! (above threshold)


% 3. assemble data structure similar to P from ND2Proc_XY.m in order to use track linker

% x & y coord
X = [];
Y = [];

for particle = 1:length(stats)
    
    centroid = stats(particle).Centroid.*ConversionFactor;
    X = [X; centroid(1)];
    Y = [Y; centroid(2)];
    
end

p_unit.X = X;
p_unit.Y = Y;
clear particle X Y centroid


% area
area = extractfield(stats,'Area')';
p_unit.A = area.*ConversionFactor^2;
clear area


% major & minor axis
majorAxis = extractfield(stats,'MajorAxisLength')';
minorAxis = extractfield(stats,'MinorAxisLength')';
p_unit.MajAx = majorAxis.*ConversionFactor;
p_unit.MinAx = minorAxis.*ConversionFactor;
clear majorAxis minorAxis


% YFP or CFP cell?
isYFP = extractfield(stats,'yfp')';
isCFP = extractfield(stats,'cfp')';
p_unit.yfp = isYFP;
p_unit.cfp = isCFP;
clear isYFP isCFP

% frame
p_unit.Frame = tpt;


% eccentricity and angle
ecc = extractfield(stats,'Eccentricity')';
angle = extractfield(stats,'Orientation')';
p_unit.Ecc = ecc;
p_unit.Angle = angle;


p_clone(tpt) = p_unit;

%% TWO. FULL ANALYSIS
%  for each xy, loop through all timepoints and build data structure

for xy = 1:xy_final
    
    % initialize xy name
    if xy < 10
        xyName = strcat('xy0',num2str(xy));
    elseif xy < 100
        xyName = strcat('xy',num2str(xy));
    end
    disp( strcat('Analyzing xy(',num2str(xy),') of (',num2str(xy_final),')') )
    
    for tpt = 1:num_tpt % num tpts
        
        % initialize timepoint name
        if tpt < 10
            tptName = strcat('t00',num2str(tpt));
        elseif tpt < 100
            tptName = strcat('t0',num2str(tpt));
        elseif tpt < 1000
            tptName = strcat('t',num2str(tpt));
        end
        
        
        % define image name of each channel from current tpt, in current xy
        ph_name = strcat(imageName,tptName,xyName,'c1.tif');
        c_name = strcat(imageName,tptName,xyName,'c2.tif');
        y_name = strcat(imageName,tptName,xyName,'c3.tif');
        
        
        % load image into matlab with function "imread"
        ph_image = imread(ph_name);
        c_image = imread(c_name);
        y_image = imread(y_name);
        clear ph_name c_name y_name
        
        

        %figure;
        %subplot(131); imshow(imadjust(ph_image)); title('phase');         % "imadjust" automatically adjusts brightness of 16-bit image
        %subplot(132); imshow(imadjust(c_image)); title('cfp');
        %subplot(133); imshow(imadjust(y_image)); title('yfp');
        
        
        % edge detection of cells in phase image
        bw = edge(ph_image,'sobel');
        
        
        % clean up edge detection
        se = strel('disk',3); % structuring element; disk-shaped with radius of 3 pixels
        bw_dil = imdilate(bw, se); % dilate
        bw_fill = imfill(bw_dil, 'holes'); % fill
        bw_final = imerode(imerode(bw_fill,se),se); % erode twice to smooth; this is your final mask
        
        
        % visualize effect of each step of edge detection and cleanup
        %figure;
        %ax(1) = subplot(231); imshowpair(imadjust(ph_image),bw_final); title('final mask-phase overlay');
        %ax(2) = subplot(232); imshow(bw); title('step 1: edge detection');
        %ax(3) = subplot(233); imshow(bw_dil); title('step 2: dilate edges');
        %ax(4) = subplot(234); imshow(bw_fill); title('step 3: fill');
        %ax(5) = subplot(235); imshow(bw_final); title('step 4: erode or smooth');
        %linkaxes(ax(:),'xy'); % link subplots so that you can zoom-in on all subplots simultaneously
        
        
        % segment cells
        cc = bwconncomp(bw_final);
        
        
        % gather properties for each identified particle
       stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');

        
        % visualize distribution of cell areas
        %area = extractfield(stats,'Area')';
        %area_conv = area.*ConversionFactor^2;
        %figure; histogram(area_conv)
        
        clear bw bw_dil bw_fill
        
        
        
        % 2. determine which cells are YFP or CFP labeled        
        % overlay mask with YFP fluorescence image by dot-product
        y_mask_overlay = bw_final .* double(y_image); % convert your uint16 image to double class
        c_mask_overlay = bw_final .* double(c_image);
        
        % compute the fluorescence intensity value for each cell
        for cell = 1:cc.NumObjects
            
            pixel_id_of_cell = []; % initialize
            pixel_id_of_cell = cc.PixelIdxList{cell}; % pixel index where the cell of interest is
            
            y_fluo_int_cells(cell) = mean(y_mask_overlay(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
            c_fluo_int_cells(cell) = mean(c_mask_overlay(pixel_id_of_cell));
            
        end
        
        
        % look at the distribution of YFP intensities of cells to determine a threshold
        %figure; histogram(y_fluo_int_cells,50); title('yfp');
        %figure; histogram(c_fluo_int_cells,50); title('cfp');
        
        % from the histogram above, it seems like a threshold of around 100 might work
        threshold_yfp = 101; % for first frame of 2018-11-23
        threshold_cfp = 103; % these thresholds give <10% cells ID as both, ~30% as neither
        
        
        % let's see how good the threshold is!
        % in this for-loop, you will
        % 1) erase the masks of YFP-negative cells
        % 2) output the cell ID numbers of YFP-positive cells
        
        bw_final_y = bw_final; % initialize
        bw_final_c = bw_final; % initialize

        for cell = 1:cc.NumObjects
            
            % yfp
            if y_fluo_int_cells(cell) < threshold_yfp
                %bw_final_y(cc.PixelIdxList{cell}) = 0;   % erase the YFP-negative cells
                stats(cell).yfp = 0;
            else
                stats(cell).yfp = 1;
            end
            
            % cfp
            if c_fluo_int_cells(cell) < threshold_cfp
                %bw_final_c(cc.PixelIdxList{cell}) = 0;    % erase the CFP-negative cells
                stats(cell).cfp = 0;
            else
                stats(cell).cfp = 1;
            end
            
        end
        
        % let's see if the threshold worked
        %figure; imshowpair(bw_final_y,imadjust(y_image)); title('after thresholding, remaining masks')
        %figure; imshowpair(bw_final_c,imadjust(c_image)); title('after thresholding, remaining masks')
        
        % add YFP or CFP designation to data structure
        %           0 = not above threshold,
        %           1 = YES label! (above threshold)
        clear ph_image y_image c_image y_mask_overlay c_mask_overlay
        
        
        
        
        % 3. assemble data structure similar to P from ND2Proc_XY.m in order to use track linker
        
        % x & y coord
        X = [];
        Y = [];
        
        for particle = 1:length(stats)
            
            centroid = stats(particle).Centroid.*ConversionFactor;
            X = [X; centroid(1)];
            Y = [Y; centroid(2)];
            
        end
        
        p_unit.X = X;
        p_unit.Y = Y;
        clear particle X Y centroid
        
        
        % area
        area = extractfield(stats,'Area')';
        p_unit.A = area.*ConversionFactor^2;
        clear area
        
        
        % major & minor axis
        majorAxis = extractfield(stats,'MajorAxisLength')';
        minorAxis = extractfield(stats,'MinorAxisLength')';
        p_unit.MajAx = majorAxis.*ConversionFactor;
        p_unit.MinAx = minorAxis.*ConversionFactor;
        clear majorAxis minorAxis
        
        
        % YFP or CFP cell?
        isYFP = extractfield(stats,'yfp')';
        isCFP = extractfield(stats,'cfp')';
        p_unit.yfp = isYFP;
        p_unit.cfp = isCFP;
        clear isYFP isCFP
        
        
        % frame
        p_unit.Frame = tpt;
        
        
        % eccentricity and angle
        ecc = extractfield(stats,'Eccentricity')';
        angle = extractfield(stats,'Orientation')';
        p_unit.Ecc = ecc;
        p_unit.Angle = angle;
        
        
        p_clone(tpt) = p_unit;

    end
    
    %
    
    % trim particles by (1) area and (2) width
    TrimField = 'A';    % choose relevant characteristic to restrict, run several times to apply for several fields
    LowerBound = 0.8;   % lower bound for restricted field, or -Inf
    UpperBound = 8;     % upper bound for glucose only
    
    % to actually trim the set:
    p_clone_trim1 = ParticleTrim_glycogen(p_clone,TrimField,LowerBound,UpperBound);
    
    
    %
    % trim particles by width
    
    TrimField = 'MinAx';  % choose relevant characteristic to restrict, run several times to apply for several fields
    LowerBound = 1.0;     % lower bound for restricted field, or -Inf
    UpperBound = 1.4;     % upper bound, reset 2018-11-25

    % to actually trim the set:
    p_clone_trim2 = ParticleTrim_glycogen(p_clone_trim1,TrimField,LowerBound,UpperBound);
    
    
    %
    % v. track remaining particles based on coordinate distance
    
    TrackMode = 'position';       % Choice of {position, velocity, acceleration} to predict position based on previous behavior
    DistanceLimit = 5;            % Limit of distance a particle can travel between frames, in units defined by ConversionFactor
    MatchMethod = 'best';         % Choice of {best, single}
    p_Tracks = Particle_Track_glycogen(p_clone_trim2,TrackMode,DistanceLimit,MatchMethod);
    
    
    % 4. Store tracked data into single workspace, D
    D{xy} = p_Tracks;
    
    disp(strcat('analyzed xy(',num2str(xy),') of (',num2str(xy_final),')!'))
    
    
end

save(strcat('glycogen-',experiment,'-width1p4.mat'),'D')  


%% 2. quality control: clean dataset prior to growth rate calculations

%  Selection criteria: rules particles must obey to remain in analysis

%   1. Tracks should not increase by more than 30% of size at previous timepoint
%           - huge jumps in cell size are cells coming together
%           - tracks are clipped, not deleted
%           - data prior and after jump are considered separate tracks 

%   2. Tracks must be at least some number of frames long
%           - we are looking for rate of change, so some continuity is important

%   3. Tracks cannot oscillate too quickly between gains and losses
%           - jiggly tracks correspond to non-growing particles with noise

%   4. Tracks must be of reasonable size, at least SizeStrainer (1.8um)



% particle tracking data
clear
clc
experiment = '2018-11-23';

% 0. open folder for experiment of interest
load(strcat('glycogen-',experiment,'-width1p4.mat'));

% reject data matrix
rejectD = cell(4,length(D));

% criteria counter
criteria_counter = 0;

% windowsize for mu calculations
windowSize = 5;



%  CRITERIA ONE
%  clip tracks to separate data surrounding >30% jumps in cell length
criteria_counter = criteria_counter + 1;

% 0. initialize
jumpFrac = 0.3;                                                            
                                                                           
for n = 1:length(D)                                                       
    
    % 0. isolate data for current movie
    data = D{n};
    
    % 0. if no data in n, continue to next movie
    if isempty(data) == 1
        disp(strcat('Clipping 0 jumps in D2 (', num2str(n),') !'))
        continue
    end
    
    % 0. for each track, search for large positive increases: > %30 jumps in size
    jump_counter = 0;
    for m = 1:length(data)
        
        %0. determine if track is shorter than window size. skip over, if so. it will be trimmed later.
        if length(data(m).X) < windowSize
            continue
        end
        
        % 1. determine change in length between each timestep                           
        Rates = diff(data(m).MajAx);
        
        % 2. express change of as a fraction of the cell size in previous timestep
        Lengths = data(m).MajAx(1:length(Rates));
        growthFrac = Rates./Lengths;
                                                                           
        % 3. list rows in which the change exceeds size jump threshold
        jumpPoints = find(growthFrac > jumpFrac);                        

        
        % 4. if the track contains jumps...
        if isempty(jumpPoints) == 0
            
            % i. isolate structure of target track z, in prep to clip all variables (MajAx, X, Y, etc.)
            originalTrack = data(m);
            
            % ii. define timepoint of earliest jump
            clipPoint = jumpPoints(1);
            
            % iii. clip data structure at desired timepoint
            clippedTarget = structfun(@(M) M(1:clipPoint), originalTrack, 'Uniform', 0);
            
            % iv. redefine track in data set as clipped structure
            data(m) = clippedTarget;
            
            % v. store remainder of original Target for reject data matrix
            remainderTrack = structfun(@(M) M(clipPoint+1:end), originalTrack, 'Uniform', 0);
            
            jump_counter = jump_counter + 1;
            trackScraps(jump_counter,1) = remainderTrack;

            % 5. repeat for all jumpy tracks
        end
        
        clear remainderTrack clippedTarget clipPoint Lengths Rates growthFrac originalTrack jumpPoints
        
    end
    
    % 6. report!
    X = ['Clipping ', num2str(jump_counter), ' jumps in D2(', num2str(n), ')...'];
    disp(X)
    
    % when all tracks finished, 
    if jump_counter >= 1
        
        % 7. save accmulated rejects
        rejectD{criteria_counter,n} = trackScraps;
        
        % 8. and insert data after jump at end of data matrix, D2
        D2{n} = [data; trackScraps];
    
    else
        D2{n} = data;
    end
    
    % 9. repeat for all movies
    clear trackScraps X data;
    
end
    
clear  jumpFrac jump_counter Rates m n;



%  CRITERIA TWO
%  tracks must be at least window size in length (5 frames)
criteria_counter = criteria_counter + 1;

% 0. initialize new dataset before trimming
D3 = D2;

for n = 1:length(D)
    
    % 0. if no data in n, continue to next movie
    if isempty(D3{n}) == 1

        disp(strcat('Removing (0) short tracks from D3 (', num2str(n),') !'))

        continue
    end
    
    % 1. determine number of timepoints in each track m 
    for m = 1:length(D3{n})
        numFrames(m) = length(D3{n}(m).MajAx);
    end
    
    % 2. find tracks that are shorter than threshold number of frames
    subThreshold = find(numFrames < windowSize);       
    
    % 3. report!
    X = ['Removing ', num2str(length(subThreshold)), ' short tracks from D3(', num2str(n), ')...'];
    disp(X)
    
    % 4. to that loop doesn't crash if nothing is too short
    if isempty(subThreshold) == 1
        continue
    end
    
    % 5. remove structures based on row # (in reverse order)
    fingers = 0;
    for toRemove = 1:length(subThreshold)
        
        r = length(subThreshold) - fingers;                  % reverse order                  
        tracks_shortGlimpses(r,1) = D3{n}(subThreshold(r));   % store data for reject data matrix
        D3{n}(subThreshold(r)) = [];                         % deletes data
        fingers = fingers + 1;
        
    end
    
    % 6. save sub-threshold tracks into reject data matrix
    rejectD{criteria_counter,n} = tracks_shortGlimpses;
    
    % 7. repeat for all movies
    clear  numFrames m subThreshold fingers toRemove r X tracks_shortGlimpses;
    
end

 clear windowSize n; 

 
 
 
%CRITERIA THREE: tracks cannot oscillate too quickly between gains and losses
criteria_counter = criteria_counter + 1;

% 0. initiaize new dataset before trimming
D4 = D3;

% 0. define threshold change in length considered a division
dropThreshold = -0.75;

% 0. define threshold under which tracks are too jiggly
jiggleThreshold = -0.5;

for n = 1:length(D)
    
    % 0. if no data in n, continue to next movie
    if isempty(D4{n}) == 1
        disp(strcat('Removing (0) jiggly tracks from D4 (', num2str(n),') !'))
        continue
    end
    
    % 1. for each track, collect % of non-drop negatives per track length (in frames)
    nonDropRatio = NaN(length(D{n}),1);
    
    for m = 1:length(D4{n})
        
        % 2. isolate length data from current track
        lengthTrack = D4{n}(m).MajAx;
        
        % 3. find change in length between frames
        diffTrack = diff(lengthTrack);
        
        % 4. convert change into binary, where positives = 0 and negatives = 1
        binaryTrack = logical(diffTrack < 0);
        
        % 5. find all drops (negatives that exceed drop threshold)
        dropTrack = diffTrack < dropThreshold;
        
        % 6. find the ratio of non-drop negatives per track length
        nonDropNegs = sum(dropTrack - binaryTrack);
        squiggleFactor = nonDropNegs/length(lengthTrack);
        
        % 7. store ratio for subsequent removal
        nonDropRatio(m) = squiggleFactor;
    
    % 8. repeat for all tracks
    end
    
    % 9. determine which tracks fall under jiggle threshold
    belowThreshold = find(nonDropRatio <= jiggleThreshold);
    
    if ~isempty(belowThreshold)
        
        % 10. report!
        X = ['Removing ', num2str(length(belowThreshold)), ' jiggly tracks from D4(', num2str(n), ')...'];
        disp(X)
        
        % 11. remove jiggly structures based on row # (in reverse order)
        counter = 0;
        for toRemove = 1:length(belowThreshold)
            
            r = length(belowThreshold) - counter;            % reverse order
            jigglers(r,1) = D4{n}(belowThreshold(r));        % store data for reject data matrix
            D4{n}(belowThreshold(r)) = [];                   % deletes data
            counter = counter + 1;
            
        end
        
        % 13. save sub-threshold tracks into reject data matrix
        rejectD{criteria_counter,n} = jigglers;
    else
        
        % 14. repeat for all movies
        clear nonDropRatio lengthTrack diffTrack dropTrack nonDropNegs squiggleFactor belowThreshold
        clear toRemove binaryTrack r X jigglers
        
        continue
        
    end
    
    
    % 14. repeat for all movies
    clear nonDropRatio lengthTrack diffTrack dropTrack nonDropNegs squiggleFactor belowThreshold
    clear toRemove binaryTrack r X jigglers
  
end

clear n gainLossRatio jiggleThreshold jigglers counter dropThreshold;



%  CRITERIA FOUR
%  maximum particle size must be greater than 1.8um
criteria_counter = criteria_counter + 1;

D5 = D4;
SizeStrainer = 1.8;

for n = 1:length(D)   
    
    % 1. so that loop doesn't crash if no data remaining in n
    if isempty(D5{n}) == 1
        X = ['Removing 0 small particles from D5(', num2str(n), ')...'];
        disp(X)
        continue
    end
    
    for i = 1:length(D5{n})
        lengthTrack(i) = max(D5{n}(i).MajAx);                          
    end          
    
    % finds tracks that don't exceed __ um
    tooSmalls = find(lengthTrack < SizeStrainer);                          
    
    % report!
    X = ['Removing ', num2str(length(tooSmalls)), ' small particles from D5(', num2str(n), ')...'];
    disp(X)
    
    % so loop doesn't crash if nothing is too small
    if isempty(tooSmalls) == 1
        continue
    end
    
    % remove too-small structures based on row # (in reverse order)
    countSmalls = 0;
    for s = 1:length(tooSmalls)
        t = length(tooSmalls) - countSmalls;
        tracks_tooSmalls(t,1) = D5{n}(tooSmalls(t));      %  recording to add into reject data matrix
        D5{n}(tooSmalls(t)) = [];
        countSmalls = countSmalls + 1;
    end
    
    % save tracks that are too small into reject data matrix
    rejectD{criteria_counter,n} = tracks_tooSmalls;
    clear lengthTrack lengthTrack_dbl i tooSmalls countSmalls s t X tracks_tooSmalls;

    
end 

clear SizeStrainer n i m tooSmalls X;

save(strcat('glycogen-',experiment,'-width1p4-jiggle-0p5.mat'), 'D', 'D2', 'D3', 'D4', 'D5', 'rejectD')
disp('Quality control: complete!')


%% 3. compute population growth rates between YFP and CFP populations

clear
clc

% define experiment of interest
date = '2018-11-23';

% define growth rates of interest
specificGrowthRate = 'log2';

% define time binning parameters
specificBinning = 10; % 10 min time bins
binsPerHour = 60/specificBinning;

% load measured experiment data
%experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date);
%cd(experimentFolder)
filename = strcat('glycogen-',date,'-width1p4-jiggle-0p5.mat');
load(filename,'D5');


% assemble data matrix from all xy positions
xy_start = 1;
xy_end = 16;
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end);


% isolate volume (Va), drop, track number, and fluorescent label data
volumes = glycogen_data(:,4);        % col 4 = calculated va_vals (cubic um)
isDrop = glycogen_data(:,2);         % col 2 = isDrop, 1 marks a birth event
trackNum = glycogen_data(:,9);       % col 9 = track number (not ID from particle tracking)
isYFP = glycogen_data(:,10);         % col 10 = 1 if above YFP threshold
isCFP = glycogen_data(:,11);         % col 11 = 1 if above CFP threshold


% calculate growth rate
growthRates = calculateGrowthRate_glycogen(volumes,isDrop,trackNum);



% truncate data to non-erroneous (e.g. bubbles) timestamps
maxTime = 8; % in hours
frame = glycogen_data(:,8);      % col 8 = frame in image sequence
timeInSeconds = frame * 120;     % 2 min per frame
timeInHours = timeInSeconds/3600;


if maxTime > 0
    glycogenData_bubbleTrimmed = glycogen_data(timeInHours <= maxTime,:);
    growthRates_bubbleTrimmed = growthRates(timeInHours <= maxTime,:);
else
    glycogenData_bubbleTrimmed = glycogen_data;
    growthRates_bubbleTrimmed = growthRates;
end
clear maxTime frame timeInSeconds timeInHours


% bin growth rate into time bins based on timestamp
frame = glycogenData_bubbleTrimmed(:,8);      % col 8 = frame in image sequence
timeInSeconds = frame * 120;     % 2 min per frame
timeInHours = timeInSeconds/3600;
bins = ceil(timeInHours*binsPerHour);


% isolate selected specific growth rate and remove nans from data analysis
specificColumn = 3;
xmin = -0.5;
xmax = 2.5;
growthRate_log2 = growthRates_bubbleTrimmed(:,specificColumn);

growthRt_noNaNs = growthRate_log2(~isnan(growthRate_log2),:);
bins_noNaNs = bins(~isnan(growthRate_log2),:);
isYFP_noNaNs = isYFP(~isnan(growthRate_log2),:);
isCFP_noNaNs = isCFP(~isnan(growthRate_log2),:);


% separate YFP and CFP populations
isBoth = isYFP_noNaNs+isCFP_noNaNs;
growthRt_final = growthRt_noNaNs(isBoth < 2);
bins_final = bins_noNaNs(isBoth < 2);
isYFP_final = isYFP_noNaNs(isBoth < 2);
isCFP_final = isCFP_noNaNs(isBoth < 2);

growthRt_yfp = growthRt_final(isYFP_final == 1);
growthRt_cfp = growthRt_final(isCFP_final == 1);



% bin growth rate by time
bins_yfp = bins_final(isYFP_final == 1);
bins_cfp = bins_final(isCFP_final == 1);

binned_yfp = accumarray(bins_yfp,growthRt_yfp,[],@(x) {x});
binned_cfp = accumarray(bins_cfp,growthRt_cfp,[],@(x) {x});


% calculate mean, standard dev, counts, and standard error
y_bin_means = cellfun(@mean,binned_yfp);
y_bin_stds = cellfun(@std,binned_yfp);
y_bin_counts = cellfun(@length,binned_yfp);
y_bin_sems = y_bin_stds./sqrt(y_bin_counts);

c_bin_means = cellfun(@mean,binned_cfp);
c_bin_stds = cellfun(@std,binned_cfp);
c_bin_counts = cellfun(@length,binned_cfp);
c_bin_sems = c_bin_stds./sqrt(c_bin_counts);



% plot growth rate over time
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
axis([0,max(timeInHours)+0.1,xmin,xmax])
xlabel('Time (hr)')
ylabel('Growth rate')
title(strcat(date,': (',specificGrowthRate,')'))
legend('YFP mutant', 'CFP WT')

