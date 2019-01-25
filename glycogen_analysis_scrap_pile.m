% glycogen analysis scrap pile

% Goal: this file contains pieces of test code, used when debugging various
%       elements of glycogen data analysis, which aims to determine growth
%       rates of WT and glycogen breakdown mutants under a pulsling nutrient environment.
%
%       input data is a timelapse with three channels: phase, YFP and CFP
%       stratgey:
%                0. initialize experiment parameters
%                1. particle identification, characterization and tracking
%                2. quality control: clean dataset prior to growth rate calculations
%                3. fun part! compute population growth rates between YFP and CFP populations


% BIIIIG thank yous to Cherry Gao for helping me get started!

% last update: Jen, 2019 January 25
% commit: rename (previously glycogen_analysis_windows) and edit comments for sharing

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
xy_final = 2; % total num of xy positions in analysis
num_tpt = 241; % total num of timepoints in analysis


% 1. particle identification, characterization and tracking
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
xy = 2;
tpt = 1;

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
figure(1); imshow(imadjust(ph_image)); title('phase');         % "imadjust" automatically adjusts brightness of 16-bit image
figure(2); imshow(imadjust(c_image)); title('cfp');
figure(3); imshow(imadjust(y_image)); title('yfp');

%% continue test analysis

% edge detection of cells in phase image
close all
bw = edge(ph_image,'sobel');

bg = imgaussfilt(ph_image,0.8);
bw_g = edge(bg,'sobel');

figure(1); imshowpair(bw,bw_g,'montage'); title('step 1: edge detection');

bw1 = edge(bg,'sobel');
figure(2); imshow(bw1); title('edge, smoothed')

%% comparing edge detection methods

bw1 = edge(bw_g,'sobel');

bw2 = edge(bw_g,'Prewitt');

figure(2); imshowpair(bw1,bw2,'montage')

bw3 = bw1+bw2;
figure(3); imshow(bw3)

bw = bw1;

%% adjusting edge detection for fullness

% clean up edge detection
se = strel('disk',1); 
%se = strel('disk',3); % structuring element; disk-shaped with radius of 3 pixels
bw_dil = imdilate(bw, se); % dilate
figure(3); imshow(bw_dil); title('step 2: dilate edges');

bw_fill = imfill(bw_dil, 'holes'); % fill
figure(4); imshow(bw_fill); title('step 3: fill');

bw_final = imerode(imerode(bw_fill,se),se); % erode twice to smooth; this is your final mask
figure(5);  imshow(bw_final); title('step 4: erode or smooth');

figure(6); imshowpair(imadjust(ph_image),bw_final); title('final mask-phase overlay');

%% cool cherry trick: linking subplots
% visualize effect of each step of edge detection and cleanup
figure;
ax(1) = subplot(231); imshowpair(imadjust(ph_image),bw_final); title('final mask-phase overlay');
ax(2) = subplot(232); imshow(bw); title('step 1: edge detection');
ax(3) = subplot(233); imshow(bw_dil); title('step 2: dilate edges');
ax(4) = subplot(234); imshow(bw_fill); title('step 3: fill');
ax(5) = subplot(235); imshow(bw_final); title('step 4: erode or smooth');
linkaxes(ax(:),'xy'); % link subplots so that you can zoom-in on all subplots simultaneously

%% segmentation and particle property measurements

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

%% measuring fluorescence signal from phase mask overlays

% 2. determine which cells are YFP or CFP labeled
% overlay mask with YFP fluorescence image by dot-product
y_mask_overlay = bw_final .* double(y_image); % convert your uint16 image to double class
c_mask_overlay = bw_final .* double(c_image);

% compute the fluorescence intensity value for each cell
for c = 1:cc.NumObjects
    
    pixel_id_of_cell = []; % initialize
    pixel_id_of_cell = cc.PixelIdxList{c}; % pixel index where the cell of interest is
    
    y_fluo_int_cells(c) = mean(y_mask_overlay(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
    c_fluo_int_cells(c) = mean(c_mask_overlay(pixel_id_of_cell));
    

    % assign mean fluorescence intensity to cell
    stats(c).yfp = y_fluo_int_cells(c); % yfp
    stats(c).cfp = c_fluo_int_cells(c); % cfp
    
end

% look at the distribution of YFP intensities of cells to determine a threshold
figure; histogram(y_fluo_int_cells,50); title('yfp'); xlabel('mean intensity')
figure; histogram(c_fluo_int_cells,50); title('cfp'); xlabel('mean intensity')

%% compiling data structure

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
yfp = extractfield(stats,'yfp')';
cfp = extractfield(stats,'cfp')';
p_unit.yfp = yfp;
p_unit.cfp = cfp;
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
        yfp = extractfield(stats,'yfp')';
        cfp = extractfield(stats,'cfp')';
        p_unit.yfp = yfp;
        p_unit.cfp = cfp;
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


%% YFP and CFP thresholds: my constructed method, from observations of reduced freq xy02 data

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

%% YFP-based particle sets, image 9

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

%% CFP-based particle sets, image 9

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

%% YFP-based particle sets, image 17

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

%% CFP-based particle sets, image 17

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

%% save figures

for f = 1:8
    figure(f)
    plotName = strcat('bar-figure',num2str(f));
    saveas(gcf,plotName,'epsc')
    close(gcf)
end



%% 3. compute population growth rates between YFP and CFP populations

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


% 4. truncate data to non-erroneous (e.g. bubbles) timestamps
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



% 9. test threshold, throw error if any points are identified as both
isBoth = isCFP+isYFP;
if sum(isBoth == 2) > 0
    error('threshold fail! some cells positive for both fluorophores')
end

% this step doesn't matter as long as threshold doesn't allow double-positives
growthRt_final = growthRt_noNaNs(isBoth < 2);
glycogenData_final = glycogenData_noNaNs(isBoth < 2,:);
bins_final = bins_noNaNs(isBoth < 2);
isYFP_final = isYFP(isBoth < 2);
isCFP_final = isCFP(isBoth < 2);



% 10. separate growth rates by fluorophore

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


% 11. bin growth rate by time
binned_yfp = accumarray(bins_yfp,growthRt_yfp,[],@(x) {x});
binned_cfp = accumarray(bins_cfp,growthRt_cfp,[],@(x) {x});


% 12. calculate mean, standard dev, counts, and standard error
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
axis([0,8.5,-1,1])
xlabel('Time (hr)')
ylabel('Growth rate')
title(strcat(date,': (',specificGrowthRate,')'))
legend('YFP WT', 'CFP mutant')




%% 4. compute population growth rates as change in total area

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
specificBinning = 30; % in minutes
binsPerHour = 60/specificBinning;


% 0. define fluorescence intensity threshold
threshold = 103.4;


% 1. assemble data matrix
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
clear xy_start xy_end


% 2. truncate data to non-erroneous (e.g. bubbles) timestamps
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



% 4. isolate YFP and CFP intensities, area and time
cfp = glycogenData_bubbleTrimmed(:,13);         % col 13 = mean CFP intensity
yfp = glycogenData_bubbleTrimmed(:,14);         % col 14 = mean YFP intensity
area = glycogenData_bubbleTrimmed(:,6);         % col 6 = measured surface area



% 5. convert intensities to (+) or (-) fluorophore
isCFP = cfp > threshold;
isYFP = yfp > threshold;



% 6. test threshold, throw error if any points are identified as both
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


% 8. bin growth rate by time
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


% 11. normalize by initial sum
y_norm2 = y_norm/y_norm(1);
c_norm2 = c_norm/c_norm(1);


% plot growth rate over time
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




