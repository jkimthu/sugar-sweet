% glycogen_analysis

% Goal: builds a dataset from particle tracking WT and glycogen breakdown
%       mutants under a pulsing nutrient environment.
%
%       input data is a timelapse with three channels: phase, CFP and YFP


% Stratgey:
%                0. initialize experiment parameters
%                1. particle identification, characterization and tracking
%                2. quality control: clean dataset prior to growth rate calculations


% Acknowledgements: Cherry Gao, Vicente Fernandez and Jeff Guasto!


% last update: Jen, 2019 Feb 26
% commit: first analysis of 2019-02-25 experiment


% ok let's go!

%% 0. initialize experiment parameters

clear
clc

%  initialize conversion from pixels to um based on camera and magnification
%  scope5 Andor COSMOS = 6.5um pixels / 60x magnification
ConversionFactor = 6.5/60; % units = um/pixel


% initialize experiment folder
experiment = '2019-02-25';


% open folder for experiment of interest
%imageFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',experiment,'_xy02_full'); % _xy02 for folder
%imageFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',experiment,'_xy02'); % _xy02 for folder
%cd(imageFolder); % move to the folder containing images 


% initialize image names
imageName = strcat('glycogen-',experiment);


% initialize total number of movies
xy_final = 30; % total num of xy positions in analysis


% initialize array of desired timepoints
selected_tpt = linspace(1,161,161);


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


% for each xy, loop through all timepoints and build data structure
for xy = 1:xy_final
    
    % initialize xy name
    if xy < 10
        xyName = strcat('xy0',num2str(xy));
    elseif xy < 100
        xyName = strcat('xy',num2str(xy));
    end
    
    
    for i = 1:length(selected_tpt)
        
        tpt = selected_tpt(i);
        
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
        
        
        % gaussian smoothing of phase image
        ph_smoothed = imgaussfilt(ph_image,0.8);
        
        
        % edge detection of cells in phase image
        bw = edge(ph_smoothed,'sobel');
        
        
        % clean up edge detection
        se = strel('disk',1); % structuring element; disk-shaped with radius of 3 pixels
        bw_dil = imdilate(bw, se); % dilate
        bw_fill = imfill(bw_dil, 'holes'); % fill
        bw_final = imerode(imerode(bw_fill,se),se); % erode twice to smooth; this is your final mask
        
        
        % segment cells
        cc = bwconncomp(bw_final);
        
        
        % gather properties for each identified particle
        stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
      
        clear bw bw_dil bw_fill
        
        
        
        % 2. determine which cells are YFP or CFP labeled        
        % overlay mask with YFP fluorescence image by dot-product
        y_mask_overlay = bw_final .* double(y_image); % convert your uint16 image to double class
        c_mask_overlay = bw_final .* double(c_image);
        
        % compute the fluorescence intensity value for each cell
        for c = 1:cc.NumObjects;
            
            pixel_id_of_cell = []; % initialize
            pixel_id_of_cell = cc.PixelIdxList{c}; % pixel index where the cell of interest is
            
            y_fluo_int_cells(c) = mean(y_mask_overlay(pixel_id_of_cell)); % compute the mean intensity of pixels in the cell
            c_fluo_int_cells(c) = mean(c_mask_overlay(pixel_id_of_cell));
            
        end
        clear bw_final
        


        % add mean fluorescent intensity of each cell to data structure
        for c = 1:cc.NumObjects;
            
            % yfp
            stats(c).yfp = y_fluo_int_cells(c);
            
            % cfp
            stats(c).cfp = c_fluo_int_cells(c);
           
        end
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
        YFP = extractfield(stats,'yfp')';
        CFP = extractfield(stats,'cfp')';
        p_unit.yfp = YFP;
        p_unit.cfp = CFP;
        clear YFP CFP
        
        % frame
        p_unit.Frame = i;
        
        
        % eccentricity and angle
        ecc = extractfield(stats,'Eccentricity')';
        angle = extractfield(stats,'Orientation')';
        p_unit.Ecc = ecc;
        p_unit.Angle = angle;
        
        p_clone(i) = p_unit;

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
    LowerBound = 0.7;     % lower bound for restricted field, or -Inf
    UpperBound = 1.4;     % upper bound for all conditions

    % to actually trim the set:
    p_clone_trim2 = ParticleTrim_glycogen(p_clone_trim1,TrimField,LowerBound,UpperBound);
    
    
    % remove rows without data
    
    % v. track remaining particles based on coordinate distance
    
    TrackMode = 'position';       % Choice of {position, velocity, acceleration} to predict position based on previous behavior
    DistanceLimit = 5;            % Limit of distance a particle can travel between frames, in units defined by ConversionFactor
    MatchMethod = 'best';         % Choice of {best, single}
    p_Tracks = Particle_Track_glycogen(p_clone_trim2,TrackMode,DistanceLimit,MatchMethod);
    
    
    
    % vi. link tracks together and store in data matrix!
    
    PT = TrackLinker(p_Tracks, 'position', 'position', 5, 3, 2);
    TL=TrackLength(PT);
    PT(TL<8)=[];
    A=arrayfun(@(Q) max(Q.Area),PT);

    
    % 4. Store tracked data into single workspace, D
    D{xy} = PT;
    %D{xy} = p_Tracks;
    
    disp(strcat('analyzed xy(',num2str(xy),') of (',num2str(xy_final),')!'))
    
    
end

save(strcat('glycogen-',experiment,'-allXYs.mat'),'D')  


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
experiment = '2019-02-25';

% 0. open folder for experiment of interest
%newFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',experiment,'_xy02_full')
%newFolder = strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',experiment,'_xy02');%,'  (t300)');
%cd(newFolder);

load(strcat('glycogen-',experiment,'-allXYs.mat'));

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
                                                                           
for n = 1:length(D);                                                       
    
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
    
    %     % when all tracks finished,
    %     if jump_counter >= 1
    %
    %         % 7. save accmulated rejects
    %         rejectD{criteria_counter,n} = trackScraps;
    %
    %         % 8. and insert data after jump at end of data matrix, D2
    %         D2{n} = [data; trackScraps];
    %
    %     else
    D2{n} = data;
    %     end
    
    % 9. repeat for all movies
    clear trackScraps X data;
    
end
    
clear  jumpFrac jump_counter Rates m n;



%  CRITERIA TWO
%  tracks must be at least window size in length (5 frames)
criteria_counter = criteria_counter + 1;

% 0. initialize new dataset before trimming
D3 = D2;

for n = 1:length(D);
    
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

for n = 1:length(D);   
    
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

save(strcat('glycogen-',experiment,'-allXYs-jiggle-0p5.mat'), 'D', 'D2', 'D3', 'D4', 'D5', 'rejectD')
disp('Quality control: complete!')




