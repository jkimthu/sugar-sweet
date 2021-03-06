% dynamicOutlines - glycogen - cfp

% Goal: this version of dynamicOutlines displays colored outlines of 
%       tracked cells on cfp images


% Strategy:
%
%     0. initialize tracking and image data
%     1. isolate ellipse data from movie (stage xy) of interest
%     2. identify tracks in rejects. all others are tracked.
%     3. for each image, initialize current image
%            4. define major axes, centroids, angles
%            5. draw ellipses from image, color based on presence of tracked cell (no color designation yet)                
%           6. display and save
%     7. woohoo!


% last edit: jen, 2019 Feb 10
% commit: cfp tracking of 2019-02-06 dataset, xy1, with YFP and CFP colored


% OK LEZ GO!

%% A. Initialize experiment data

clc
clear


% 0. initialize data
date = '2019-02-06';
%cd(strcat('D:\',date))
load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');

% 0. initialize channel and xy movie to analyze
channel = 'c2';  % phase = c1; c2 = CFP; c3 = YFP)

% 0. initialize threshold intensity value
threshold = 103.4; % determined to not ID both fluorophores in a single cell

%%
% 0. initiate loop through movies of interest
for xy = 1
    
    if xy >= 10
        xy_nomen = strcat('xy',num2str(xy));
    else
        xy_nomen = strcat('xy0',num2str(xy));
    end
    
    
    % 1. compile experiment data matrix
    dt_min = 3;
    xy_start = xy;
    xy_end = xy;
    xyData = buildDM_glycogen(D5, xy_start, xy_end,dt_min);
    clear xy_start xy_end
    
    
    % 2. initialize image data
    conversionFactor = 6.5/60;      % scope5 Andor COSMOS = 6.5um pixels / 60x magnification
    
    
    % 3. create directory of image names in chronological order
    %imgDirectory = dir(strcat('glycogen-combo1-',date,'t*.tif'));
    imgDirectory = dir(strcat('*',xy_nomen,channel,'.tif'));
    names = {imgDirectory.name};
    
    
    % 4. identify tracks present in each frame
    totalFrames = length(imgDirectory); % total frame number
    trackNums = [];
    for fr = 1:totalFrames
        
        tracksInCurrentFrame = xyData(xyData(:,9) == fr,:);      % col 9  = frame # in buildDM_glycogen
        trackNums{fr,1} = tracksInCurrentFrame(:,12);            % col 12 = TrackNum, which differs from orginial dO script (which uses trackID)
        
    end
    clear fr tracksInCurrentFrame totalFrames
    
   
    
    % 5. overlay colored cell outlines over each image file
    for img = 1:length(names)
        
        % i. initialize current image
        cla
        I=imread(names{img});
        filename = strcat('dynamicOutlines-glycogen-xy',num2str(xy),'-cfp-frame',num2str(img),'.tif');
        
        figure(1)
        imshow(I, 'DisplayRange',[100 250]); % 2019-02-06
        % imtool(I), displays image in grayscale with range
        % lowering right # increases num sat'd pxls
        
        
        % ii. if no particles to display, save and skip
        if isempty(trackNums{img}) == 1
            saveas(gcf,filename)
            continue
            
        else
            % iii. else when tracked lineages are present, isolate data for each image
            dm_currentImage = xyData(xyData(:,9) == img,:);    % col 9 = frame #
            
            isCFP = dm_currentImage(:,13)> threshold;          % col 13 = CFP (above threshold? yes/no)
            isYFP = dm_currentImage(:,14)> threshold;          % col 14 = YFP (above threshold? yes/no)
            
            isBoth = isCFP + isYFP;
            
            majorAxes = dm_currentImage(:,2);       % col 2 = lengths
            minorAxes = dm_currentImage(:,4);       % col 4 = widths
            
            centroid_X = dm_currentImage(:,7);      % col 7 = x coordinate of centroid
            centroid_Y = dm_currentImage(:,8);      % col 8 = y coordinate of centroid
            angles = dm_currentImage(:,11);         % col 11 = angle of rotation of fit ellipses
            
            trackNum = dm_currentImage(:,12);       % col 12 = track #
            
            
            
            % iv. for each tracked particle in current image, draw ellipse 
            
            for particle = 1:length(trackNum)
                
                [x_rotated, y_rotated] = drawEllipse(particle,majorAxes, minorAxes, centroid_X, centroid_Y, angles, conversionFactor);
                lineVal = 0.5;
                
                
                % color ellipse based on fluorescent signal
                if isBoth(particle) == 1
                    color = rgb('Crimson');
                elseif isYFP(particle) == 1
                    color = rgb('GoldenRod');
                elseif isCFP(particle) == 1
                    color = rgb('Cyan');
                else
                    color = rgb('Indigo');
                end
                
                
                hold on
                plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
                text((centroid_X(particle)+2)/conversionFactor, (centroid_Y(particle)+2)/conversionFactor, num2str(trackNum(particle)),'Color',color,'FontSize',10);
                xlim([0 2048]);
                ylim([0 2048]);
                
            end
            
            
        end
        title(num2str(img))
        
        % 12. save
        saveas(gcf,filename)
        
    end
    
    
    
end



