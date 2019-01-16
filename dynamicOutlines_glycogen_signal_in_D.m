% dynamicOutlines - glycogen - signal in D

% Goal: this version of dynamicOutlines displays colored outlines of 
%       all tracked cells (pre-quality control) on phase images


% Strategy:
%
%     0. initialize tracking and image data
%     1. isolate ellipse data from movie (stage xy) of interest
%     2. identify tracks in rejects. all others are tracked.
%     3. for each image, initialize current image
%            4. define major axes, centroids, angles, CFP and YFP signal
%            5. draw ellipses from image, color based on presence of signal
%
%                  CFP = 1; YFP = 0           = Cyan
%                  CFP = 0; YFP = 1;          = Gold
%                  CFP = 1; YFP = 1;          = Red
%                  CFP = 0; YFP = 0;          = Purple
%                      
%           6. display and save
%     7. woohoo!


% last edit: jen, 2019 January 16
% commit: plots all particles found through particle tracking, D, over CFP
%         channel


% OK LEZ GO!

%% A. Initialize experiment data

clc
clear


% 0. initialize data
date = '2018-11-23';
cd(strcat('D:\',date))
load(strcat('glycogen-',date,'-width1p4-jiggle-0p5.mat'),'D');


% 0. initialize channel and xy movie to analyze
channel = 'c2'; %c1 = phase; c2 = CFP; c3 = YFP

for xy = 16
    
    if xy >= 10
        xy_nomen = strcat('xy',num2str(xy));
    else
        xy_nomen = strcat('xy0',num2str(xy));
    end
    
    
    % 1. compile experiment data matrix
    xy_start = xy;
    xy_end = xy;
    xyData = buildDM_glycogen(D, xy_start, xy_end);
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
    
    
    %% B. Determine fluorescent signal profiles per particle, per image
    
    
    % 5. overlay colored cell outlines over each image file
    for img = 1:length(names)
        
        % i. initialize current image
        cla
        I=imread(names{img});
        filename = strcat('dynamicOutlines-glycogen-xy',num2str(xy),'-frame',num2str(img),'-',channel,'.tif');
        
        figure(1)
        imshow(I, 'DisplayRange',[100 150]); % 2018-11-23
        % imtool(I), displays image in grayscale with range
        % lowering right # increases num sat'd pxls
        
        
        % ii. if no particles to display, save and skip
        if isempty(trackNums{img}) == 1
            saveas(gcf,filename)
            continue
            
        else
            % iii. else when tracked lineages are present, isolate data for each image
            dm_currentImage = xyData(xyData(:,9) == img,:);    % col 9 = frame #
            
            isCFP = dm_currentImage(:,13);          % col 13 = CFP (above threshold? yes/no)
            isYFP = dm_currentImage(:,14);          % col 14 = YFP (above threshold? yes/no)
            signalSum = isCFP + isYFP;
            
            majorAxes = dm_currentImage(:,2);       % col 2 = lengths
            minorAxes = dm_currentImage(:,4);       % col 4 = widths
            
            centroid_X = dm_currentImage(:,7);      % col 7 = x coordinate of centroid
            centroid_Y = dm_currentImage(:,8);      % col 8 = y coordinate of centroid
            angles = dm_currentImage(:,11);         % col 11 = angle of rotation of fit ellipses
            
            trackNum = dm_currentImage(:,12);       % col 12 = track #
            
            
            
            % iv. for each particle of interest in current image,
            %     draw ellipse colored based on presence of fluorescent labels
            %
            %                  CFP = 1; YFP = 0    = Cyan
            %                  CFP = 0; YFP = 1    = Gold
            %                  Both signals        = Red
            %                  Neither             = Purple
            
            
            for particle = 1:length(trackNum)
                
                [x_rotated, y_rotated] = drawEllipse(particle,majorAxes, minorAxes, centroid_X, centroid_Y, angles, conversionFactor);
                lineVal = 0.5;
                
                % if both signals > threshold
                if signalSum(particle) == 2
                    color = rgb('Crimson');
                    
                    
                    % if neither signal > threshold
                elseif signalSum(particle) == 0
                    color = rgb('Indigo');
                    
                    
                    % if only CFP > threshold
                elseif isCFP(particle) == 1
                    color = rgb('Cyan');
                    
                    
                    % if only YFP > threshold
                elseif isYFP(particle) == 1
                    color = rgb('GoldenRod');
                    
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



