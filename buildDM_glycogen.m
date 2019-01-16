% buildDM_glycogen (windows)

% goal: build a data matrix of track parameters (each column) over time
%       (rows) for all cells in a given xy position. option to specificy xy
%       positions and streamline data concatenation.

% adapted from original buildDM, but streamlined for glycogen proejct analysis

% last updated: jen, 2019 January 15

% commit: streamline original for glycogen analysis


function [dm] = buildDM_glycogen(D5,xy_start,xy_end)
%% initialize all values
  
tn_counter = 0;
dropThreshold = -0.75; % consider greater negatives a division event


Time = [];              % 2. Time
lengthVals = [];        % 3. lengthVals
isDrop = [];            % 4. isDrop
curveFinder = [];       % 5. curveFinder
widthVals = [];         % 6. widthVals
vaVals = [];            % 7. vaVals
surfaceArea = [];       % 8. surfaceArea
x_pos = [];             % 9. x coordinate of centroid
y_pos = [];             % 10. y coordinate of centroid
orig_frame = [];        % 11. orig_frame
stage_num = [];         % 12. stage_num
eccentricity = [];      % 13. eccentricity
angle = [];             % 14. angle of rotation of fit ellipse
trackNum = [];          % 15. trackNum  =  total track number (vs ID which is xy based)
CFP = [];               % 16. CFP
YFP = [];               % 17. YFP

%% loop through all xy positions and all tracks for data concatenation

for n = xy_start:xy_end % n = each inidividual xy position from experiment (movie)
    
    for m = 1:length(D5{n}) % m = each individual cell track from current movie
   
        %% frame number in original image
        frameTrack = D5{n}(m).Frame;
        orig_frame = [orig_frame; frameTrack];
        
        %% track number
        tn_counter = tn_counter + 1;
        tnTrack = ones(length(frameTrack),1)*tn_counter;
        trackNum = [trackNum; tnTrack];
        
        %% time
        timeTrack = (frameTrack-1)*2;  % experiment 2018-11-23 imaged every 2 min            
        Time = [Time; timeTrack];                                         
        clear frameTrack
        
        %% lengths
        lengthTrack = D5{n}(m).MajAx;%(7:lengthCurrentTrack+6);              % collect lengths (um)
        lengthVals = [lengthVals; lengthTrack];                            % concatenate lengths
        
        %% drop?
        dropTrack = diff(lengthTrack);
        trackDrops = dropTrack < dropThreshold;                                % converts different to a Boolean based on dropThreshold
        trackDrops = [0; trackDrops];                                              % * add zero to front, to even track lengths
        isDrop = [isDrop; trackDrops];
        
        
        %% widths
        widthTrack = D5{n}(m).MinAx;%(7:lengthCurrentTrack+6);               % collect widths (um)
        widthVals = [widthVals; widthTrack];                               % concatenate widths
        
        %% volume as a cylinder with hemispherical caps
        
        %v_cylinder = pi * lengthTrack .* (widthTrack/2).^2;                % approx. volume as a cylinder = pi * r^2 * h
        %v_ellipse = 4/3 * pi * lengthTrack/2 .* (widthTrack/2).^2;         % approx. volume as an ellipse
        vol_smallCylinder = pi * (widthTrack/2).^2 .* (lengthTrack - widthTrack);
        vol_sphere = 4/3 * pi * (widthTrack/2).^3;
        v_anupam = vol_smallCylinder + vol_sphere;                         % approx. volume as cylinder with spherical caps
        
        vaVals = [vaVals; v_anupam];
        
        clear v_ellipse v_cylinder vol_sphere vol_smallCylinder
        
        %% surface area
        sa_rectangle = (lengthTrack - widthTrack) .* widthTrack;
        sa_sphere = 4 * pi .* (widthTrack/2);
        sa_total = sa_rectangle + sa_sphere;
        
        surfaceArea = [surfaceArea; sa_total];
        
        clear sa_total sa_sphere sa_rectangle
        
        %% x positions in original image
        xTrack = D5{n}(m).X;%(7:lengthCurrentTrack+6);
        x_pos = [x_pos; xTrack];
        clear xTrack
        
        %% y positions in original image
        yTrack = D5{n}(m).Y;%(7:lengthCurrentTrack+6);
        y_pos = [y_pos; yTrack];
        clear yTrack
        
        %% eccentricity of ellipses used in particle tracking
        eccTrack = D5{n}(m).Ecc;%(7:lengthCurrentTrack+6);
        eccentricity = [eccentricity; eccTrack];
        clear eccTrack
        
        %% angle of ellipses used in particle tracking
        angTrack = D5{n}(m).Ang;%(7:lengthCurrentTrack+6);
        angle = [angle; angTrack];
        clear angTrack
        
        %% CFP and YFP
        cfpTrack = D5{n}(m).CFP;
        CFP = [CFP; cfpTrack];
        
        yfpTrack = D5{n}(m).CFP;
        YFP = [YFP; yfpTrack];
        
        clear cfpTrack yfpTrack
        
    end % for m
    
    disp(['Tracks (', num2str(m), ') assembled from movie (', num2str(n), ') !'])
    
end % for n


% compile data into single matrix
dm = [Time lengthVals isDrop curveFinder widthVals vaVals surfaceArea x_pos y_pos orig_frame stage_num eccentricity angle trackNum CFP YFP];
% 1. Time
% 2. lengthVals
% 3. isDrop
% 4. curveFinder
% 5. widthVals  
% 6. vaVals
% 7. surfaceArea
% 8. x_pos
% 9. y_pos
% 10. orig_frame
% 11. stage_num
% 12. eccentricity
% 13. angle
% 14. trackNum  =  total track number (vs ID which is xy based)
% 15. CFP
% 16. YFP

end