function ParticleOut=Particle_Track_glycogen(particles, trackmode, Rmax, matchmethod)

PlotFlag=0;
% Jeffrey Guasto
% 8/5/2011
% MIT

% last edit: Jen Nguyen 2019 Jan 22
% commit: minor capitalization change while troubleshooting particle tracking 


% "trackmode" = use nearest neighbors or predictive methods ['position', 'velocity', or 'acceleration']
% "Rmax" = maximum search radius to consider particle matches [0, inf]
% "matchmethod" = ['best', 'single']
%
% OUTPUT: ParticleTrack = [1.x, 2.y, 3.R_mean, 4.area, 5.ecc, 6.theta, 
%                              7.time, 8.number, 9.u_x, 10.u_y, 11.a_x, 12.a_y]
% INPUT:  
%   contains ParticleLocs: [1.x, 2.y, 3.R_mean, 4.area, 5.ecc, 6.theta, 7.time]

tic

% input/output
if nargin < 2  || isequal(trackmode,'default');     trackmode = 'position'; end
if nargin < 3  || isequal(Rmax,'default');          Rmax = inf; end
if nargin < 4  || isequal(matchmethod,'default');   matchmethod = 'best'; end

% pre-compute & initialize
N_frames = length(particles);
ParticleTracks = cell(2,1);
    
% loop over frames
h = waitbar(0,['Tracking Particles']);
for n = 1:N_frames
    if n/round(N_frames/100)==round(n/round(N_frames/100))
    waitbar(n/N_frames,h)
    end

    % find particles in current frame
    N_now = length(particles(n).X);
    % First loop: make cell array "ParticleTracks"; 
    % Discontinuous frames: deactivate all old tracks, start new tracks with all new particles
    
    % original from Jeff and Vicente
    %particleMat=[particles(n).X,particles(n).Y,particles(n).A,particles(n).MaxInt,ones(size(particles(n).X))*particles(n).Frame,particles(n).Ecc,particles(n).MajAx,particles(n).MinAx,particles(n).Ang];
    particleMat=[particles(n).X,particles(n).Y,particles(n).A,particles(n).MajAx,particles(n).MinAx,particles(n).yfp,particles(n).cfp,ones(size(particles(n).X)).*particles(n).Frame,particles(n).Ecc,particles(n).Angle];
    % 1 = X
    % 2 = Y
    % 3 = Area
    % 4 = Major axis
    % 5 = Minor axis
    % 6 = YFP
    % 7 = CFP
    % 8 = Frame
    % 9 = Eccentricity
    % 10 = Angle
    
    if n == 1; 
        ParticleTracks([1:N_now]) = num2cell(particleMat,2);    
        active = 1:N_now;
        continue        
    elseif particles(n).Frame-particles(n-1).Frame ~= 1
        N_tracks = length(ParticleTracks);
        ParticleTracks(N_tracks+[1:N_now])=num2cell(particleMat,2);
        active = (N_tracks+1):(N_tracks+N_now);
        continue  
    elseif N_now==0
        N_tracks = length(ParticleTracks);
        active=[];
        continue
    end


    % subsequent loops:    
    % estimate future position based on current position, velocity, and acceleration for active tracks      
    pos_est  = zeros(length(active),2);
    pos_past  = zeros(length(active),2);
    for m = 1:length(active)        
        pos_past(m,:)=ParticleTracks{active(m)}(end,1:2);
        N_past = length(ParticleTracks{active(m)}(:,1));
        if     N_past == 1 || isequal(trackmode,'position')
            pos_est(m,:) = ParticleTracks{active(m)}(end,1:2);
        elseif N_past == 2 || isequal(trackmode,'velocity')
            pos_est(m,:) = ParticleTracks{active(m)}(end,1:2) + ...
                          (ParticleTracks{active(m)}(end,1:2)-ParticleTracks{active(m)}(end-1,1:2));
        elseif  N_past > 2 && isequal(trackmode,'acceleration')
            pos_est(m,:) = ParticleTracks{active(m)}(end,1:2) + ...
                          (ParticleTracks{active(m)}(end,1:2)-ParticleTracks{active(m)}(end-1,1:2)) + ...
                           0.5*(ParticleTracks{active(m)}(end,1:2)-2*ParticleTracks{active(m)}(end-1,1:2)+ParticleTracks{active(m)}(end-2,1:2));
        end
    end

      
    % compute costs from prediction to next possible particles
    pos = [particles(n).X,particles(n).Y];
    [x_est, x] = meshgrid(pos_est(:,1), pos(:,1));
    [y_est, y] = meshgrid(pos_est(:,2), pos(:,2)); 
    cost = (x-x_est).^2 + (y-y_est).^2;  


    % compute distance from current position to next possible position
    [x_past, x] = meshgrid(pos_past(:,1), pos(:,1));
    [y_past, y] = meshgrid(pos_past(:,2), pos(:,2)); 
    Rng = (x-x_past).^2 + (y-y_past).^2; 
      
          
    % find best matches for estimates & possible future positions          
    if isequal(matchmethod,'best')
        [r_pos,c_est] = find(cost <= Rmax^2 & Rng <= Rmax^2);  
        ind_rc = sub2ind(size(Rng),r_pos,c_est);
        v = cost(ind_rc);
        matches = [];
        for jj = 1:length(active)
            ind_match = find(v == min(v));
            if length(ind_match)
                ind_match = ind_match(1);
                matches = [matches ; r_pos(ind_match) c_est(ind_match)];
                ind_del = find(r_pos == r_pos(ind_match) | c_est == c_est(ind_match)); %all row and column indices
                r_pos(ind_del) = []; c_est(ind_del) = []; v(ind_del) = [];
            end
        end
        if isempty(matches)
            N_tracks = sum(~cellfun(@isempty,ParticleTracks));
            ParticleTracks(N_tracks+[1:N_now])=num2cell(particleMat,2);
            active = (N_tracks+1):(N_tracks+N_now);     
            continue  
        end
        r_match = matches(:,1); c_match = matches(:,2);
        active_matched = active(matches(:,2));                  % matched estimated positions
        now_matched = matches(:,1);                             % matched new particles        
        active_nomatch = setdiff(1:length(active),matches(:,2));% tracks now inactive
        now_nomatch = setdiff(1:N_now,matches(:,1));            % unmatched new particles, start new tracks          
    elseif isequal(matchmethod,'single')
        [r_pos,c_est] = find(cost <= Rmax^2 & Rng <= Rmax^2);        
        [b, m1] = unique(c_est, 'first');
        [b, m2] = unique(c_est, 'last');
        c_match = find(ismember(c_Est, b(m1==m2)));   
        [b, m1] = unique(r_pos, 'first');
        [b, m2] = unique(r_pos, 'last');
        r_match = find(ismember(r_pos, b(m1==m2)));
        ind_match = intersect(c_match,r_match); 
        active_matched = active(c_est(ind_match));                  % matched tracks
        now_matched = r_pos(ind_match);                             % matched new particles        
        active_nomatch = setdiff(1:length(active),c_est(ind_match));% tracks now inactive
        now_nomatch = setdiff(1:N_now,r_pos(ind_match));           % unmatched new particles, start new tracks
    end
        
    % append matches to appropriate cells
    ParticleTracks(active_matched)=cellfun(@(x,y) cat(1,x,y),ParticleTracks(active_matched),num2cell(particleMat(now_matched,:),2),'UniformOutput',0);

    % deactivate particles
    active(active_nomatch) = []; % tracks now inactive
    
    % start new tracks
    N_tracks = sum(~cellfun(@isempty,ParticleTracks));

    if length(now_nomatch)
        active = [active (N_tracks+1):(N_tracks+length(now_nomatch))];
        ParticleTracks(N_tracks+1:N_tracks+length(now_nomatch))=num2cell(particleMat(now_nomatch,:),2);
    end
end % ii
  
  
    % pack tracks
    N_tracks = length(ParticleTracks);    
    tr_length=cellfun(@numel,ParticleTracks);
    [s,sInd]=sort(tr_length,'descend');
    ParticleTracks=ParticleTracks(sInd);
    ParticleTracks=cellfun(@(x,y) cat(2,x,ones(size(x,1),1)*y),ParticleTracks,num2cell([1:length(ParticleTracks)])','UniformOutput',0);
    
    % informational
    close(h);
    mytime = toc;
    disp(['  Number of Tracks: ', num2str(N_tracks)])
    disp(['  Avg. Track Length: ', num2str(mean(tr_length)),' frames']);
    disp(['  Elapsed Time: ', num2str(mytime), ' seconds'])
    disp('  ')
  
    ParticleOut=cellfun(@(A) struct('X',A(:,1),'Y',A(:,2),'Area',A(:,3),'MajAx',A(:,4),'MinAx',A(:,5),'YFP',A(:,6),'CFP',A(:,7),'Frame',A(:,8),'Ecc',A(:,9),'Angle',A(:,10)),ParticleTracks); 
end