% calculateGrowthRate_glycogen


% goal: this function centralizes the calculation of 4 plausible growth
%       rates. the only growth rate not calculated here is mu.


% main difference from original function: no time or curve input. dt is set
% at 120 sec (2 min), the imaging frequency.



% last updated: jen, 2018 October 10

% commit: edit for use in first glycogen analysis, 2018-10-03 data


% Go go let's go!

%%
function [growthRates] = calculateGrowthRate_glycogen(volumes,isDrop,trackNum)

% input data:
%        volumes     =  calculated va_vals (cubic um)
%        isDrop      =  1 marks a birth event, 0 is normal growth


% 0. define dt
dt = 120; % timestep in seconds


% 1. calculate dVdt
dV_noNan = diff(volumes);
dV = [NaN; dV_noNan];
dVdt_raw = dV/dt * 3600;                % final units = cubic um/hr



% 2. calculate dVdt_norm (normalized by initial volume)
dV_norm = [NaN; dV_noNan./volumes(1:end-1)];
dVdt_norm = dV_norm/dt * 3600;          % final units = 1/hr   
        

                
% 3. calculate dVdt_log = d(log V)/dt
dV_log_noNan = diff(log(volumes));
dV_log = [NaN; dV_log_noNan];
dVdt_log = dV_log/dt * 3600;           % final units = cubic um/hr
dVdt_log2 = dVdt_log/log(2);


% 4. calculate dVdt_lognorm = d(log V)/dt normalized by initial volume
dV_lognorm = [NaN; dV_log_noNan./volumes(1:end-1)];
dVdt_lognorm = dV_lognorm/dt * 3600;         % final units = 1/hr



% 5. replace all growth rates at division events with NaN
growthRates = [dVdt_raw, dVdt_norm, dVdt_log2, dVdt_lognorm];
growthRates(isDrop == 1,:) = NaN;



% 6. replace all growth rates at track transitions with NaN
isTransition = [NaN; diff(trackNum)];
growthRates(isTransition > 0,:) = NaN; % doesn't make a different if data comes only from full curves
        

% 7. output array with all growth rates, of columns in following order:
%     (i) dVdt_raw; (ii) dVdt_norm; (iii) dVdt_log2; (iv) dVdt_lognorm
end







