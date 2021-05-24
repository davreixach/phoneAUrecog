function [data] = test(m,state,time,subject,STATES_VEC,SUBJ_VEC)
%TEST performs test for specified state, time and subject
%   [data] = test(state,time,subject)
%
% Authors:
% Bermeo & Reixach, ETSEIB(UPC) 2018

m.Logging = 1; % Take data
pause(1);
m.Logging = 0; 
discardlogs(m);

m.Logging = 1; % Take data
for i=1:time
    pause(1)
    fprintf('%i ',time-i+1)
end
fprintf('\n')
m.Logging = 0;

[accel,accel_time]=accellog(m); % data
[angvel,angve_time]=angvellog(m);

data.Acceleration         = accel;
data.Acceleration_time    = accel_time;
data.AngularVelocity      = angvel;
data.AngularVelocity_time = angve_time;
data.time_accel = linspace(0,time,length(accel_time));
data.time_angvel = linspace(0,time,length(angve_time));
data.type_test = STATES_VEC(state);
data.id_type_test = (state) ;
data.time = time;
data.subject = SUBJ_VEC(subject);

discardlogs(m)
end

