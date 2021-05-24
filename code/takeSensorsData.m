close all; clear all; clc;
%% Create server


%% Taking data

%% Save data
save('dataset-12-5-18.mat','data','states','subjects')
%% Delete server
clear m
connector off
%% Reference
% https://blogs.mathworks.com/community/2014/10/06/acquire-data-from-device-sensors-with-matlab-mobile/
q = 1;
figure(1);
subplot(2,1,1);
plot(data(q).Acceleration_time,data(q).Acceleration);hold on; grid on;
xlabel('time[s]');
title('Acceleration Sensor')
subplot(2,1,2);
plot(data(q).AngularVelocity_time,data(q).AngularVelocity);hold on; grid on;
xlabel('time[s]');
title('Angular Velocity Sensor');
