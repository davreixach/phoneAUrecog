clear all; close all; clc;
%% Load dataSet
load 'SPdata-1305.mat'
%% Data-preprocessing
numClasses = 7;
classes = [1:numClasses];
% David Data(1-4)
[featuresDataObj,Xdata,Ydata]=getFeacturesData(datad);
% Miguel Data(1-4)
[featuresDataObj1,Xdata1,Ydata1]=getFeacturesData(datam);
Xdata = [Xdata,Xdata1];
Ydata = [Ydata,Ydata1];
% David and Miguel Data(5-6)
[featuresDataObj1,Xdata1,Ydata1]=getFeacturesData(datae);
Xdata = [Xdata,Xdata1];
Ydata = [Ydata,Ydata1];

Xdata=dataNormalization(Xdata);
[Xtrain,Ytrain,Xtest,Ytest,Xdata1]=dividationDataSet(Xdata,Ydata,classes,0.6);

ytrain=transYtrain(Ytrain,numClasses);
ytest=transYtrain(Ytest,numClasses);


x = Xtrain;
t = ytrain;

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.

% Create a Pattern Recognition Network
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize, trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
percentErrors = sum(tind ~= yind)/numel(tind);

% View the Network
view(net)
view2(net)
%% Test
y = net(Xtest);
% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
figure, plotconfusion(ytest,y)
%figure, plotroc(t,y)

%% Appendix - Functions
function [featuresDataObj,Xdata,Ydata]=getFeacturesData(data)
    featuresDataObj(1).acc_means             = [0 0 0];  
    featuresDataObj(1).acc_standar_deviation = [0 0 0];
    featuresDataObj(1).acc_avg_res_means     = 0;
    featuresDataObj(1).acc_avg_res_stnds     = 0;
    featuresDataObj(1).anv_means             = [0 0 0];  
    featuresDataObj(1).anv_standar_deviation = [0 0 0];
    featuresDataObj(1).anv_avg_res_means     = 0;
    featuresDataObj(1).anv_avg_res_stnds     = 0;
    Xdata=[];
    N = length(data);
    for (i=1:N)
        featuresDataObj(i).acc_means              = mean(data(i).Acceleration);  
        featuresDataObj(i).acc_standard_deviation = std(data(i).Acceleration);
        featuresDataObj(i).acc_avg_res            = ...
            sqrt(sum(sum(data(i).Acceleration).^2))/(3*length(data(i).Acceleration));
        featuresDataObj(i).anv_means              = mean(data(i).AngularVelocity); 
        featuresDataObj(i).anv_standar_deviation  = std(data(i).AngularVelocity);
        featuresDataObj(i).anv_avg_res            =  ...
            sqrt(sum(sum(data(i).AngularVelocity).^2))/(3*length(data(i).AngularVelocity));
        featuresDataObj(i).y=data(i).id_type_test;
        Xdata(:,i)=[ featuresDataObj(i).acc_means';...
                            featuresDataObj(i).acc_standard_deviation';...
                            featuresDataObj(i).acc_avg_res ;...
                            featuresDataObj(i).anv_means' ;...
                            featuresDataObj(i).anv_standar_deviation' ;...
                            featuresDataObj(i).anv_avg_res ;...
                            ];
        Ydata(i)=data(i).id_type_test;
    end
end

function [Xtrain,Ytrain,Xtest,Ytest,Xdata1]=dividationDataSet(Xdata,Ydata,classes,porTrain)
        Xdata1=[];
        Xtrain=[];
        Ytrain=[];
        Xtest =[];
        Ytest =[];
        for (i =classes)
            Xdata1{i}=Xdata(:,Ydata==i);
            N = size(Xdata1{i},2);
            M = floor(N*porTrain);
            Xdata1{i}=Xdata1{i}(:,randperm(N));
            Xtrain=[Xtrain,Xdata1{i}(:,1:M)];
            Ytrain=[Ytrain;categorical(ones(M,1)*i)];
            Xtest=[Xtest,Xdata1{i}(:,(M+1):N)];
            Ytest=[Ytest;categorical(ones(N-M,1)*i)];
        end
        
end

function [yt,yp]=transYtYp(Ytest,YPred,numClasses)
    yt=zeros(numClasses,size(YPred,1));
    yp =zeros(numClasses,size(YPred,1));
    for i=1:size(YPred,1)
        yt(double(Ytest(i)),i)=1;
        yp(double(YPred(i)),i)=1;
    end
end

function [yt]=transYtrain(Ytrain,numClasses)
    yt=zeros(numClasses,size(Ytrain,1));
    for i=1:size(Ytrain,1)
        yt(double(Ytrain(i)),i)=1;
    end
end



function datan=dataNormalization(data)
        [M,N]=size(data);
        datan = [];
        for i=1:M
            max_ = max(data(i,:));
            min_ = min(data(i,:));
            datan(i,:) = (data(i,:)-min_)/(max_-min_);
        end
end

function plotData(data,nfigure);
    figure(nfigure)
    subplot(2,1,1);
    plot(data.Acceleration_time,data.Acceleration); grid on;
    xlabel('time[s]');
    ylabel('[m/s^2]');
    title('Acceleration Sensor');
    legend('x[m/s^2]','y[m/s^2]','z[m/s^2]');
    subplot(2,1,2);
    plot(data.AngularVelocity_time,data.AngularVelocity); grid on;
    xlabel('time[s]');
    ylabel('[rad/s]');
    title('Angular Velocity Sensor');
    legend('x[rad/s]','y[rad/s]','z[rad/s]')
end

function  outH = view2(net)
    jframe = view(net);
    %# create it in a MATLAB figure
    hFig = figure('Menubar','none', 'Position',[100 100 565 166]);
    jpanel = get(jframe,'ContentPane');
    [~,h] = javacomponent(jpanel);
    set(h, 'Parent', hFig, 'units','normalized', 'position',[0 0 1 1])
    %# close java window
    jframe.setVisible(false);
    jframe.dispose();
    % capture figure as a frame
    F = getframe(hFig);
    % close current figure
    close(hFig)
    % display the captured image data using imshow
    hFig = figure('Menubar','none', 'Position',[100 100 565 166]);
    axes('units','normalized', 'position',[0 0 1 1])
    imshow(F.cdata)
    axis off
    % return handle if needed
    if nargout == 1
        outH = hFig;
    end
end
