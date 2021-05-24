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


%% Train
inputSize = 14;
outputSize = 100;
outputMode = 'last';

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(outputSize,'OutputMode',outputMode)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

maxEpochs = 150;
miniBatchSize = 27;
options = trainingOptions('sgdm', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize);

net = trainNetwork(num2cell(Xtrain,1)',(Ytrain),layers,options);

%% Test
YPred = classify(net,num2cell(Xtest,1)');
acc = sum((YPred) == (Ytest))./numel(Ytest)
[yt,yp]=transYtYp(Ytest,YPred,numClasses);
figure(1)
plotconfusion(yt,yp);
%% PCA
figure();
gplotmatrix(Xdata',[],Ydata); grid on;


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