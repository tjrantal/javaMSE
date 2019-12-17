
fclose all;
close all;
clear all;
clc;
addpath('functions');


javaaddpath('../..//build/libs/javaMSE-1.0.0.jar');	%java class for orientation 
sRate = 100;	%Resample data to this sample rate
dataPath = 'data/';
sensorPaths = {'imuLog'};
imuFileName = 'IMU_192630002994_2019-12-17_111559.txt';	%Change this manually to correspond to your file

%Read data
imuData = readLog([dataPath sensorPaths{1} '/' imuFileName]);

%Sort the data based on sensor stamp and package index (data sampled with
%Movesense sensor and saved in separate threads in Android -> order can get
%jumbled)
[ignore sortOrder] = sortrows(imuData.data(:,2:3),[1 2]);
imuData.data = imuData.data(sortOrder,:);

%Check for missing data
uStamps = unique(imuData.data(:,2));
figure,plot(diff(uStamps));

if max(diff(uStamps)) < 80
   
else
    disp('Missing data found!');
    keyboard;
end



%Assume constant sample rate at 100 Hz (requested from the sensor)

aRes = sqrt(sum(imuData.data(:,4:6).^2,2))./9.81; %Divide by 9.81 to get acceleration in g

%Look for the walking bout
figure,plot(aRes)   %Not how the calibration is a bit off, standing still should result in 1 g

[buffered ignore] = buffer(aRes,5*sRate,0,'nodelay');
[bInd ignore] = buffer(1:length(aRes),5*sRate,0,'nodelay');


sds = std(buffered,[],1);
sdInd = bInd(1,:);
hold on;
plot(sdInd,sds,'r*');

gaitPeaks = getPeaks(sds,0.1);
pLengths = cellfun(@(x) length(x),gaitPeaks);
[ignore gaitBout] = max(pLengths);
gaitIndices = sdInd(gaitPeaks{gaitBout}(1))+2.5*sRate:(sdInd(gaitPeaks{gaitBout}(end))+2.5*sRate-1);
gaitResultant = aRes(gaitIndices);
plot(gaitIndices,gaitResultant,'g');


%RCME
 
if length(gaitResultant) < 60*100
    gaitbuffered = gaitResultant;
    indbuffered = gaitIndices;
else
    [gaitbuffered, discard] = buffer(gaitResultant,60*sRate,60*sRate/2,'nodelay'); %50% overlap
    [indbuffered, discard] = buffer(gaitIndices,60*sRate,60*sRate/2,'nodelay'); %50% overlap
end

includedEpochCnt = 0;
result = struct();
for ep = 1:size(gaitbuffered,2)
    includedEpochCnt = includedEpochCnt +1;
    result(ep).values = getBoutAnalysisJava(gaitbuffered(:,ep),[0.7 1.52]./2);
    plot(indbuffered(:,ep),gaitbuffered(:,ep));
end

