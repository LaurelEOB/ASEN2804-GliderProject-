clear; clc; close all;

%% ThrustTest_Single Summary
% First, some of the summary from Thrust:
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, process them, and fit the data into a standard
% formatting for output
%
% ThrustTest_Single is meant to be a testing function where you can test
% the conditioning steps that your group develops. This version of the code
% loads in only one data set at a time so that the workspace is less
% cluttered. It is suggested that students make plots of their thrust data
% often throughout their development in this code section to visually check
% that their conditioning is working as desired. Once conditioning is
% working on one set of data, students are engcouraged to try other single
% data sets. Once groups are satisfied, the full funciton has the same form
% as this, so the conditioning code can simply be dragged and dropped into
% the full "Thrust.m" funciton.
% 
% Note that there is one major step missing from this testing funciton.
% That is averaging over multiple tests from the same setup. This will need
% to be added in the full "Thrust.m" function


%% Outputs:
% ThrustCurves:
%   A table containing 0.5 seconds of thrust data for each of the cases
%   available, this data will have formatting such that there are 501
%   evenly spaced thrust data points (rows) for each test (columns). The
%   ordering of the columns will go from max to min water volume in the 2L
%   bottle and then max to min in the 1.25L bottle
%
% Time:
%   A 1D array corresponding to the times of the thrust data points in the
%   ThrustCurves table
%
% <User defined variables for statistics>

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),1);

%% List what configuration is being used
bottleSize = 2000; % [ml]
%waterSize = ;% [ml]
% Be sure that the above matches up with the file name specified below. In
% the full "Thrust" script, all of this data will be loaded automatically
% for you
testName = 'Variable Water Volume/LA_Test_W1000_B2000'; % Should be a string of the path to the data (including the data file name)

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////

%% Load data
% This should not have to be modified
fileName = testName; % again weird indexing is due to string arrays, we have to ask for all the characters in a row
data = readmatrix(fileName); % load the data
data = data(:,3)*4.448; % take only the third column and converting from lbf to N
dataTime = 0:(1/f):(1/f)*(length(data)-1);



%% Finds max thrust and crops data to 0.5 seconds

[maxThrust, ~] = max(data); % Finds max thrust index to help find the thrust start and end

% Find the start of the thrust impulse
ifFoundStart = false;
for i = find(data == maxThrust)-100:length(data) % Starts close to max thrust and moves towards max thrust
    if (data(i) >= 30) && (ifFoundStart ~= true)
        sectionStartIndex = i;
        ifFoundStart = true;
    end
end
sectionStartIndex = sectionStartIndex-8;
sectionStartTime = dataTime(sectionStartIndex);
sectionEndIndex = round( (dataTime(sectionStartIndex)+0.5)/(1/f) ); % Adds 0.5s onto the start time



% Crops the data to the 0.5 seconds of wanted thrust data
data = data(sectionStartIndex:sectionEndIndex);
dataTime = dataTime(sectionStartIndex:sectionEndIndex);

dataTime = dataTime-sectionStartTime; % Modifies the time to the 0.5 seconds
data(1) = 0;

RawData = data;



%% Data Conditioning
% Find the end of the thrust impulse
ifFoundStart = false;
for i=length(data):-1:find(data == maxThrust) % Starts at the end and moves towards max thrust
    if (data(i) >= 15) && (ifFoundStart ~= true)
        thrustEndIndex = i;
        ifFoundStart = true;
    end
end
thrustEndIndex = round( (dataTime(thrustEndIndex)+0.01)/(1/f) );

% Sets all data after the end of the impulse to zero
for i = thrustEndIndex:length(data) 
        data(i) = 0;
end

% Modifies the data based on a offset
dataOffset = cat(1,zeros(round(.03*f),1), (0:(5/(.25*f)):5)',(zeros(round(.22*f)-1,1)));
ConditionedData = data-dataOffset;

% Sets negative thrust to zero
for i=1:length(ConditionedData) 
    if ConditionedData(i) < 0
        ConditionedData(i) = 0;
    end
end


figure('Position', [40 350 500 400]); hold on; grid on; grid minor;
plot(dataTime',RawData);
plot(dataTime',ConditionedData);
% scatter(dataTime',RawData,5);
% scatter(dataTime',ConditionedData,5);
legend("Raw Data", "Conditioned Data");
title({"Static Test, Thrust over Time     Test: " + testName(32:33)});
xlabel("Time [s]");
ylabel("Thrust [N]");



%% Data Fitting
offset = -3;
coefficients = polyfit(dataTime(1:find(data == maxThrust)+offset), ConditionedData(1:find(data == maxThrust)+offset),2);
xFit = linspace(min(dataTime(1:find(data == maxThrust)+offset)), max(dataTime(1:find(data == maxThrust)+offset)));
yFit = polyval(coefficients, xFit);

coefficients2 = polyfit(dataTime(find(data == maxThrust)+offset:thrustEndIndex), ConditionedData(find(data == maxThrust)+offset:thrustEndIndex),3);
xFit2 = linspace(min(dataTime(find(data == maxThrust)+offset:thrustEndIndex)), max(dataTime(find(data == maxThrust)+offset:thrustEndIndex)));
yFit2 = polyval(coefficients2, xFit2);

% coefficients = polyfit(dataTime, ConditionedData,10);
% xFit = linspace(min(dataTime), max(dataTime));
% yFit = polyval(coefficients, xFit);

ylineEnd = zeros(1,length(data)-length(yFit)-length(yFit2)+1);
yline = cat(2,yFit,yFit2,ylineEnd);
xline = cat(2,xFit,xFit2,linspace(dataTime(1,thrustEndIndex), dataTime(1,end),length(ylineEnd)));
yline(1) = 0;

for i=1:length(yline) 
    if yline(i) < 0
        yline(i) = 0;
    end
end

plot(xline,yline,'k','LineWidth',1);



%% Sample onto the standard output array format

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////


