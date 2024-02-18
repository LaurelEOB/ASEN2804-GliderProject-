function [ThrustCurves, peakThrust, durationThrust, Time] = Thrust()
close all;
%% Thrust Summary
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, cendition their data, and fit that data into a standard
% formatting for output. Note that statistics are also requested for the
% student deliverable, but how to pass those out will be left up to the
% students as they are not needed to be passed into any later functions.
% Despite this, the first two outputs of the funciton are not permitted to
% have their form modified.

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
% <User defined variable(s) for statistics>
%

%% Define data locations
% This is hard coded!!!
fileLoc_2L = 'Variable Water Volume/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Variable Water Volume/'; % path to the data files, be sure to include a trailing slash
%fileLoc_2L = 'Static Test Stand Data/2000mL Bottle/Variable Volume/'; % path to the data files, be sure to include a trailing slash
%fileLoc_1pt25L = 'Static Test Stand Data/1250mL Bottle/'; % path to the data files, be sure to include a trailing slash

%% Read in all of the avilable data and find what data there is
testInfo_2L = getThrustTestNames(fileLoc_2L);
configs_2L = unique(testInfo_2L.waterVol);
numConfigs_2L = length(configs_2L);

testInfo_1pt25L = getThrustTestNames(fileLoc_1pt25L);
configs_1pt25L = unique(testInfo_1pt25L.waterVol);
numConfigs_1pt25L = length(configs_1pt25L);

numConfigs = numConfigs_2L + numConfigs_1pt25L;

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),numConfigs);

ThrustCurvesNames = {};

%% Loop over all of the configurations
for N = 1:numConfigs % use upper case N to distiguish that it is counting something different from the aerodynamic modeling loops
    %% Dertemine what configuration to use for this iteration in the loop
    if N <=  numConfigs_2L % determine if we should be reading 2L or 1.25L data
        bottleSize = '2000'; % [ml]
        waterSize = configs_2L(N);
        testIndexes = find(testInfo_2L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_2L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    else
        bottleSize = '1250'; % [ml]
        waterSize = configs_1pt25L(N-numConfigs_2L);
        testIndexes = find(testInfo_1pt25L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_1pt25L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    end


    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    % Notice that there is little to no guidance in place for this
    % function. This is on purpose as there are many different and equally
    % valid ways to process data (not to say that any way is valid though).
    % The lack of guidance is therefore to encourage you to think about,
    % discuss, and potentially debate as a group the best set of steps to
    % extract just the meaningful part of the thrust profile

    for j = 1:numTests
        %% Load data
            % The folloowing three lines will pull all of the files in each
            % test setup for you and give the array "data" which should be
            % conditioned. You should not need to modify any of this section of
            % code
            fileName = testNames(j, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
            data = readmatrix(fileName); % load the data
            data = data(:,3)*4.448; % take only the third column and converting from lbf to N
            dataTime = 0:(1/f):(1/f)*(length(data)-1);
            
            if min(data) <= -100
                fileName = testNames(j-1, :);
                data = readmatrix(fileName);
                data = data(:,3)*4.448;
                dataTime = 0:(1/f):(1/f)*(length(data)-1);
            end
            
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
                if (data(i) >= 20) && (ifFoundStart ~= true)
                    thrustEndIndex = i;
                    ifFoundStart = true;
                end
            end
            thrustEndIndex = round( (dataTime(thrustEndIndex)+0.01)/(1/f) );
            
            % Sets all data after the end of the impulse to zero
            for i = thrustEndIndex:length(data) 
                    data(i) = 0;
            end
            if thrustEndIndex > length(dataTime)
                thrustEndIndex = length(dataTime);
            end
            thrustEndTime = dataTime(thrustEndIndex);
            thrustDuration = thrustEndTime-0.03;
            thrustOffset = 10;
            % Modifies the data based on a offset
            dataOffset = cat(1,zeros(round(.03*f),1), (0:(thrustOffset/(thrustDuration*f)):thrustOffset)',(ones(round((0.5-thrustEndTime)*f)-1,1))*thrustOffset);
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
            title("Static Test, Thrust over Time");
            xlabel("Time [s]");
            
            dataArray(:,j) = ConditionedData(:,1);
    
    

    end


    

    %% Averaging
    averagedData = mean(dataArray');
    figure('Position', [500 350 500 400]);
    plot(dataTime,averagedData);
    title("")
    [maxThrustAvg, ~] = max(averagedData);
    % Note that averaging should accour before data fitting. Technically
    % either can be done, but the output will be much more smooth if the
    % fit is applied to an average
    %% Data Fitting
%% Data Fitting
            ifFoundStart2 = false;
            for k=length(averagedData):-1:find(averagedData == maxThrustAvg) % Starts at the end and moves towards max thrust
                if (averagedData(k) >= 15) && (ifFoundStart2 ~= true)
                    thrustEndIndexAvg = k;
                    ifFoundStart2 = true;
                end
            end
            thrustEndIndexAvg = round( (dataTime(thrustEndIndexAvg)+0.01)/(1/f) );

    offset = 3;
    coefficients = polyfit(dataTime(1:find(averagedData == maxThrustAvg)+offset), averagedData(1:find(averagedData == maxThrustAvg)+offset),2);
    xFit = linspace(min(dataTime(1:find(averagedData == maxThrustAvg)+offset)), max(dataTime(1:find(averagedData == maxThrustAvg)+offset)));
    yFit = polyval(coefficients, xFit);
    
    coefficients2 = polyfit(dataTime(find(averagedData == maxThrustAvg)+offset:thrustEndIndexAvg), averagedData(find(averagedData == maxThrustAvg)+offset:thrustEndIndexAvg),3);
    xFit2 = linspace(min(dataTime(find(averagedData == maxThrustAvg)+offset:thrustEndIndexAvg)), max(dataTime(find(averagedData == maxThrustAvg)+offset:thrustEndIndexAvg)));
    yFit2 = polyval(coefficients2, xFit2);
    
    % coefficients = polyfit(dataTime, ConditionedData,10);
    % xFit = linspace(min(dataTime), max(dataTime));
    % yFit = polyval(coefficients, xFit);
    
    ylineEnd = zeros(1,length(averagedData)-length(yFit)-length(yFit2)+1);
    yline = cat(2,yFit(1,1:end-1),yFit2(1,1:end-1),ylineEnd);
    xline = cat(2,xFit(1,1:end-1),xFit2(1,1:end-1),linspace(dataTime(1,thrustEndIndexAvg), dataTime(1,end),length(ylineEnd)));
    yline(1) = 0;
    
    for i=1:length(yline) 
        if yline(i) < 0
            yline(i) = 0;
        end
    end
    hold on;
    plot(xline,yline,'k','LineWidth',1);

    C = unique(xline,"first");

    %% Sample onto the standard output array format
    newXaxis = linspace(0,0.5,501);
    tableData = interp1(xline, yline, newXaxis);
    thrustOut = tableData;
    figure('Position', [1000 350 500 400]);
    plot(newXaxis,tableData);

    


% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = thrustOut;
    %ThrustData(:, N) = 
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
    %ThrustData{N} = 

    peakThrust(N) = max(averagedData);
    durationThrust(N) = dataTime(thrustEndIndexAvg);

end


ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames;

end
