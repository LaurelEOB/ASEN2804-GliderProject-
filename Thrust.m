function [ThrustCurves, Time] = Thrust()
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
fileLoc_2L = 'Static Test Stand Data/2000mL Bottle/Variable Volume/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Static Test Stand Data/1250mL Bottle/'; % path to the data files, be sure to include a trailing slash

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

    for i = 1:numTests
    %% Load data
        % The folloowing three lines will pull all of the files in each
        % test setup for you and give the array "data" which should be
        % conditioned. You should not need to modify any of this section of
        % code
        fileName = testNames(i, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

    %% Data Conditioning



        
    %% Averaging



    end



    % Note that averaging should accour before data fitting. Technically
    % either can be done, but the output will be much more smooth if the
    % fit is applied to an average
    %% Data Fitting




    %% Sample onto the standard output array format




% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = thrustOut;
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
end
ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames;
end
