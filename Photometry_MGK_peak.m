%% Code by Munir Gunes Kutlu and Banu Kutlu for Calipari Lab
% Contact gunes.kutlu@vanderbilt.edu

clear all; close all; %Clear all variables and close all figures
BlockDir=uigetdir('Select TDT Photometry Block'); %Get folder name
cd(BlockDir); %Change directory to photometry folder
[Tank,Block,~]=fileparts(cd); %List full directory for Tank and Block
data=TDTbin2mat(BlockDir); %Use TDT2Mat to extract data.

%% SYNAPSE- Extract Relevant Data from Data file
%Create Variables for each Photometry Channel and timestamps
Ch490=data.streams.x90R.data; %GCaMP
Ch405=data.streams.x05R.data; %Isosbestic Control
Ts = ((1:numel(data.streams.x90R.data(1,:))) / data.streams.x90R.fs)'; % Get Ts for samples based on Fs
StartTime=30000; %Set the starting sample(recommend eliminating a few seconds for photoreceiver/LED rise time).
EndTime=length(Ch490)-1000; %Set the ending sample (again, eliminate some).
Fs=data.streams.x90R.fs; %Variable for Fs
Ts=Ts(StartTime:EndTime); % eliminate timestamps before starting sample and after ending.
Ch490=Ch490(StartTime:EndTime);
Ch405=Ch405(StartTime:EndTime);


%% Function to get DeltaF/F using regression fit
%Function to get DF/F for whole session
Delta490=DeltaF_no_smooth(Ch490,Ch405);
%Delta490=data.streams.x90R.data;

%% Find peaks
[peakLoc, peakMag] = peakfinder(Delta490);
peaks = peakLoc;

%% Taking TTL times and Delta490 values around TTLs

TTL_size = numel(peaks); %change this
TTL_times = zeros(TTL_size:1);
TTL_temp = zeros(TTL_size:1);
interval_pre = -2000; %milliseconds
interval_post = 2000; %milliseconds
interval_count = abs(interval_pre) + abs(interval_post);
TTL_signal = zeros(interval_count, TTL_size);
for TTL_index = 1:TTL_size %% timestamp
    ind = peaks(TTL_index, 1); % change this %returns the position best fit between T and Po0 onset
    
    position = 0;
    for interval_ind = interval_pre:interval_post
        position = position+1;
        TTL_signal(position, TTL_index) = Delta490(ind + interval_ind);
    end
    
    TTL_temp(TTL_index) = Delta490(ind); %writes the Delta490 value for the best T+Po0 position
    TTL_times(TTL_index) = ind;
end

%% Calculate z-score
baseline_pre = -2000; %milliseconds
baseline_post = 0; %milliseconds
baseline_size = abs(baseline_pre) + abs(baseline_post);
baseline = zeros(size(TTL_signal,1));
z_all = zeros(size(TTL_signal));

for interval_index = 1:TTL_size
    baseline = TTL_signal(1:baseline_size, interval_index);
    b_mean = mean(baseline); % baseline period mean
    b_stdev = std(baseline); % baseline period stdev
    for TTL_index = 1:size(TTL_signal(:, interval_index)) % Z score per bin
        z_all(TTL_index, interval_index) = (TTL_signal(TTL_index, interval_index) - b_mean)/b_stdev;
    end
end

plot(z_all)


%Save TTL_signal to CSV file
csvwrite('peaks.csv',TTL_signal);
csvwrite('peaks(Zscore).csv',z_all);


plot(TTL_signal,'DisplayName','TTL_signal')






