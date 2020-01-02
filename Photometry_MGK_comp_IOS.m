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
StartTime=1000; %Set the starting sample(recommend eliminating a few seconds for photoreceiver/LED rise time).
EndTime=length(Ch490)-500; %Set the ending sample (again, eliminate some).
Fs=data.streams.x90R.fs; %Variable for Fs
Ts=Ts(StartTime:EndTime); % eliminate timestamps before starting sample and after ending.
Ch490=Ch490(StartTime:EndTime);
Ch405=Ch405(StartTime:EndTime);


%% Function to get DeltaF/F using regression fit
%Function to get DF/F for whole session
Delta490=DeltaF(Ch490,Ch405);
%Delta490=data.streams.x90R.data;

%% Generate correct and missed trials & tone and ITI LNPs
lnp_ch = data.epocs.Po0_.onset;
tone_ch = data.epocs.Po6_.onset;
lnp_ch_size = numel(lnp_ch);
correct_trial_onset = zeros(lnp_ch_size,1);
missed_trial_onset = zeros(lnp_ch_size,1);
correct_lnp_onset = zeros(lnp_ch_size,1);
correct_lnp_timestamps = [];
cto_index = 1;
mto_index = 1;
for tone_index = 1:numel(tone_ch) %% timestamp
    tone_timestamp = tone_ch(tone_index, 1);
    % tone_timestamp_interval = [tone_timestamp tone_timestamp+10];
    I = find(lnp_ch>=tone_timestamp & lnp_ch<=(tone_timestamp+10));
    if isempty(I)
        % add corresponding tone channel value
        missed_trial_onset(mto_index) = tone_timestamp;
        
        mto_index = mto_index + 1;
    else
        % add corresponding tone channel value
        correct_trial_onset(cto_index) = tone_timestamp;
        % add corresponding ttl channel value
        correct_lnp_timestamp = lnp_ch(I);
        correct_lnp_onset(cto_index) = correct_lnp_timestamp(1);
        % collect all correct_lnp_timestamps to later diff against ttl_ch to generate iti_lnp_onset
        correct_lnp_timestamps = [correct_lnp_timestamps;correct_lnp_timestamp];
        unique(correct_lnp_timestamps(:,1));
        
        cto_index = cto_index + 1;
    end
end
        
I = find(correct_trial_onset > 0);
correct_trial_onset = correct_trial_onset(I);
I = find(missed_trial_onset > 0);
missed_trial_onset = missed_trial_onset(I);
I = find(correct_lnp_onset > 0);
correct_lnp_onset = correct_lnp_onset(I);
% Diff correct_lnp_timestamps against ttl_ch to generate iti lnp
iti_lnp_onset = setdiff(lnp_ch, correct_lnp_timestamps);



%% Generate tone and ITI headentries
he_ch = data.epocs.Po4_.onset;
%tone_ch = data.epocs.Po6_.onset;
c_tone_ch = correct_trial_onset;
he_ch_size = numel(he_ch);
correct_he_onset = zeros(he_ch_size,1);
correct_he_timestamps = [];
cho_index = 1;
for tone_index = 1:numel(c_tone_ch) %% timestamp
    tone_timestamp = c_tone_ch(tone_index, 1);
   
     I = find(he_ch>=tone_timestamp & he_ch<=(tone_timestamp+30));
    
    if ~isempty(I)
        % add corresponding he channel value
        correct_he_timestamp = he_ch(I);
        correct_he_onset(cho_index) = correct_he_timestamp(1);
        % collect all correct_he_timestamps to later diff against ttl_ch to generate iti_he_onset
        correct_he_timestamps = [correct_he_timestamps;correct_he_timestamp];
        unique(correct_he_timestamps(:,1));
        
        cho_index = cho_index + 1;
    end
end
        
I = find(correct_he_onset > 0);
correct_he_onset = correct_he_onset(I);
% Diff correct_he_timestamps against he_ch to generate iti he
iti_he_onset = setdiff(he_ch, correct_he_timestamps);


%% Taking TTL times and Delta490 values around TTLs
%correct_trial_onset
%missed_trial_onset
%correct_lnp_onset
%iti_lnp_onset
%correct_om_onset
%missed_om_onset
%correct_he_onset
%iti_he_onset
%om_he_onset
%data.epocs.Po0_.onset - LNP
%data.epocs.Po2_.onset - RNP
%data.epocs.Po4_.onset - Headentry
%data.epocs.Po6_.onset - Tone

TTL_type = data.epocs.Po4_.onset; %change this
TTL_size = numel(TTL_type); 
TTL_times = zeros(TTL_size:1);
TTL_temp = zeros(TTL_size:1);
interval_pre = -2000; %milliseconds
interval_post = 18000; %milliseconds
interval_count = abs(interval_pre) + abs(interval_post);
TTL_signal = zeros(interval_count, TTL_size);
for TTL_index = 1:TTL_size %% timestamp
    [c, ind] = min(abs(Ts-TTL_type(TTL_index, 1))); %returns the position best fit between T and Po0 onset
    
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
csvwrite('all_he_onset (low smooth).csv',TTL_signal);
csvwrite('all_he_onset (low smooth) (Z Score).csv',z_all);
plot(TTL_signal,'DisplayName','TTL_signal')


%%Append TTL_signal to an existing CSV file
%directory = 'C:\Users\kutlumg\Dropbox\Public\Calipari Lab\Seeking vs. Avoidance Project\Fiber Photometry\D1 FP Data (Sucrose Test)\';
%filename  = 'TTL_signal.csv';
%fileDest  = fullfile(directory,filename);
%data = csvread(fileDest);
%newdata = [data TTL_signal];
%csvwrite(fileDest, newdata);v
%%close all


%calculate latency to nosepoke for correct trials
%latency = correct_active_NP_timestamps - correct_trial_onset
%csvwrite('latency to NP.csv',latency);
