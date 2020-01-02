%BLOCKPATH = ;

%data = TDTbin2mat(BLOCKPATH); 

%%
%Delta=Delta490';
Delta=Ch490;
%Delta=Ch405;
%time = (1:length(data.streams.x90R.data))/data.streams.x90R.fs;
time = Ts;
%figure; plot(time, data.streams.x90R.data); hold on; 
figure; plot(time, Delta); hold on; 
%figure; plot(time, z_all); hold on; 

% StimSync epoc event
STIMSYNC = 'Po6_';
pc0_on = data.epocs.(STIMSYNC).onset;
pc0_off = data.epocs.(STIMSYNC).offset;
pc0_x = reshape(kron([pc0_on, pc0_off], [1, 1])', [], 1);
%Make a time series waveform of epoc values and plot them.

sz = length(pc0_on);
d = data.epocs.(STIMSYNC).data';
pc0_y = reshape([zeros(1, sz); d; d; zeros(1, sz)], 1, []);
hold on; plot(pc0_x, 1*(pc0_y) +420, 'g', 'LineWidth', 2);
%hold on; plot(pc0_x, 0.001*(pc0_y) , 'g', 'LineWidth', 1);