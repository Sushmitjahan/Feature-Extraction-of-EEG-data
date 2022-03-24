%% 

% plot a power spectrum of resting-state EEG data for single channel.

load MDDS1ECR.mat

erp = MDDS1EC.data;


% pick a channel to calculate power spectra
chan2plot = 'F3';

singlechandata= erp(strcmpi({MDDS1EC.chanlocs.labels},chan2plot),:);

% create a time vector that starts from 0
npnts = length(singlechandata);
time = (0:npnts-1)/MDDS1EC.srate;


% plot the time-domain signal
figure(13), clf
plot(time,singlechandata)
xlabel('Time (s)'), ylabel('Voltage (\muV)')
zoom on


% static spectral analysis
hz = linspace(0,MDDS1EC.srate/2,floor(npnts/2)+1);
ampl = 2*abs(fft(singlechandata)/npnts);
powr = ampl.^2;

figure(14), clf, hold on
plot(hz,ampl(1:length(hz)),'k','linew',2)
plot(hz,powr(1:length(hz)),'r','linew',2)

xlabel('Frequency (Hz)')
ylabel('Amplitude or power')
legend({'Amplitude';'Power'})

%% 
% All channel power spectra
%  


% convert to double-precision
Allchandata = double(MDDS1EC.data);


% the second dimension (time)
chanpowr = ( 2*abs( fft(Allchandata,[],2) )/MDDS1EC.pnts );

% then average over trials
%chanpowr = mean(chanpowr,3);

% vector of frequencies
hz = linspace(0,MDDS1EC.srate/2,floor(MDDS1EC.pnts/2)+1).^2;


% do some plotting
% plot power spectrum of all channels
figure(15), clf
plot(hz,chanpowr(:,1:length(hz)),'linew',2)
xlabel('Frequency (Hz)'), ylabel('Power (\muV)')

%set(gca,'xlim',[0 30],'ylim',[0 100])

%% extracting alpha power

% boundaries in hz
alphabounds = [ 8 12 ];

% convert to indices
freqidx = dsearchn(hz',alphabounds');


% extract average power
alphapower = mean(chanpowr(:,freqidx(1):freqidx(2)),2);


% plot the time-domain signal
figure(13), clf
plot(alphapower)
xlabel('Time (s)'), ylabel('Voltage (\muV)')
zoom on




% and plot
figure(16), clf
topoplotIndie(alphapower,MDDS1EC.chanlocs,'numcontour',0);
set(gca,'clim',[0 6])
colormap hot



%% 
% 
%   VIDEO: Welch's method on resting-state EEG data
% 
% 


% load data
load EEGrestingState.mat
N = length(eegdata);

% time vector
timevec = (0:N-1)/srate;

% plot the data
figure(26), clf
plot(timevec,eegdata,'k')
xlabel('Time (seconds)'), ylabel('Voltage (\muV)')

%% one big FFT (not Welch's method)

% "static" FFT over entire period, for comparison with Welch
eegpow = abs( fft(eegdata)/N ).^2;
hz = linspace(0,srate/2,floor(N/2)+1);

%% "manual" Welch's method

% window length in seconds*srate
winlength = 1*srate;

% number of points of overlap
nOverlap = round(srate/2);

% window onset times
winonsets = 1:nOverlap:N-winlength;

% note: different-length signal needs a different-length Hz vector
hzW = linspace(0,srate/2,floor(winlength/2)+1);

% Hann window
hannw = .5 - cos(2*pi*linspace(0,1,winlength))./2;

% initialize the power matrix (windows x frequencies)
eegpowW = zeros(1,length(hzW));

% loop over frequencies
for wi=1:length(winonsets)
    
    % get a chunk of data from this time window
    datachunk = eegdata(winonsets(wi):winonsets(wi)+winlength-1);
    
    % apply Hann taper to data
    datachunk = datachunk .* hannw;
    
    % compute its power
    tmppow = abs(fft(datachunk)/winlength).^2;
    
    % enter into matrix
    eegpowW = eegpowW  + tmppow(1:length(hzW));
end

% divide by N
eegpowW = eegpowW / length(winonsets);


%%% plotting
figure(27), clf, subplot(211), hold on

plot(hz,eegpow(1:length(hz)),'k','linew',2)
plot(hzW,eegpowW/10,'r','linew',2)
set(gca,'xlim',[0 40])
xlabel('Frequency (Hz)')
legend({'"Static FFT';'Welch''s method'})
title('Using FFT and Welch''s')

%% MATLAB pwelch

subplot(212)

% create Hann window
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;

% number of FFT points (frequency resolution)
nfft = srate*100;

pwelch(eegdata,hannw,round(winsize/4),nfft,srate);
set(gca,'xlim',[0 40])


%% 
% 
%   VIDEO: Welch's method on v1 laminar data
% 
% 

clear
load v1_laminar.mat
csd = double(csd);

% specify a channel for the analyses
chan2use = 7;


% create Hann window
hannw = .5 - cos(2*pi*linspace(0,1,size(csd,2)))./2;

% Welch's method using MATLAB pwelch
[pxx,hz] = pwelch(squeeze(csd(chan2use,:,:)),hannw,round(size(csd,2)/10),1000,srate);

figure(28), clf
subplot(211)
plot(timevec,mean(csd(chan2use,:,:),3),'linew',2)
set(gca,'xlim',timevec([1 end]))
xlabel('Time (s)'), ylabel('Voltage (\muV)')

subplot(212)
plot(hz,mean(pxx,2),'linew',2)
set(gca,'xlim',[0 140])
xlabel('Frequency (Hz)')
ylabel('Power (\muV^2)')


%% done.
