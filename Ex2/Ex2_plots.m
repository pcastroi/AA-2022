%% Clear stuff
% clear variables
% clear mex
% close all
% clc
%% Install subfolders
addpath meas
addpath tools

%% Load Stuff
filename='G8s1r1.wav';
[h,fs] = audioread(filename);

%% Recover IR
h_norm = h ./ max(max(abs(h))); % Normalize IRs to the overall maximum

N = 2^8;          % Window size
M = N;             % Number of DFT points
R = 2^4;           % Step size
w = blackman(N);   % Analysis window

DRdB = 80;

%% Choose the Impulse response
for i=1:width(h_norm)
h = h_norm(:,i);

[H,tSec,fHz] = stft(h,fs,w,R,M);

% Plot STFTs
a=plotSTFT(tSec,fHz,H,fs,false,DRdB);title(['Impulse Response']);

% Trunacte and select which IRs to process
% Truncate the IR (if needed) to remove most of the part that is just noise,
% keeping a short part to allow estimating the noise floor.

%% Calculate the EDC and reverberation time

% Choose an appropriate truncation time for the EDC calculation
trunctime = 0.6; % no truncation if trunctime = length(h)/fs

% Calculate the EDC
[EDC_log, t] = calcEDC(h, fs, trunctime);

% Plot the EDC
ETC_log = 10*log10(h.^2);
figure
hold on
plot(t,ETC_log) %choose truncation time by looking at this plot 
          %(when there is no decay anymore, it is only backgroud noise)

% Choose appropriate fitting points for the RT60 calculation
L1 = -5;
L2 = -25;
L3 = -35;

% Select which EDC to process
% Calculate  the reverberation time
getReverbTime(EDC_log, fs, L1, L3)
end