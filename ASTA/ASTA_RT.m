%% Ar.A - ASTA RT
clear variables
clear mex
close all
clc

addpath('G8-ASTA\')
addpath('tools\')

f_bands = [63 125 250 500 1000 2000 4000 8000];
colors = ['#96ca2d';'#1c3ffd';'#d90000'];
colorsrgb = [150 202 45;28 63 253;217 0 0];
colors_alpha = [150 202 45 255*0.35;28 63 253 255*0.35;217 0 0 255*0.35]./255;

rpos=[11.5,6.2;12.8,13.2;6.8,15.4;4.8,20;5.9,25.7;11.4,29.3;13.5,29.2;17.6,29.6;8.5,33.5;4.3,26.5];
spos=[6,6.5;11.8,20;17,33.5];
% Change of mirrored coordinates
rpos(:,1)=24.5-rpos(:,1);
spos(:,1)=24.5-spos(:,1);

% (:,:,S1/S2/S3)
% r1.r2...r10
EDT = [1.70 1.73 1.59 2.55 2.82 2.25 1.57 0.64;0.34 1.66 1.48 3.19 3.47 3.27 1.19 0.39;0.16 1.85 1.77 2.95 3.55 2.66 1.95 0.57;0.67 1.95 1.26 2.78 3.44 2.76 1.92 0.57;0.23 1.66 1.69 3.20 3.51 2.87 1.80 0.24;1.41 1.05 1.39 2.69 3.43 2.66 1.49 1.00;0.76 1.80 0.80 2.35 2.95 2.61 1.48 0.76;0.60 1.00 1.21 2.29 2.88 2.71 1.63 0.93;0.21 1.47 1.55 2.68 3.12 2.68 1.81 0.55;5.42 0.81 1.52 2.12 2.76 2.65 1.49 0.78];
EDT(:,:,2)=[0.37 0.96 0.90 2.55 3.52 2.70 1.16 0.58;0.34 1.05 1.35 2.79 3.12 2.22 0.66 0.54;NaN 1.14 1.12 2.77 2.88 2.05 1.33 0.57;0.37 1.64 1.69 3.01 3.18 2.82 1.94 0.98;0.72 1.34 1.07 3.01 3.24 2.37 0.70 0.49;NaN 1.29 1.43 2.59 3.07 2.34 1.39 0.95;0.31 0.78 1.41 2.97 3.44 2.75 1.11 0.51;0.81 1.12 1.23 2.20 2.50 1.88 1.32 0.72;0.66 0.93 1.57 2.86 3.14 2.89 0.77 0.41;NaN 1.11 1.48 2.42 3.06 2.70 1.81 1.54];
EDT(:,:,3)=[NaN 1.18 1.00 2.31 2.82 2.18 0.95 0.80;0.10 1.46 1.22 2.25 2.68 2.19 1.32 0.59;NaN 1.41 1.66 2.63 3.46 2.74 1.69 1.10;1.08 1.67 1.53 3.05 3.26 2.24 1.51 0.36;0.98 0.48 0.78 2.28 2.49 2.18 1.19 0.42;0.09 1.25 1.32 1.96 2.40 2.18 1.57 0.89;0.07 0.78 1.52 3.20 3.45 2.77 0.70 0.52;0.10 1.31 1.82 2.62 3.19 2.53 1.38 1.08;0.69 2.16 1.18 2.34 2.69 0.64 0.30 0.24;0.15 1.30 1.69 2.91 3.29 2.92 1.29 0.36];

T20 = [2.26 1.31 1.50 2.77 3.28 3.05 1.98 1.11;NaN 1.42 1.64 2.90 3.22 2.78 1.96 1.18;0.13 1.56 1.73 3.01 3.53 3.15 2.02 1.04;NaN 1.32 1.78 2.86 3.58 2.91 2.05 1.12;NaN 1.11 1.65 2.97 3.56 3.32 1.99 1.10;NaN 1.36 1.65 2.95 3.60 3.07 1.98 1.15;NaN 1.43 1.53 3.03 3.94 3.15 1.98 0.95;NaN 1.28 1.63 2.63 3.63 NaN NaN 1.29;NaN 1.33 1.78 3.16 3.76 3.20 2.05 1.20;NaN 1.22 1.58 2.94 3.38 NaN NaN NaN];
T20(:,:,2) = [NaN 1.52 1.81 3.09 3.85 3.03 2.10 1.21;0.14 1.42 1.81 2.74 3.50 2.94 2.02 1.09;NaN 1.29 1.62 3.00 3.24 2.84 2.03 1.22;NaN 1.37 1.75 2.64 3.24 3.28 2.16 1.20;0.48 1.16 1.63 2.88 3.22 3.12 2.06 1.16;NaN 1.23 1.66 2.75 3.44 3.01 2.19 1.21;NaN 1.32 1.70 2.88 6.06 5.14 1.95 0.95;NaN 1.33 1.57 2.52 3.19 2.75 NaN 1.24;1.05 1.46 1.81 2.81 3.22 2.86 1.92 1.30;NaN 1.26 1.63 2.79 3.86 NaN NaN NaN];
T20(:,:,3) = [NaN 1.32 1.73 2.82 3.31 2.90 2.04 1.12;0.06 1.25 1.57 2.60 3.50 3.04 1.92 1.00;NaN 1.35 1.82 3.35 3.86 3.09 2.12 1.17;NaN 1.08 1.61 2.83 3.20 3.17 2.04 1.05;NaN 1.45 1.50 2.65 3.48 3.22 2.13 1.20;NaN 1.29 1.46 2.71 3.20 3.06 1.88 0.97;NaN 1.26 1.76 2.86 3.32 3.20 2.08 1.19;0.13 1.42 1.54 2.58 3.34 3.17 2.13 1.13;NaN 0.96 1.39 2.77 3.21 2.95 1.89 0.72;0.08 1.53 1.55 2.67 3.47 3.14 2.02 1.10];

D50 = [0.92 0.89 0.68 0.66 0.62 0.67 0.73 0.86;0.49 0.73 0.86 0.79 0.74 0.70 0.87 0.98;0.95 0.57 0.31 0.51 0.49 0.55 0.78 0.90;0.73 0.74 0.41 0.44 0.48 0.49 0.69 0.92;0.70 0.58 0.41 0.50 0.54 0.46 0.69 0.94;0.15 0.61 0.40 0.16 0.19 0.24 0.44 0.63;0.17 0.57 0.40 0.50 0.58 0.50 0.55 0.77;0.06 0.47 0.50 0.27 0.32 0.38 0.44 0.69;0.98 0.60 0.63 0.44 0.49 0.53 0.76 0.88;0.17 0.51 0.42 0.33 0.30 0.38 0.34 0.39];
D50(:,:,2) = [0.75 0.64 0.54 0.36 0.25 0.30 0.36 0.50;0.78 0.72 0.58 0.48 0.59 0.65 0.80 0.93;NaN 0.55 0.50 0.31 0.39 0.39 0.52 0.69;0.23 0.77 0.62 0.33 0.39 0.45 0.68 0.89;0.08 0.62 0.54 0.34 0.44 0.58 0.75 0.81;NaN 0.62 0.47 0.19 0.19 0.28 0.35 0.59;0.21 0.63 0.63 0.66 0.59 0.66 0.78 0.94;0.45 0.60 0.50 0.24 0.23 0.26 0.34 0.61;0.35 0.75 0.56 0.51 0.58 0.64 0.85 0.93;NaN 0.64 0.46 0.24 0.24 0.37 0.54 0.71];
D50(:,:,3) = [NaN 0.27 0.25 0.20 0.12 0.16 0.20 0.39;1.00 0.52 0.23 0.42 0.46 0.41 0.41 0.59;NaN 0.30 0.28 0.15 0.16 0.13 0.20 0.46;0.78 0.63 0.48 0.59 0.55 0.58 0.71 0.91;0.45 0.55 0.68 0.41 0.46 0.47 0.60 0.72;1.00 0.53 0.26 0.36 0.46 0.65 0.80 0.87;1.00 0.80 0.56 0.53 0.53 0.67 0.83 0.89;1.00 0.64 0.35 0.20 0.43 0.58 0.78 0.91;0.46 0.60 0.81 0.82 0.78 0.90 0.96 0.98;0.92 0.85 0.54 0.48 0.68 0.74 0.90 0.96];

% SPL BN estimation, 3 positions, fband(2:7)
Vent_out = [79.36 75.65 74.72;70.8 64.67 67.18;71.89 67.26 69.26;69.12 63.76 63.76;68.33 62.47 61.47;64.93 58.56 57.2];
Vent_in = [75.15 75.68 74.55;68.41 68.36 66.43;70.36 71.52 69.61;69.38 71.03 68.78;68.17 70.13 68.09;67.42 70.58 69.01];

% BN estimation (10 positions)
BN = [56.7 54.8 55.4 55.4 53.1 57.6 58.1 65.3 60.1 58.7;54.1 51.8 54 54.9 53.6 55.2 59.5 63.3 56.3 58.9;45.6 46.1 45.1 45.3 47.5 45.1 50.5 54.9 50.1 50.5;51.2 51.8 52.4 51.7 51.1 53 54 55.3 54 53.1;50.8 49.9 50.8 49.8 51.1 49.8 50.9 54 50.4 51;50 48.5 49.6 49 50.2 48.6 50.2 52.8 50.5 49.6;44.1 41.3 43.7 42.7 44.5 42.4 44.6 49.1 45.3 45.7;30.6 26.7 30 26.7 31.2 26.4 32.7 38.1 33 32.9];

% Odeon estimates NaN differently than just ignoring NaN values
EDT_avg=mean([nanmean(EDT(:,:,1));nanmean(EDT(:,:,2));nanmean(EDT(:,:,3))]);
T20_avg=mean([nanmean(T20(:,:,1));nanmean(T20(:,:,2));nanmean(T20(:,:,3))]);
D50_avg=mean([nanmean(D50(:,:,1));nanmean(D50(:,:,2));nanmean(D50(:,:,3))]);

EDT_std=mean([nanstd(EDT(:,:,1));nanstd(EDT(:,:,2));nanstd(EDT(:,:,3))]);
T20_std=mean([nanstd(T20(:,:,1));nanstd(T20(:,:,2));nanstd(T20(:,:,3))]);
D50_std=mean([nanstd(D50(:,:,1));nanstd(D50(:,:,2));nanstd(D50(:,:,3))]);

% Plots
figure
errorbar(f_bands,EDT_avg,EDT_std)
xticks(f_bands)
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title(['EDT = ' num2str(mean(EDT_avg(:))) ' s'])
xlim([50 10000])
xlabel('Frequency (Hz)')
ylabel('Time (s)')
grid on

% Plots
figure
for i=1:3
    for j=1:10
        hold on
        plot(f_bands,(T20(j,:,i)),'Color',colors_alpha(i,:));
    end
end
e=errorbar(f_bands,T20_avg,T20_std);
e.LineWidth=2;
e.Color='k';
xticks(f_bands)
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title(['Average T_{20} = ' num2str(mean(T20_avg(:))) ' s'])
xlim([50 10000])
xlabel('Frequency (Hz)')
ylabel('Time (s)')
grid on

figure
for i = 1:3
    for j=1:10
    hold on
    plot3(f_bands,i*ones(length(f_bands)),T20(j,:,i),'Color',colors(i,:))
    end
end
view([45 45])
xticks(f_bands)
yticks ([1 2 3])
yticklabels({'Pos 1','Pos 2','Pos 3'})
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title('Measured T_{20} with separated source positions')
xlim([50 10000])
xlabel('Frequency (Hz)')
ylabel('Source Positions')
zlabel('Time (s)')
grid on

figure
for i = 1:10
    for j=1:3
    hold on
    plot3(f_bands,i*ones(length(f_bands)),T20(i,:,j),'Color',colors(j,:))
    end
end
view([45 45])
xticks(f_bands)
yticks ([1 2 3 4 5 6 7 8 9 10])
yticklabels({'Pos 1','Pos 2','Pos 3','Pos 4','Pos 5','Pos 6','Pos 7','Pos 8','Pos 9','Pos 10'})
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title('Measured T_{20} with separated receiver positions')
xlim([50 10000])
xlabel('Frequency (Hz)')
ylabel('Receiver Positions')
zlabel('Time (s)')
grid on

figure
str = string(1:10);
for i=1:length(str)
    str(i)='R' + str(i);
end
scatter(rpos(:,1),rpos(:,2),[],"black",'filled','o')
hold on
% textscatter(rpos-1,str,"TextDensityPercentage",100)
for i=1:3
scatter(spos(i,1),spos(i,2),[],colorsrgb(i,:)./255,'filled','^')
end
% textscatter(spos+1,["S1";"S2";"S3"],"TextDensityPercentage",100)
xline(0,LineStyle="--",Color='red')
xline(24.5,LineStyle="--",Color='red')
yline(0,LineStyle="--",Color='red')
yline(40.5,LineStyle="--",Color='red')
xlim([0 24.5])
ylim([0 40.5])
xlabel('Length (m)')
ylabel('Length (m)')
title('Source and Receiver positions in ASTA')
grid on

figure
for i=1:size(BN,2)
plot(f_bands,BN(:,i),'Color',colors_alpha(randi([1 3]),:))
hold on
end
for i=1:length(f_bands)
BN_avg(i)=mean(BN(i,:));
BN_std(i)=std(BN(i,:));
end
e=errorbar(f_bands,BN_avg,BN_std);
e.LineWidth=2;
e.Color='k';
xlim([93.75 6000])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xticks(f_bands)
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
legend('Measurement 1','Measurement 2','Measurement 3','Average',Location='best')
% SPL in DBA!!! DOOOOO THIS
title(['Average Background Noise SPL = ' num2str(mean(BN(:))) ' dB'])
grid on

%% BG NOISE
%
%
[x1,fs] = readAudio('Aud_BN_inhole_moregain.wav');
x2 = readAudio('Aud_BN_outhole_moregain.wav');

% Select linear or exponential frequency increase
isExp1 = false;
isExp2 = true;

N = length(x1); % Number of samples
Ws = 100000; % Window size

% Rectangular window
w = ones(1,Ws);
% Step size
R = round(N/2^10); % Hop Length depends of the N samples for each frame
M = round(N/2^10); % Frame Size depends of the N samples for each frame

% 
% 
% Calculate the magnitude spectra

f0 = 10; % Start frequency
f1 = fs/2; % Stop frequency

t1=linspace(0,length(x1)/fs,length(x1));
t2=linspace(0,length(x2)/fs,length(x2));

Y1 = fft(x1(:,1),length(x1));
YdB1 = 20 * log10(abs(Y1)); % Magnitude spectrum in dB
Y2 = fft(x2(:,1),length(x2));
YdB2 = 20 * log10(abs(Y2)); % Magnitude spectrum in dB

fHz1 = linspace(f0,f1,length(YdB1)); 
fHz2 = linspace(f0,f1,length(YdB2)); 

% Calculate the STFTs
[X1,T1,F1] = stft(x1(:,1),fs,w,R,M);
[X2,T2,F2] = stft(x2(:,1),fs,w,R,M);

% 
%
figure
d = tiledlayout(3,2,'Padding','Compact');
% Plot the time-domain signals
nexttile
plot(t1,x1(:,1));
grid on;
xlabel('Time (s)')
ylabel('Amplitude')
title('Background Noise @ in-hole')
xlim([t1(1) t1(end)])

nexttile
plot(t2,x2(:,1));
grid on;
xlabel('Time (s)')
ylabel('Amplitude')
title('Background Noise @ out-hole')
xlim([t2(1) t2(end)])

% Plot the magnitude spectra
nexttile
for i=1:3
plot(f_bands(2:7),Vent_in(:,i),'Color',colors_alpha(i,:))
hold on
end
for i=1:6
VNi_avg(i)=mean(Vent_in(i,:));
VNi_std(i)=std(Vent_in(i,:));
end
e=errorbar(f_bands(2:7),VNi_avg,VNi_std);
e.LineWidth=2;
e.Color='k';
xlim([93.75 6000])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xticks(f_bands(2:7))
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title(['Average SPL = ' num2str(mean(Vent_in(:))) ' dB'])
grid on

nexttile
for i=1:3
plot(f_bands(2:7),Vent_out(:,i),'Color',colors_alpha(i,:))
hold on
end
for i=1:6
VNo_avg(i)=mean(Vent_out(i,:));
VNo_std(i)=std(Vent_out(i,:));
end
e=errorbar(f_bands(2:7),VNo_avg,VNo_std);
e.LineWidth=2;
e.Color='k';
xlim([93.75 6000])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xticks(f_bands(2:7))
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
title(['Average SPL = ' num2str(mean(Vent_out(:))) ' dB'])
grid on

% Plot the STFTs
DRdB = 50;

% Spectrogram in dB
XdB1 = 20 * log10(abs(X1));
XdB2 = 20 * log10(abs(X2));

% Select dynamic range
range1 = [max(XdB1(:)) - DRdB max(XdB1(:))];
range2 = [max(XdB2(:)) - DRdB max(XdB2(:))];

% Select colormap
colormap(colormapVoicebox);

nexttile
imagesc(T1(:)',F1(:),XdB1,range1);
axis xy;
ylim([0 fs/2])
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('STFT')

nexttile
imagesc(T2(:)',F2(:),XdB2,range2);
axis xy;
ylim([0 fs/2])
colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('STFT')

figure
abscoeff
