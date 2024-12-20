%% Ar.A - Ex 5
clear all
close all
clc

data=importdata('Sound_insulation_data.txt');
f=data.data(:,1);
Ln=data.data(:,2);
R=data.data(:,3);

figure
plot(f,Ln)
title('Impact Sound Level L_n')
xlabel('Frequency [Hz]')
ylabel('Impact Sound Level [dB]')
grid on
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
xticks(f)
xticklabels(f)
xlim([40 6000])

figure
plot(f,R)
title('Transmission Loss R')
xlabel('Frequency [Hz]')
ylabel('Sound Transmission Loss R [dB]')
grid on
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
xticks(f)
xticklabels(f)
xlim([40 6000])
