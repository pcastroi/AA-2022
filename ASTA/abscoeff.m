close all; clc; clear all;

%speed of sound
c= 343;
rho=1.225;
 
%reverberation times, for frequencies:
T20 = [0.613888888888889,1.320333333333333,1.646333333333333,2.838666666666666,3.539666666666667,3.130546296296296,2.026666666666667,1.122407407407408];

%the volume, approx:
V = 10356;
 
%the total area, approx:
area = 2461.83;
 
%the area of the floor:
floorarea = 989.04;
original_asphalt = [0.05 0.05 0.08 0.01 0.06 0.08 0.08 0.08]; % From some asphalt


% ODEON estimated alphas
absasph = [0.999 0.376 0.581 0.32 0.329 0.393 0.495 0.285];
absfab = [0.918 0.375 0.135 0.098 0.046 0.022 0.048 0.326];

ODEON_asp = [0.36 0.36 0.44 0.31 0.29 0.39 0.25 0.25];
 
EquiAbs = 6*log(10)*4*V./(c*T20);
 
LeftoeerAbs = EquiAbs-floorarea*ODEON_asp;

absfabsabine = LeftoeerAbs./(area-floorarea);
 
freq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]; %the corresponding frequencies

% Estimating alpha from mass law (theta = 45 degrees -  to resemble diffuse field incidence)

theta=pi/4;
mass=3; % 2 Kg

absfabmass=1./(1+((2*pi.*freq.*mass)./(2*rho*c)*cos(theta)).^2);


figure
semilogx(freq,absfabsabine,freq,absfabmass)
grid on
xlabel('frequency (Hz)')
ylabel('Absorption coefficient \alpha')
title('Estimated \alpha''s')
legend('\alpha_{fabric} Sabine''s equation','\alpha_{fabric} Mass Law',Location='best')
xticks(freq)

figure
semilogx(freq,absfab,freq,absasph)
grid on
xlabel('frequency (Hz)')
ylabel('Absorption coefficient \alpha')
title('Estimated \alpha''s')
legend('\alpha_{fabric} genetic material optimizer ODEON','\alpha_{asphalt} genetic material optimizer ODEON',Location='best')
xticks(freq)