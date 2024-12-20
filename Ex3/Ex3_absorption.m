close all,
clear all
% Reflection and Absorption coefficients and surface impedance measurements
% [band,freq,realH43,imgH43]=read_pulse('lowh21.txt');
% [band,freq,realH34,imgH34]=read_pulse('lowh12.txt');

 [band,freq,realH43,imgH43]=read_pulse('H12.txt');
 [band,freq,realH34,imgH34]=read_pulse('H21.txt');

h21=realH43+i*imgH43;
h12=realH34+i*imgH34;


% Define constants:
%freq = []; 		% frequency vector (Hz)
rho = 1.21;		% density of air (kg/m^3)
c = 343;		% speed of sound in air at 23 Celsius (m/s)
s = 0.04;			% microphone spacing (m)
d = 0.07;           % dimension of the tube
Zair = rho*c;		% characteristic impedance of air (kg/m^2/s)
ko = (2*pi*freq)/c;	% wavenumber in air (m^-1)
x1 = (40e-2-3.5e-2);			% distance between the sample and the farther microphone

hcal=sqrt(h12.*h21);
ko=2.04*pi.*freq/c;
k2prime=1.94e-2.*sqrt(freq)/c/d;
k=ko-i*k2prime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H12=h21./hcal;
% Reflection coefficient
R = ( H12 - exp(-j.*k.*s) )./(exp(j.*k.*s) - H12).*exp(2.*j.*k.*x1); 
% H12 is Transfer function measured between two mics
% Absorption coefficient
alpha = 1 - abs(R).^2;
% Surface impedance
Zs = Zair*((1+R)./(1-R));
% Normalized Surface Impedance
Zs_n = ((1+R)./(1-R));


% Plots

figure(2)
subplot(211)
plot(freq,real(R),'b','LineWidth',2)
axis([0 2600 -1 1])
title('Reflection Coefficient - Real part','FontSize',12)
xlabel('Frequency (Hz)'), ylabel('Reflection Coefficient')
grid on
subplot(212)
plot(freq,imag(R),'b','LineWidth',2)
axis([0 2600 -1 1])
title('Reflection Coefficient - Imaginary part','FontSize',12)
xlabel('Frequency (Hz)'), ylabel('Reflection Coefficient')
grid on


figure(3)
subplot(2,1,1)
plot(freq,real(Zs),'b','LineWidth',2)
xlim([0 2600])
title('Surface Impedance – Real part','FontSize',12)
xlabel('Frequency (Hz)')
grid on
subplot(2,1,2)
plot(freq,imag(Zs),'b','LineWidth',2)
xlim([0 2600])
title('Surface Impedance – Imaginary part','FontSize',12)
xlabel('Frequency (Hz)')
grid on

figure(4)
subplot(2,1,1)
plot(freq,real(Zs_n),'b','LineWidth',2)
xlim([0 2600])
title('Normalized Surface Impedance – Real part','FontSize',12)
xlabel('Frequency (Hz)')
grid on
subplot(2,1,2)
plot(freq,imag(Zs_n),'b','LineWidth',2)
xlim([0 2600])
title('Normalized Surface Impedance – Imaginary part','FontSize',12)
xlabel('Frequency (Hz)')
grid on

figure(1)
plot(freq,alpha,'b','LineWidth',2), hold on

axis([0 2600 0 1])
title('Normal incidence absorption Coefficient','FontSize',12)
xlabel('Frequency (Hz)'), ylabel('Absorption Coefficient')
grid on



