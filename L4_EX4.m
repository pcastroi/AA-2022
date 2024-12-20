%% Lecture 4 Ar.A - Exercise 4
% Clear stuff
clear
close all
clc

% Parameters
c=343;
rho=1.2;
rho_g=2.5*10^3;
E_g=60*10^9;
h=3*10^-3;
fc=13/h;
damp=0.02;
lx=1.1;
ly=1.7;
tetha=60*pi/180;

f=[100 500 5000];
k=2*pi.*f./c;
e=(2*lx*ly)/(lx+ly);

% Q1
rad_eff=((cos(tetha)^2-0.6*2*pi./(k.*e)).^2+pi*(0.6*2*pi./(k.*e)).^2+(2*pi./(k.*e).^2).^4).^(-1/4);
rad_eff_db=10*log10(rad_eff);

% Q2
rad_eff_dif=1/2*(0.2+log(k.*e));
rad_eff_dif_db=10*log10(rad_eff_dif);

% Q3
tau_0=(2*rho*c./(2*pi*f*rho_g*h)).^2;
tau=(tau_0.*rad_eff/cos(tetha)).*(1./(damp^2+(1-(f./fc.*sin(tetha)^2).^2).^2));
R=-10*log10(tau);