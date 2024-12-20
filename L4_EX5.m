%% Lecture 4 Ar.A - Exercise 5
% Clear stuff
clear
close all
clc

% Parameters
mu_g=0.24;
damp_g=0.005;
c=343;
rho=1.21;
rho_g=2.5*10^3;
E_g=60*10^9;
h=3*10^-3;
f=100;
m=rho_g*h;
S=[10 1];
k=2*pi*f/c;
e=sqrt(S); % approximation

% Q1
fc=c^2/(pi*h)*sqrt(3*rho_g*(1-mu_g^2)/E_g);
R0=20*log10(2*pi*f*m/(2*rho*c));

% Q2
Ra=R0-10*log10(0.23*R0);

rad_eff_dif=1/2*(0.2+log(k.*e));
Rb=R0-10*log10(2*rad_eff_dif(1));

Rc=R0-10*log10(2*rad_eff_dif(2)+pi*fc*(c^2/(S(2)*fc^2))^2/(2*f*damp_g));

Rd=Rc-10*log10(S(2)/S(1));