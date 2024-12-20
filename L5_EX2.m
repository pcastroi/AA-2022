%% Lecture 5 Ar.A - Exercise 2
% Clear stuff
clear
close all
clc

%% BA - 2.1
% Q1
lx=1.1;
ly=1.7;
h=3*10^-3;
f=[100 500 1000];
fc=13/h;
c=343;
rho=1.22;
U=2*(lx+ly);

f11=c^2/(4*fc)*((1/lx)^2+(1/ly)^2);
rad_eff_res=zeros(1,3);
for i=1:length(rad_eff_res)
    if i==1
        rad_eff_res(i)=c^2/(lx*ly*fc^2);
    else
        rad_eff_res(i)=U*c/(pi^2*lx*ly*fc)*sqrt(f(i)/fc);
    end
end
rad_eff_res_db=10*log10(rad_eff_res);

% Q2
n=0.01;
rho_g=2.5*10^3;
m=rho_g*h;

R0=20*log10(2*pi*f*m/(2*rho*c));
Rr=R0-10*log10(pi*rad_eff_res.^2*fc/(2*f*n));

%% BA - 2.2
% Q2
lx=4.2;
ly=2.4;
h=22*10^-3;
kc=22;
fc=kc/h;
n=0.013;
rho_w=0.6*10^3;
m=rho_w*h;
U=2*(lx+ly);

f=100:3150;
k=2*pi*f/c;
e=4*lx*ly/U;
rad_eff_for_dif=1/2*(0.2+log(k*e));

figure
subplot(2,1,1)
plot(f,rad_eff_for_dif)
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('$\sigma_{for,diff}$','Interpreter','latex')
subplot(2,1,2)
semilogx(f,rad_eff_for_dif)
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('$\sigma_{for,diff}$','Interpreter','latex')

% Q3
R0=20*log10(2*pi*f*m./(2*rho*c));
Rf=R0+10*log10(((1-(f/fc).^2).^2)+n^2)-10*log10(2*rad_eff_for_dif);
figure
subplot(2,1,1)
plot(f,Rf)
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('$R_r$','Interpreter','latex')
subplot(2,1,2)
semilogx(f,Rf)
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('$R_r$','Interpreter','latex')

%% BA - 2.3
% Q1
w=2*pi*500;
T=0.22;
n=55/(4*w*T);

%% BA - 2.4
% Q1
h1=0.2;
lx1=3;
ly1=7;
h2=0.1;
lx2=3;
ly2=2;

mod_den_12_Q1=lx1*ly1*h2/(lx2*ly2*h1);

% Q2
mod_den_12_Q2=lx1*ly1/(lx2*ly2);

%% BA - 2.5
% Q1
f=[250 500 1000 2000];
R=[25 30 34 36];
S0=2;
e0=1/2*sqrt(S0);
S_win=8;
e=1/2*sqrt(S_win);
S_floor=30;
h_room=2.8;
T=0.5;
d=2;
L1=68;
tetha_Q1=60*pi/180;
tetha_Q2=90*pi/180;

k=2*pi*f/c;
A2=0.16*S_floor*h_room/T;
rad_eff_dif=1/2*(0.2+log(k.*e0));
rad_eff_angle_Q1=((cos(tetha_Q1)^2-0.6*2*pi./(k.*e)).^2+pi*(0.6*2*pi./(k.*e)).^2+(2*pi./(k.*e).^2).^4).^(-1/4);
rad_eff_angle_Q2=((cos(tetha_Q2)^2-0.6*2*pi./(k.*e)).^2+pi*(0.6*2*pi./(k.*e)).^2+(2*pi./(k.*e).^2).^4).^(-1/4);
L2_Q1=L1-R-10*log10(rad_eff_dif)+10*log10(rad_eff_angle_Q1)+10*log10(S_win/A2);
L2_Q2=L1-R-10*log10(rad_eff_dif)+10*log10(rad_eff_angle_Q2)+10*log10(S_win/A2);


figure
subplot(2,2,1)
plot(f,L2_Q1)
grid on
title('Q1')
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('Magnitude [dB]','Interpreter','latex')
subplot(2,2,3)
semilogx(f,L2_Q1)
grid on
title('Q1')
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('Magnitude [dB]','Interpreter','latex')
subplot(2,2,2)
plot(f,L2_Q2)
grid on
title('Q2')
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('Magnitude [dB]','Interpreter','latex')
subplot(2,2,4)
semilogx(f,L2_Q2)
grid on
title('Q2')
xlabel('Frequency [Hz]','Interpreter','latex')
ylabel('Magnitude [dB]','Interpreter','latex')
