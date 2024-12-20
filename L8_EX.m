%% Lecture 8 Ar.A - Exercises
% Clear stuff
clear
close all
clc

%% BA - 20
rho=2300;
Kc=18;

% Q 2.1
f=500;
h=250*10^-3;
m=rho*h;
fc=Kc/h;
loss_f=1/sqrt(f);
rad_eff=1/sqrt(1-fc/f);
% Mass Law
R0=20*log10(2*pi*f*m/(2*1.2*343));

% f>fc -> Only resonant transmission (R=Rr)
R_21=R0-10*log10(pi*rad_eff^2*fc/(2*loss_f*f));
fprintf(['Question 2.1:\n   R = ' num2str(R0) ' - ' num2str(10*log10(pi*rad_eff^2*fc/(2*loss_f*f))) ' = ' num2str(R_21) ' dB \n']);

% Q 2.4
S=20;
S2=200;
f0=200;
R_24=R_21-10*log10(S2/S);
dL=40*log10(f/f0);
Ln=38-R_24+30*log10(f)-dL;
fprintf(['Question 2.4:\n   R = ' num2str(R_21) ' - ' num2str(10*log10(S2/S)) ' = ' num2str(R_24) ' dB \n   With Floor covering -> Ln = ' num2str(Ln) ' dB \n   Without Floor covering -> Ln = ' num2str(Ln+dL) ' dB \n']);

%% BA - 21
close all

% Q 3.1
f=[100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150];
R=[29 25 26 31 29 30 32 33 24 23 31 45 48 50 44 41];
B=450;
rhof=770;
rhoc=27.3;
hf=0.013;
hc=0.055;
h=2*hf+hc;
mf=hf*rhof;
mc=hc*rhoc;
m=2*mf+mc;
Ect=6*10^6;
G=2.8*10^6;
fdil=1/(2*pi)*sqrt(4*Ect/(hc*(2*mf+mc/3)));
cs=sqrt(G*h/m);
% cs << c = 343 | Then fcrit is decided by speed of bending waves (cbf)
cbf=sqrt(2*pi*f)*(B/(mf+1/2*mc))^(1/4);
fcrit=343^2/(2*pi)*sqrt([mf mc m]/B);
fprintf(['Question 3.1:\n   Resonance Frequency = f0 = fdil = ' num2str(fdil) ' Hz \n   Critical Frequency = fcrit = [fcritgypsum = ' num2str(fcrit(1)) ' fcritcore = ' num2str(fcrit(2)) ' fcrittotal = ' num2str(fcrit(3)) '] Hz \n']);
figure;plot(f,R);xline(fcrit,Label='fcrit');xline(fdil,Label='fdil');xlabel('Frequency, Hz');ylabel('Sound Reduction Index, R [dB]');title('Q 3.1'); grid on;

% Q 3.2
rhom=2300;
hm=0.1;
mm=rhom*hm;
fcm=19/hm; % Table 1.2
% Mass Law
R0m=20*log10(pi*f*mm/(1.2*343));
for i=1:length(f)
    if f(i)<fcm % R is calculated as forced transmission
        Rm(:,i)=R0m(i)-5+20*log10(1-(f(i)/fcm).^2);
    elseif f(i)>fcm % R is calculated as resonant transmission
        Rm(:,i)=R0m(i)-10*log10(pi*fcm/(2*sqrt(f(i))));
    end
end
figure;plot(f,Rm);xline(fcm,Label='fc');xlabel('Frequency, Hz');ylabel('Sound Reduction Index, R [dB]');title('Q 3.2'); grid on;

% Q 3.3
clear all
m1=230;
m2=10;
d=0.055;

