%% Ar.A - Exercises
% Clear stuff
clear
close all
clc

%% BA 4.1
fprintf(['BA 4.1: \n']);
% Q1
c=344;
B=6.3*10^9*(240*10^-3)^3/(12*(1-0.02^2));
fc=343^2/(2*pi)*sqrt((240*10^-3*1700)/B);
for m=1:2
    for n=1:2
        fmn(m,n)=343^2/(4*fc)*((m/3)^2+(n/4)^2);
    end
end
cs=sqrt((6.3*10^9/(2*(1+0.02)))/1700);
fs=fc*(cs/c)^2;
fprintf(['Question 1:\n   f11 = ' num2str(fmn(1,1)) ', f12 = ' num2str(fmn(1,2)) ', f21 = ' num2str(fmn(2,1)) ', f22 = ' num2str(fmn(2,2)) ', fs = ' num2str(fs) ', fc = ' num2str(fc) ' [Hz] \n']);

% Q2
f_bands=[50 63 80 100 125 160 200 250];
% 1/3 Octave -> 2^(1/6) ; Octave -> 2^(1/2)
f_low=f_bands./2^(1/6); 
f_up=f_bands.*2^(1/6);
f_diff=f_up-f_low;
% Estimated in Q1
N_modes=[0 0 1 1 0 1 1 1];
Mod_dens=N_modes./f_diff;
% Theoretical modal density (for HF)
Mod_dtheo=ones(1,length(f_bands))*(pi*3*4)/c^2*fc;
Answ=Mod_dens;
% When the f>fc, the modal density can be approximated to be
% theoretical one
Answ([end-1 end])= Mod_dtheo([end-1 end]);
fprintf(['Question 2:\n   Answer = ' num2str(Answ) ' \n']);

% Q3
% Only taking the real part for rad_eff > fc
rad_eff_mtfc=real(1./sqrt(1-(fc./f_bands)));
rad_eff_fc=sqrt(pi*2*(3+4)*f_bands/(16*c));
% The rad_eff < fc is found at the edges (0 when f > fc)
rad_eff_ltfc=2*(3+4)*c/(pi^2*3*4*fc).*sqrt(f_bands./fc);
rad_eff_ltfc([end-2:end])=[0 0 0];
% Between fc/2 and fc, rad_eff_fc is dominating over rad_eff<fc
rad_eff_rev=rad_eff_fc;
% First 2 values are below fc/2, = c^2/(3*4*fc^2)
rad_eff_rev([1:2])=[c^2/(3*4*fc^2) c^2/(3*4*fc^2)];
e=4*3*4/(2*(3+4));
rad_eff_for=0.5*(0.2+log(2*pi.*f_bands*e/343));

R0=20*log10(pi*f_bands*1700*0.24/(1.2*c));
Rr=R0-10*log10(c^2*rad_eff_rev.^2/(2*0.02*3*4*f_bands)*Answ);
Rf=R0+10*log10((1-f_bands.^2/fc^2).^2+0.02^2)-10*log10(2.*rad_eff_for);
R=-10*log10(10.^(-Rf./10)+10.^(-Rr./10));
fprintf(['Question 3:\n   R0 = ' num2str(R0) ' \n   Rr = ' num2str(Rr) ' \n   Rf = ' num2str(Rf) ' \n   R = ' num2str(R) ' \n']);

figure
plot(f_bands,R,Color='b',Marker='o')
hold on
plot(f_bands,Rf,Color='r',Marker='square',LineStyle='-.')
plot(f_bands,Rr,Color='g',Marker='+')
xlabel('f')
ylabel('R')
title('BA 4.1')
xticks(f_bands)
legend('R','R_f','R_r',Location='best')
grid on

%% BA 4.2
Ex=7.8*10^9;
Ey=10^9;
den=540;
mu=0.1;
h=22*10^-3;

fcx=343^2/(2*pi)*sqrt(12*den*(1-mu^2)/(Ex*(h)^2));
fcy=343^2/(2*pi)*sqrt(12*den*(1-mu^2)/(Ey*(h)^2));
fprintf(['BA 4.2: \n fcx = ' num2str(fcx) ' [Hz] \n fcy = ' num2str(fcy) ' [Hz] \n']);

%% BA 4.3
h_gyp=13*10^-3;
den_gyp=0.84*10^3;
h_poly=55*10^-3;
den_poly=100;
E_poly=4*10^6;
mu=0.4; 

f_dil=1/(2*pi)*sqrt(4*E_poly/(3*(1-2*mu)*h_poly*(2*den_gyp*h_gyp+den_poly*h_poly/3)));
fprintf(['BA 4.3: \n fdil = ' num2str(f_dil) ' [Hz] \n']);

%% BA 4.4
h_air=150*10^-3;
den_air=1.2;

% Q1 - 1 plate on each side
f0_q1=1/(2*pi)*sqrt(den_air*c^2/h_air*2/(h_gyp*den_gyp));
fd=c/(2*pi*h_air);
fc=34/h_gyp;
fprintf(['BA 4.4: \n']);
fprintf(['Question 1:\n   f0_q1 = ' num2str(f0_q1) ' [Hz] , fd = ' num2str(fd) ' [Hz] , fc = ' num2str(fc) ' [Hz] \n']);

% Q2 - 2 plates on each side
f0_q2=1/(2*pi)*sqrt(den_air*c^2/h_air*1/(h_gyp*den_gyp));
fprintf(['Question 2:\n   f0_q2 = ' num2str(f0_q2) ' [Hz] , fd = ' num2str(fd) ' [Hz] , fc = ' num2str(fc) ' [Hz] \n']);

% Q3 - 1 plate on 1 side, 2 plates on the other
f0_q3=1/(2*pi)*sqrt(den_air*c^2/h_air*3/(2*h_gyp*den_gyp));
fprintf(['Question 3:\n   f0_q3 = ' num2str(f0_q3) ' [Hz] , fd = ' num2str(fd) ' [Hz] , fc = ' num2str(fc) ' [Hz] \n']);

% Q4 - 1 plate on each side + mineral wool instead of air
% The stiffness of the spring with the mineral wool will be 1.4 times lower than if it's just air.
f0_q4=1/(2*pi)*sqrt(den_air*c^2/(1.4*h_air)*2/(h_gyp*den_gyp));
fprintf(['Question 4:\n   f0_q4 = ' num2str(f0_q4) ' [Hz] - The stiffness of the spring with the mineral wool will be 1.4 times lower than if it''s just air. \n']);

%% BA 4.5
e=0.2;
% Q1 - 13mm gypsum boards
fc=34/h_gyp;
fcp=2*den_gyp*h_gyp*fc/(2*den_gyp*h_gyp);
corr_q1=20*log10(e*fcp)-45;
fprintf(['BA 4.5: \n']);
fprintf(['Question 1:\n   Correction term Q1 = ' num2str(corr_q1) ' [dB] - fc and fcp are = when is a double construction with same materials \n']);

% Q2 - 150mm concrete
h_conc=h_air;
fc=19/h_conc;
corr_q2=20*log10(e*fc)-45;
fprintf(['Question 2:\n   Correction term Q2 = ' num2str(corr_q2) ' [dB] - It cannot be negative so: Correction term Q2 = 0 [dB] \n']);

%% BA 4.6
f_bands2=[100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150];
Kc=12;
h_glass=12*10^-3;
den_glass=2.5*10^3;
fc=Kc/h_glass;
S_glass=1.21*1.34;
U_glass=2*(1.21+1.34);

% Q1 - 1 glass panel
R=zeros(1,length(f_bands2));
R0=zeros(1,length(f_bands2));
Rt=zeros(1,length(f_bands2));
Rr=zeros(1,length(f_bands2));
% For f<fc
for i=1:10
    R0(1,i)=20*log10(pi*f_bands2(i)*h_glass*den_glass/(c*den_air));
    Rt(1,i)=R0(1,i)-10*log10(0.2+log(2*pi*f_bands2(i)/c*sqrt(S_glass)))+10*log10((1-(f_bands2(i)/fc)^2)^2);
end

% For f>=fc
for i=11:length(f_bands2)
    R0(1,i)=20*log10(pi*f_bands2(i)*h_glass*den_glass/(c*den_air));
    Rr(1,i)=R0(1,i)-10*log10(pi*fc/(2*10^-2*f_bands2(i)));
end
R1=Rt+Rr;
fprintf(['BA 4.6: \n']);
fprintf(['Question 1:\n   R0 = ' num2str(R0) ' \n   Rr = ' num2str(Rr) ' \n   Rt = ' num2str(Rt) ' \n   R1 = ' num2str(R1) ' \n']);

% Q2 - Double glass panel
d_gap=0.012;
f0=c/(2*pi)*sqrt(den_air*2/(d_gap*h_glass*den_glass));
fd=c/(2*pi*d_gap); % Value higher than highest frequency
R2=[27 28 29 29 30 32 34 36 37 38 38 39 39 40 43 46];

% For f<f0
for i=1:2
    R(1,i)=10*log10(4*(10^(R1(1,i)/10))); % +6 dB
end

% For f0<f<315 Hz
for i=3:6
    R(1,i)=R1(1,i)+R2(1,i)+10*log10(4)+20*log10(f_bands2(i)/fd); % +6 dB +20log(f/fd) 
end

% For f>315 Hz
for i=7:length(f_bands2)
    R(1,i)=R1(1,i)+R2(1,i)+10*log10(0.5*d_gap*U_glass/S_glass); % +10log(alpha*d*U/S)
end

fprintf(['Question 2:\n   R = ' num2str(R) ' \n']);


figure
plot(f_bands2,R,Color='b',Marker='o')
hold on
plot(f_bands2,R1,Color='r',Marker='square',LineStyle='-.')
plot(f_bands2,R2,Color='g',Marker='+')
xlabel('f')
ylabel('R')
title('BA 4.6')
xlim([63 8000])
ylim([0 100])
set(gca, 'XScale', 'log')
xticks(f_bands2)
legend('R - Double panel calculated','R_1 - Single panel calculated','R_2 - Laminated single panel measured',Location='best')
grid on