%% Lecture 1 Ar.A - Exercise 4
%% Clear stuff
clear
close all
clc

fc=[125 250 500 1000 2000 4000];
SSI=[1.08+i*6.48 1.06+i*3.08 1.15+i*1.47 1.54+i*0.53 1.68+i*0.26 1.69+i*0.01];
r=real(SSI);
x=imag(SSI);

alpha_0=zeros(1,length(fc));
alpha_rand=zeros(1,length(fc));

for i=1:length(fc)
   alpha_0(i)=4*real(SSI(i))/abs((abs(SSI(i)))^2+2*real(SSI(i))+1);
   alpha_rand(i)=8*r(i)/(r(i)^2+x(i)^2)*(1-r(i)/(r(i)^2+x(i)^2)*log((r(i)+1)^2+x(i)^2)+(r(i)^2-x(i)^2)/(x(i)*(r(i)^2+x(i)^2))*atan(x(i)/(r(i)+1)));   
end