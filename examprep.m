%% Exam prep
c=344;
% Sabine's equation - ceil diff abs than other walls
lx=10;
ly=5;
lz=3;
alpha_ceil=0.2;
alpha_walls=0.1;
Trev=55.3*(lx*ly*lz)/(lx*ly*alpha_ceil+alpha_walls*(lx*ly+2*lx*lz+2*ly*lz));

% Random incidence absorption || Normal incidence absorption
% ζ - zeta = specific surface impedance
Z = 500-1i*1000;
zeta=Z/(1.2*343);
r=real(zeta);
x=imag(zeta);
alpha_rand=8*r/(r^2+x^2)*(1-r/(r^2+x^2)*log((r+1)^2+x^2)+(r^2-x^2)/(x*(r^2+x^2))*atan(x/(r+1)));
alpha_norm=4*r/(abs(zeta)^2+2*r+1);


% Clarity (C80) and Definition (D50) exponential decays as a function of Trev
Trev=2;
C_exp=10*log10(exp(1.104/Trev)-1);
D_exp=1-exp(-0.69/Trev);

% Wavelength λ [m] - of a f=100 Hz tone
f=100;
lambda=c/f;

