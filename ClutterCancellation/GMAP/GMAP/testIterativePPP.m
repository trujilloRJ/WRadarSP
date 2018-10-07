%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   This script provides an evaluation of the step 4 in GMAP algorithm
%   performed by the function iterativePPP. It tries to recover the weather
%   samples that were removed along with the clutter samples during the
%   filtering stage in the previous step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Generation Parameters
factorN = 10;
K = 100;
M = 128;
N = factorN*M;
PRI = 2e-3;
fs = 1/PRI;
fc = 5e9;
c = 3e8;
lambda = c/fc;
vs = 0.5/PRI*lambda;

% Weather Parameters
Sp = 10;
vm = 0.1*vs;
var_v = 0.8;

% Noise Parameters
SNR = 10; % dB
An = Sp*N/(10^(SNR/10));
N0 = An/M*ones(1,N); % Noise Level

% Clutter Parameters
var_c = (0.2)^2;
CSR = 30; % dB

v = (-M/2:M/2 - 1)*vs/M;
vN = (-N/2:N/2 - 1)*vs/N;

%% WeatherGen and WelchSpectraEstimation test.

onlyW = 2; time = 1;

for k = 1:K;
    [z(k,:),Sz] =  WeatherSignalGen(factorN,M,vm,var_v,Sp,CSR,var_c,SNR,fc,PRI,onlyW,time);
end
zW = reshape(z,1,M*K);

S = WelchSpectraEstimation(zW,K,'RECTANGULAR');

figure;
plot(vN,10*log10(fftshift(Sz)));
hold on;
plot(v,10*log10(fftshift(S)));
title('Only Weather and Noise Spectre')
xlabel('Doppler velocity[m/s]')
ylabel('Spectral Density[dB]')
legend('Theoretical Spectrum','Spectre estimated using rectangular window')

% Let's remove ten samples from the center, simulating clutter
% cancellation.
Sf = S;
Sf(1:11) = 0; Sf(end-8:end) = 0;
figure;
plot(v,10*log10(fftshift(Sf)));
title('Only Weather and Noise Spectre with five samples removed')
xlabel('Doppler velocity[m/s]')
ylabel('Spectral Density[dB]')


% Let's estimate the moments directly without recovering samples
timePPP = 0;
[P_est0,vm_est0,std_v_est0,error0] = PulsePairProcessing2(Sf,N0(1),0,0,M,PRI,lambda,timePPP);

% Let's now try to recover the weather samples that were removed
p1 = 0.01;
p2 = 0.1;
map = find(fftshift(Sf)==0);
[P_est,vm_est,std_v_est,k,G,Sw,error] = iterativePPP(Sf,N0(1),PRI,lambda,map,p1,p2);

figure;
plot(v,10*log10(fftshift(S)),'-.');
hold on;
plot(v,10*log10(G),'k--');
plot(v,10*log10(fftshift(Sw)),'-.');
ylim([min(10*log10(S))-5 max(10*log10(S))+5])
title('Signal Spectrum after recovering Weather Samples');
xlabel('Doppler velocity[m/s]');
ylabel('Spectral Density[dB]');
legend('Original spectrum','Weather gaussian','Spectrum with recovered weather samples')

ANSWER = [ Sp P_est0 P_est; vm vm_est0 vm_est; var_v std_v_est0 std_v_est]




