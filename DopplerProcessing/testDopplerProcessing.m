%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   This script provides an example of the Doppler Processing techniques
%   that are performed by the functions on this folder. Data generation function
%   is located at inside the data generation folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Generation Parameters
factorN = 10;
K = 10;
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
vm = 0.2*vs;
var_v = 1;

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

onlyW = 0; time = 1;

for k = 1:K;
    [zW1(k,:),Sz] =  WeatherSignalGen(factorN,M,vm,var_v,Sp,CSR,var_c,SNR,fc,PRI,onlyW,time);
end
zW = reshape(zW1,1,M*K);

Sr = WelchSpectraEstimation(zW,K,'RECTANGULAR');
Sh = WelchSpectraEstimation(zW,K,'HAMMING');
Sb = WelchSpectraEstimation(zW,K,'BLACKMAN');

figure;
plot(vN,10*log10(fftshift(Sz)));
hold on;
plot(v,10*log10(fftshift(Sr)),v,10*log10(fftshift(Sh)),v,10*log10(fftshift(Sb)));
title('Typical Weather Radar Echo Signal Spectre')
xlabel('Doppler velocity[m/s]')
ylabel('Spectral Density[dB]')
legend('Theoretical Spectre','Spectre estimated using rectangular window','Spectre estimated using hamming window','Spectre estimated using blackman window')

%% NoiseLevelDetermination test

[N0r,THMr,kr] = NoiseLevelDetermination(Sr,N,K);
[N0h,THMh,kh] = NoiseLevelDetermination(Sh,N,K);
[N0b,THMb,kb] = NoiseLevelDetermination(Sb,N,K);

figure;
plot(v,10*log10(fftshift(Sr)),v,10*log10(fftshift(Sh)),v,10*log10(fftshift(Sb)));
hold on;
plot(v,10*log10(N0r)*ones(1,M),'--',v,10*log10(N0h)*ones(1,M),'--',v,10*log10(N0b)*ones(1,M),'--');
title('Noise Level Determination')
xlabel('Doppler velocity[m/s]')
ylabel('Spectral Density[dB]')
legend('Spectre estimated using rectangular window','Spectre estimated using hamming window','Spectre estimated using blackman window','Noise Level rectangular','Noise Level hamming','Noise Level blackman')


%% PulsePairProcessing test.

onlyW = 2; time = 1;

for k = 1:K;
    [z(k,:),Sz] =  WeatherSignalGen(factorN,M,vm,var_v,Sp,CSR,var_c,SNR,fc,PRI,onlyW,time);
end
zPPP = reshape(z,1,M*K);

S = WelchSpectraEstimation(zPPP,K,'RECTANGULAR');

figure;
plot(vN,10*log10(fftshift(Sz)));
hold on;
plot(v,10*log10(fftshift(S)));
title('Only Weather and Noise Spectre')
xlabel('Doppler velocity[m/s]')
ylabel('Spectral Density[dB]')
legend('Theoretical Spectre','Spectre estimated using rectangular window')

timePPP = 0;
[P_est,vm_est,std_v_est,error] = PulsePairProcessing2(S,N0(1),0,0,M,PRI,lambda,timePPP);

ANSWER = [ Sp P_est; vm vm_est; var_v std_v_est]




