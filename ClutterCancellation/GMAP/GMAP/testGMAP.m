%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   This script provides an example of the GMAP algorithm
%   performed by the function GMAP. It dives in each step of the algorithm except 
%   the last step wich is window selection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Generation Parameters
factorN = 10;
K = 100;
M = 64;
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

onlyW = 0; time = 1;

for k = 1:K;
    [z(k,:),Sz] =  WeatherSignalGen(factorN,M,vm,var_v,Sp,CSR,var_c,SNR,fc,PRI,onlyW,time);
end
zW = reshape(z,1,M*K);

%---------------- Power spectrum estimation. Step 1 ----------------------%
w = 'BLACKMAN';
S = WelchSpectraEstimation(zW,K,w);
S_dBg = 10*log10(fftshift(S));

f1 = figure('Name','Step1_SpectraEstimation');
plot(v,S_dBg,'-o');
ylim([(min(S_dBg)-5) (max(S_dBg)+5)]);
title('Estimated Spectrum')
xlabel('Velocity[m/s]')
ylabel('Power Spectral Density[dB]')
legend('Signal Spectrum')
allText = findall(gcf, 'type', 'text');
allAxes = findall(gcf, 'type', 'axes');
allFont = [allText; allAxes];
set(allFont,'FontSize',14);
im{1} = frame2im(getframe(f1));

%---------------- Noise level determination. Step 2 ----------------------%
[N0,THM,k] = NoiseLevelDetermination(S,500,K);

f2 = figure('Name','Step2_NoiseLevelDetermination');
plot(v,S_dBg,'-o');
hold on;
plot(v,10*log10(N0)*ones(1,M),'k--');
ylim([(min(S_dBg)-5) (max(S_dBg)+5)]);
title('Spectrum and Noise Level')
xlabel('Velocity[m/s]')
ylabel('Power Spectral Density[dB]')
legend('Signal Spectrum', 'Noise Level')
allText = findall(gcf, 'type', 'text');
allAxes = findall(gcf, 'type', 'axes');
allFont = [allText; allAxes];
set(allFont,'FontSize',14);
im{2} = frame2im(getframe(f2));

%---------------- Removing clutter components. Step 3 --------------------%
[SwF, sigmaw] = FitWindow(M,w);  
sigmav2c = (var_c + (sigmaw*vs)^2);
[Ssc,vsc,Gc,C,NOCLUTTER] = RemoveClutter(v,S,sigmav2c,THM,N0);
mapc = Gc > THM;
Gc_dB = 10*log10(Gc);

f3 = figure('Name','Step3_RemovingClutterSamples');
plot(v,S_dBg,'-o');
hold on;
plot(v,10*log10(N0)*ones(1,M),'k--');
plot(v,Gc_dB,'r--');
plot(v(mapc),S_dBg(mapc),'k-o');
ylim([(min(S_dBg)-5) (max(S_dBg)+5)]);
title('Finding Clutter Samples')
xlabel('Velocity[m/s]')
ylabel('Power Spectral Density[dB]')
legend('Signal Spectrum', 'Noise Level','Clutter Gaussian','Clutter Samples')
allText = findall(gcf, 'type', 'text');
allAxes = findall(gcf, 'type', 'axes');
allFont = [allText; allAxes];
set(allFont,'FontSize',14);
im{3} = frame2im(getframe(f3));

%---------------- Moments Estimation. Step 4 -----------------------------%
p1 = 0.01;
p2 = 0.1;
[P_est,vm_est,std_v_est,k,Gw,Sw,error] = iterativePPP(Ssc,N0,PRI,lambda,mapc,p1,p2);

% This is only for graphing purposes
Ssc(Ssc == 0) = N0;
Ssc_g = fftshift(Ssc);
Ssc_g(mapc) = 0;
Sw_dBg = 10*log10(fftshift(Sw));

f4 = figure('Name','Step4_RecoveringWeatherSamples');
plot(v,10*log10(Ssc_g),'-o');
hold on;
plot(v,10*log10(Gw),'k--');
plot(v(mapc),Sw_dBg(mapc),'-o');
ylim([min(10*log10(S))-5 max(10*log10(S))+5])
title('Recovering Weather Samples')
xlabel('Velocity[m/s]')
ylabel('Power Spectral Density[dB]')
legend('Spectrum Without Clutter', 'Weather Gaussian','Weather Samples')
allText = findall(gcf, 'type', 'text');
allAxes = findall(gcf, 'type', 'axes');
allFont = [allText; allAxes];
set(allFont,'FontSize',14);
im{4} = frame2im(getframe(f4));

ANSWER = [ Sp P_est; vm vm_est; var_v std_v_est]

%% For making an animated GIF.

% filename = 'GMAP/testAnimated.gif'; % Specify the output file name
% for idx = 1:4
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     end
% end
