function [z,Sz] = WeatherSignalGen(factorN,M,vm,var_v,P,CSR,var_c,SNR,fc,PRI,onlyW,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Zrnic, D. S., Simulations of weatherlike doppler spectra and signals
%--------------------------------------------------------------------------
% Description:
%   This function generate data that simulate a typical echo signal
%   from the same range cell of a weather radar using a method described in the reference. 
%   It assume that both clutter and the weather have gaussian autocorrelation and 
%   therefore gaussian spectra. Noise is considered white gaussian.
%--------------------------------------------------------------------------
%  Inputs:
%   factorN: Ratio between the actual number of samples and the number of samples that is used
%      to simulate the analog signal. This emulate the effect of windowing
%   M: Number of samples
%   vm: Mean radial velocity of the weather measurable from the radar [m/s]
%   var_v: Spectral width of the weather [m/s] 
%   P: Power of the weather 
%   CSR: Clutter-Weather power ratio [dB]
%   var_c: Spectral width of the clutter [m/s]
%   SNR: Noise-Weather power ratio [dB]
%   fc: Carrier frequency of the radar [Hz]
%   PRI: Pulse Repetition Interval(Sample Time) [s]
%   onlyW: 0 - generate a signal with clutter, weather and noise
%          1 - generate a signal with only weather 
%          2 - generate a signal with only weather aand noise
%   time: 0 - generate a signal suitable for doppler processing (algorithm: GMAP)
%         1 - generate a signal suitable for time processing (algorithm: GMAP-TD)
%
%   Outputs:
%    z: data
%    Sz: ideal spectra 

 
N = factorN*M;
theta = 2*pi*rand(1,N); 

c = 3e8;
lambda = c/fc;
vs = 0.5/PRI*lambda;
if time == 1
    v = (0:N-1)*vs/N;
else
    v = (-N/2:N/2-1)*vs/N; %fftshift
end

%------------------Weather Spectra-------------------------------------------
% This simulate the effect of sampling. If you put an antialiasing filter, remove
% the the right and the left side spectra (SwR & SwL).
vplus = v + vs; % right-sided velocities
vless = v - vs; % left-sided velocities.
Sw0 = P*N*exp(-1/(2*var_v).*(v - vm).^2)/sqrt(2*pi*var_v);
SwR = P*N*exp(-1/(2*var_v).*(vplus - vm).^2)/sqrt(2*pi*var_v);
SwL = P*N*exp(-1/(2*var_v).*(vless - vm).^2)/sqrt(2*pi*var_v);
Sw = Sw0 + SwL + SwR;

%------------------Clutter Spectra-------------------------------------------
Ac = P*N*10^(CSR/10);
Sc0 = Ac*exp(-1/(2*var_c).*(v).^2)/sqrt(2*pi*var_c);
ScR = Ac*exp(-1/(2*var_c).*(vplus).^2)/sqrt(2*pi*var_c);
ScL = Ac*exp(-1/(2*var_c).*(vless).^2)/sqrt(2*pi*var_c);
Sc = Sc0 + ScL + ScR;

%------------------Noise Spectra-------------------------------------------
An = P*N/(10^(SNR/10));
Sn = An/M*ones(1,N);

if onlyW == 1
    Sz = Sw; 
elseif onlyW == 2
    Sz = Sw + Sn; 
else
    Sz = Sw + Sc + Sn;
end
P = -Sz.*log(rand(1,N));
z = ifft(sqrt(P).*exp(1j*theta)); 
z = z(1:M);
Sz = 1/N*Sz;

end
