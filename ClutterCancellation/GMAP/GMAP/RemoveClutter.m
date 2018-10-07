function [Ssc,vsc,G,C,NOCLUTTER] = RemoveClutter(v,S,var_c,THM,N0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Siggia A. D., Gaussian Model Adaptive Processing(GMAP) for improved
%       ground clutter cancellation and moment calculation
%--------------------------------------------------------------------------
% Description:
%   This is the core function in the third step of GMAP. It performs
%   adaptive filtering of clutter samples by fitting a gaussian shape
%   spectrum centered in v = 0 to the clutter spectrum using the three 
%   central components.
%--------------------------------------------------------------------------
% Inputs:
%   v: velocity vector [m/s]
%   S: power spectrum [not in dB]
%   var_c: square clutter spectral width [(m/s)^2]
%   THM: threshold below which all samples are considered noise (see NoiseLevelDetermination)
%   N0: noise level
%
% Outputs:
%   Ssc: spectrum without clutter components
%   vsc: velocity samples that contain clutter components
%   G: clutter gaussian
%   C: clutter power
%   NOCLUTTER: 0 - there is clutter
%              1 - there is no clutter

NOCLUTTER = 0;
M = length(S);
res_v = v(2) - v(1);
Anum = sqrt(2*pi*var_c)*(S(end) + S(1) + S(2));
Aden = exp(-0.5*res_v.^2/var_c) + 1 + exp(-0.5*res_v.^2/var_c);
A = Anum/Aden;
C = A;

% Comparison against noise level.
NP = 5*N0;  % Noise Power equivalent to 5 samples.
if C < NP
    NOCLUTTER = 1;
end

% Clutter gaussian
G = A/(sqrt(2*pi*var_c))*exp(-0.5*v.^2/var_c);

Ssc = S;
map1 = find(fftshift(G) > THM); % indicates whre is the clutter
Ssc(map1) = 0;
map2 = find(Ssc <= THM); % indicates where is the noise
Ssc(map2) = 0;  % now Ssc only have weather
vsc = v(fftshift(G) > THM);
end