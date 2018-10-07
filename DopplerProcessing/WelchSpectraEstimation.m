function [Pf] = WelchSpectraEstimation(x,K,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Welch, Peter D., The Use of Fast Fourier Transformfor the Estimation
%       of Power Spectra: A Method Based on Time Averaging Over Short, 
%       Modified Periodograms.
%--------------------------------------------------------------------------
% Description:
%   This function estimate the spectral density of a random vector x, using
%   the method described in the reference.
%--------------------------------------------------------------------------
% Inputs:
%   x: complex random vector (IQ components echo signal)
%   K: number of segments in which the signal is divided
%   window: window to select data
%       RECTANGULAR
%       HAMMING
%       BLACKMAN
%
% Outputs:
%   Pf: Estimation of the spectral density

L = length(x)/K;
if strcmp(window,'RECTANGULAR') == 1
    w = ones(L,1);  
elseif strcmp(window,'HAMMING') == 1
    w = hamming(L);
elseif strcmp(window,'BLACKMAN') == 1 
    w = blackman(L);
end

% window energy
U = (1/L)*sum(w.^2);

if (K == 1)
    X_fseg = fft(x.'.*w);  
    Pf = ((1/(L*U)*abs(X_fseg).^2))';
else
    x_seg = reshape(x,K,L);
    wm = repmat(w',K,1);
    X_fseg = fft(x_seg.*wm,[],2);
    Pf = 1/K*sum(1/(L*U)*abs(X_fseg).^2);
end

end