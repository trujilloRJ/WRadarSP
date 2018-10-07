function [P,vm,std_v,error] = PulsePairProcessing2(Swn,N0,z,R,M,PRI,lambda,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Richards, Mark A., Fundamentals of Radar Signal Processing. Chapter 5, 
%   Section 4 
%--------------------------------------------------------------------------
% Description:
%   This function estimate the first three spectral moments of a weather
%   signal using the method described in the reference. Signal spectre MUST
%   be symmetric.
%--------------------------------------------------------------------------
% Inputs:
%   Swn: signal power spectre. It MUST be symmetric (only noise and one weather peak)
%        no fftshift.
%   N0: signal noise level
%   z: sampled signal in time domain
%   R: autocorrelation matrix of the signal
%   M: number of samples
%   PRI: pulse repetition interval (sample time in a range bin)
%   lambda: wavelenght of the carrier
%   time: 0 - spectral moments will be estimated using the spectre (Swn)
%         1 - spectral moments will be estimated using sampled signal in
%         time domain (z)
%         2 - spectral moments will be estimated using autocorrelation
%         matrix (R)
%
% Outputs:
%   P: weather power estimator
%   vm: weather radial mean velocity estimator
%   std_v: weather spectral width
%   error: 0 - no errors
%          1 - error. (R0 < |R1|), although this has no teoretical meaning 
%              it could happen for some practical applications. When this 
%              happens the spectral width estimation is negative wich don't 
%              have physical meaning


N1 = length(find(Swn)); % only has menaing for GMAP.
unaso = ones(1,M);
k = 1:M;
T = PRI;
vs = 0.5/PRI*lambda;

error = 0;

if time == 1
    % Time domain
    [rzz,lags] = xcorr(z,'unbiased');
    rz = rzz(find(lags >= 0));
    Rzz = toeplitz(conj(rz),rz)*vs; % autocorrelation matrix   
    R0 = mean(diag(Rzz) - N0*ones(M,1)); 
    R1 = mean(diag(Rzz,1));
    if mean(diag(Rzz)) < abs(R1)
        error = 1;
    end    
elseif time == 2
    R0 = mean(diag(R) - N0*ones(M,1)); 
    R1 = mean(diag(R,1));
    if mean(diag(R)) < abs(R1)
        error = 1;
    end
else
    % Frequency domain
    R0 = (sum(Swn) - N0*N1)*vs/M; 
    R1 = vs/M*sum(Swn(k).*exp(1j*2*pi*(k-unaso)/M));    
    if (sum(Swn)*vs/M) < abs(R1)
        error = 1;
    end
end

if R0 > abs(R1)
    P = R0;
elseif (time == 2) || (time == 1)
    R0 = R0 + N0*vs;
    P = R0;
else
    R0 = R0 + N0*N1*vs/M;
    P = R0;     
end

vm = lambda/(4*pi*T)*angle(R1);
std_v = sqrt(lambda^2/(8*pi^2*T^2)*log(P/abs(R1)));  


end
