function [N0,THM,k] = NoiseLevelDetermination(S,N,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Hildebrand Peter H., Objective Determination of the Noise Level in 
%   Doppler Espectra. 
%--------------------------------------------------------------------------
% Description:
%   This function estimate the noise level in a doppler spectrum using the
%   method described in the reference. Basically, it cut spectral components
%   above a threshold until it recognize a white spectrum by comparing the
%   statistical properties of the spectrum with the statistical properties
%   of an ideal white noise spectrum in each iteration. The similarity is
%   measured by computing two likelihood ratios(R1 & R2).
%--------------------------------------------------------------------------
% Inputs:
%   S: signal power spectrum 
%   N: maximun number of iterations
%   K: number of periodograms that where averaged when estimating the power
%      spectrum (see WelchSpectraEstimation)
%
% Outputs:
%   N0: noise level
%   THM: threshold below which all samples are considered noise 
%   k: iteration that indicates maximun likelihood

n = 1;
while(isempty(S) == 0)   
    if(n == N+1)
        disp('número máximo de iteraciones')
        break
    end    
    F = length(S);
    fn = 0:(F-1);
    A = sum(S);
    
    sigma2 = sum(fn.^2.*S)/A - (sum(fn.*S)/A)^2;
    sigmaN = (F^2-1)/12;
    R(1,n) = sigmaN/sigma2;
    
    P = mean(S);
    Q = var(S);
    R(2,n) = P^2/(Q*K);
    
    N0v(n) = P;
    
    % Threshold selection. You can uncomment the one that best suits you.
    % For fewer samples I think is best removing one spectral component at
    % a time.
    
    % Ths(n) = max(S);
    Ths(n)=1/2*(max(S)+min(S));
    
    S = S(S<Ths(n));
    n = n+1;
end

err = abs(R - ones(size(R)));
erra = err(1,:).*err(2,:);
k = (erra == min(erra));    %where the least square error is found.
N0 = N0v(k);
THM = Ths(k);
end