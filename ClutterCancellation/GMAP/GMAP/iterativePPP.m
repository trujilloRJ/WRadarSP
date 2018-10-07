function [P,vm,std_v,k,G,S,error] = iterativePPP(Swn,N0,PRI,lambda,map,p1,p2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Siggia A. D., Gaussian Model Adaptive Processing(GMAP) for improved
%       ground clutter cancellation and moment calculation
%--------------------------------------------------------------------------
% Description:
%   This is the core function in the fourth step of GMAP. It tries to recover
%   the samples of weather that were removed during the filtering stage.
%   Therefore, it estimate the first three spectral moments and then use
%   them to build a gaussian shape spectrum that represents the weather spectrum. It does
%   this task until the estimation converges.
%--------------------------------------------------------------------------
% Inputs:
%   Swn: signal power spectrum without clutter components
%   N0: noise level [not in dB]
%   PRI: pulse repetition interval [s]
%   lambda: carrier wavelength 
%   map: mask that indicate where it should recover the weather samples
%   p1, p2: coefficients that indicate the stopping condition.
%
% Outputs:
%   P: weather power estimator
%   vm: weather radial mean velocity estimator
%   std_v: weather spectral width
%   k: number of iterations
%   G: weather gaussian
%   S: signal spectrum after recovering the weather samples.
%   error: 0 - no errors
%          1 - error. (R0 < |R1|), although this has no teoretical meaning 
%              it could happen for some practical applications. When this 
%              happens the spectral width estimation is negative wich don't 
%              have physical meaning

time = 0;
M = length(Swn);
vs = 0.5/PRI*lambda;
v = (-M/2:(M/2 - 1))*vs/M;
pP = p1 + 1; pv = p2 + 1; %to ensure the fist iteration
k = 0;
e = 0;
error = 0;

while((pP > p1) || (pv > p2))   
    [P_old,vm_old,sigmav2_old,error] = PulsePairProcessing2(Swn,N0,0,0,M,PRI,lambda,time);    
    if error == 1
        e = 1;
        break;
    end
    
    % Recovering the weather samples.
    G = P_old/(sqrt(2*pi*(sigmav2_old)^2))*exp(-0.5*(v - vm_old).^2/(sigmav2_old)^2);    
    Sf = fftshift(Swn);
    Sf(map) = G(map);
    Swn = fftshift(Sf);
    
    [P_new,vm_new,sigmav2_new,error] = PulsePairProcessing2(Swn,N0,0,0,M,PRI,lambda,time);
    
    pP = abs(10*log10(P_new/P_old)); pv = abs(vm_new - vm_old);
    k = k+1;
    if (k > 100) || (error == 1) % maximun number of iterations or sigmav2 negative
        break;
    end;   
end

switch e
    case 0
       P = P_new; vm = vm_new; std_v = sigmav2_new; S = Swn;
    case 1
       P = P_old; vm = vm_old; std_v = sigmav2_old; S = Swn;     
end

end