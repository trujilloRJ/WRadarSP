function [P,vm,std,error] = GMAP(z,K,std_c,PRI,lambda,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%   Siggia A. D., Gaussian Model Adaptive Processing(GMAP) for improved
%       ground clutter cancellation and moment calculation
%--------------------------------------------------------------------------
% Description:
%   This is the main function that implements GMAP in frequency domain, an
%   algorithm that remove clutter components and try to recover the
%   weather components that were removed in the filter stage. It assumes
%   that both the weather and the clutter have gaussian shape spectrum, and
%   that the clutter spectral width is given.
%--------------------------------------------------------------------------
% Inputs:
%   z: complex random vector (IQ components echo signal)
%   K: number of segments in which the signal will be divided for
%      estimating its spectrum. 
%   std_c: clutter spectral width [m/s]
%   PRI: pulse repetition interval [s] (sample time for samples from the 
%       same range bin)
%   lambda: carrier wavelength [m]
%   w: window to select data
%       RECTANGULAR
%       HAMMING
%       BLACKMAN
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

time = 0; % all processing will be in frequency domain.
p1 = 0.001;
p2 = 0.001;
vs = 0.5/PRI*lambda;
M = length(z)/K;
if M < 32 
    Mfit = 32;
end
v = (-M/2:(M/2 - 1))*vs/M;
error = 0;
t = 0;
N = 500;

while(t <= 1)
%---------------- Power spectrum estimation. Step 1 ----------------------% 
    S = WelchSpectraEstimation(z,K,w);
    
%---------------- Noise level determination. Step 2 ----------------------%
    [N0,THM,k] = NoiseLevelDetermination(S,N,K);    
    
%---------------- Removing clutter components. Step 3 --------------------%
    % It's necessary to adjust the clutter spectral width because of the 
    % spectral leakage produced by the window that is used.
    [SwF, sigmaw] = FitWindow(Mfit,w);  
    sigmav2c = ((std_c)^2 + (sigmaw*vs)^2);
    [Ssc,vsc,G,C,NOCLUTTER] = RemoveClutter(v,S,sigmav2c,THM,N0);
    % If there isn't clutter it's not necessary to remove any clutter
    % component and the estimation is performed with a rectangular window
    
%---------------- Moments Estimation. Step 4 -----------------------------%    
    if NOCLUTTER == 1
        w = 'RECTANGULAR';
        S = WelchSpectraEstimation(z,K,w);
        [P,vm,std,error] = PulsePairProcessing2(S,THM,0,0,M,PRI,lambda,time);
        Sw = S;
        break;
    else
        NROZERO = length(find(Ssc));
        if NROZERO == 0             % No weather detected
            disp('No se detecto fenomeno');
            Sw = S;                 
            P = 0; vm = 0; std = 0;
        elseif NROZERO == 1;        % Weather is represented by only one sample
            P = Ssc(find(Ssc))*vs/M;
            Sw = Ssc;
            vm = find(Ssc)*vs/M;
            if vm > vs/2
                vm = vm - vs;
            end
            std = 0;
        else
            map = find(Gc > THM);
            [P,vm,std,k,G,Sw,error] = iterativePPP(Ssc,N0,PRI,lambda,map,p1,p2);
        end
        
%---------------- Choosing the appropriate window. Step 5 -----------------------------%    
        if P == 0
            CSR_est = 10*log10(C/(N0*length(Ssc))); % If there si no weather, it is replaced by noise. 
            break;
        else
            if C > P
                CSR_est = 10*log10(C/P); %dB
            else
                CSR_est = 0;
            end
        end
        % Window selection
        if t == 0
            if CSR_est < 15
                if CSR_est > 2.5
                    break; % Hamming is OK
                else
                    P0 = P; vm0 = vm; sigmav20 = std;
                    w = 'RECTANGULAR';   % Repeat with rectangular
                end
            else
                CSR0 = CSR_est;
                P0 = P; vm0 = vm; sigmav20 = std;
                w = 'BLACKMAN'; % Repeat with blackman
            end
        else
            if strcmp(w,'RECTANGULAR') == 1
                if CSR_est > 1
                    P = P0; vm = vm0; std = sigmav20; % Haming is OK
                end
            else
                if CSR0 < 40
                    if CSR_est < 20
                        P = P0; vm = vm0; std = sigmav20; % Haming is OK
                    end
                end
            end
        end
        t = t+1;
    end

end