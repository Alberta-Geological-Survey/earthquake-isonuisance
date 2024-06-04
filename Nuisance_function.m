% Calculation of nuisance from ground motion

% This function calculates the probability of nuisance caused by a given
% ground motion. This code generates the probability of nuisance for CDI 3,
% though code to generate nuisance for other CDI levels is included but
% commented out.

% 'GM' is an array of ground motions with form of [ground motion,
% realization].
% The variable 'motion' should either be 'PGA' or 'PGV'. It determines the
% equations used to calculate nuisance.

% This function perturbs the nuisance. The perturbation is consistent along
% the first dimension of 'GM', but varies along the second dimension. Thus,
% a different perturbation is applied to each realization.

% Code modified from:
% Schultz, R., Beroza, G. C., & Ellsworth, W. L. (2021). A strategy for
% choosing red-light thresholds to manage hydraulic fracturing induced
% seismicity in North America. Journal of Geophysical Research: Solid
% Earth, 126

% Based on nuisance eqations from:
% Schultz, R, V Quitoriano, DJ Wald, GC Beroza. (2021). Quantifying
% nuisance ground motion thresholds for induced earthquakes - Earthquake
% Spectra.

function Mu = Nuisance_function(GM, motion)

    % Define some variables to perterb the nuisance
    dE1 = normrnd(0,1,1,size(GM,2));
    dE2 = normrnd(0,1,1,size(GM,2));

    if strcmp(motion, 'PGA') % For PGA values 

        %%%%%%%% For CDI 2
        
%         % Define variables
%         B2=[4.9965 3.2018];
%         dB2=[+1.98e-2 +4.86e-3 +6.21e-3];
%         
%         % Calculate the nuisance factors
%         B21 = B2(1) + (dE1 * sqrt(dB2(1)));
%         B22 = B2(2) + (dE2 * sqrt(dB2(2))) + (dE1 * (dB2(3) / sqrt(dB2(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM2 = log10(GM);
%         
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B21 + (B22 .* GM2))) ./ (1 + (exp(B21 + (B22 .* GM2))));
        

        %%%%%%%% For CDI 3

        % Define more variables
        B3=[2.1147 2.2630];
        dB3=[+1.74e-2 +2.10e-3 +4.30e-3];

        % Calculate the nuisance factors
        B31 = B3(1) + (dE1 * sqrt(dB3(1)));
        B32 = B3(2) + (dE2 * sqrt(dB3(2))) + (dE1 * (dB3(3) / sqrt(dB3(1))));

        % Convert ground motion to logarithmic scale
        GM3 = log10(GM);

        % Calculate the probabilities of nuisance
        Mu = (exp(B31 + (B32 .* GM3))) ./ (1 + (exp(B31 + (B32 .* GM3))));
        

        %%%%%%%%% For CDI 4
        
%         % Define variables
%         B4=[0.2266 1.8061];
%         dB4=[+9.36e-3 +9.16e-4 +1.09e-3];
%         
%         % Calculate the nuisance factors
%         B41 = B4(1) + (dE1 * sqrt(dB4(1)));
%         B42 = B4(2) + (dE2 * sqrt(dB4(2))) + (dE1 * (dB4(3) / sqrt(dB4(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM4 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B41 + (B42 .* GM4))) ./ (1 + (exp(B41 + (B42 .* GM4))));


        %%%%%%%%% For CDI 5
        
%         % Define variables
%         B5=[-0.8126 1.7358];
%         dB5=[+6.06e-3 +1.60e-4 -6.67e-5];
%         
%         % Calculate the nuisance factors
%         B51 = B5(1) + (dE1 * sqrt(dB5(1)));
%         B52 = B5(2) + (dE2 * sqrt(dB5(2))) + (dE1 * (dB5(3) / sqrt(dB5(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM5 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B51 + (B52 .* GM5))) ./ (1 + (exp(B51 + (B52 .* GM5))));
        

        %%%%%%%%% For CDI 6
        
%         % Define variables
%         B6=[-1.1618 1.7703];
%         dB6=[+4.58e-3 +9.63e-5 -1.28e-4];
%         
%         % Calculate the nuisance factors
%         B61 = B6(1) + (dE1 * sqrt(dB6(1)));
%         B62 = B6(2) + (dE2 * sqrt(dB6(2))) + (dE1 * (dB6(3) / sqrt(dB6(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM6 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B61 + (B62 .* GM6))) ./ (1 + (exp(B61 + (B62 .* GM6))));
        
        
    elseif strcmp(motion, 'PGV') % For PGV values 
        %%%%%%%% For CDI 2
        
%         % Define variables
%         B2=[11.2336 3.8771];
%         dB2=[+6.05e-1 +4.47e-2 +1.63e-1];
%         
%         % Calculate the nuisance factors
%         B21 = B2(1) + (dE1 * sqrt(dB2(1)));
%         B22 = B2(2) + (dE2 * sqrt(dB2(2))) + (dE1 * (dB2(3) / sqrt(dB2(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM2 = log10(GM);
%         
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B21 + (B22 .* GM2))) ./ (1 + (exp(B21 + (B22 .* GM2))));
        
        
        % %%%%%%%% For CDI 3

        % Define variables
        B3=[5.9376 2.5875];
        dB3=[+2.48e-1 +1.84e-2 +6.68e-2];

        % Calculate the nuisance factors
        B31 = B3(1) + (dE1 * sqrt(dB3(1)));
        B32 = B3(2) + (dE2 * sqrt(dB3(2))) + (dE1 * (dB3(3) / sqrt(dB3(1))));

        % Adjust ground motion for equation
        GM3 = log10(GM);

        % Calculate the probabilities of nuisance
        Mu = (exp(B31 + (B32 .* GM3))) ./ (1 + (exp(B31 + (B32 .* GM3))));

        % %%%%%%%% For CDI 4
        
%         % Define variables
%         B4=[3.4387 2.1331];
%         dB4=[+8.10e-2 +4.23e-3 +1.78e-2];
%         
%         % Calculate the nuisance factors
%         B41 = B4(1) + (dE1 * sqrt(dB4(1)));
%         B42 = B4(2) + (dE2 * sqrt(dB4(2))) + (dE1 * (dB4(3) / sqrt(dB4(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM4 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B41 + (B42 .* GM4))) ./ (1 + (exp(B41 + (B42 .* GM4))));
        

        %%%%%%%% For CDI 5

%         % Define variables
%         B5=[2.3712 2.0828];
%         dB5=[+2.37e-2 +6.63e-4 +3.16e-3];
%         
%         % Calculate the nuisance factors
%         B51 = B5(1) + (dE1 * sqrt(dB5(1)));
%         B52 = B5(2) + (dE2 * sqrt(dB5(2))) + (dE1 * (dB5(3) / sqrt(dB5(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM5 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B51 + (B52 .* GM5))) ./ (1 + (exp(B51 + (B52 .* GM5))));


        %%%%%%%% For CDI 6
        
        % Define variables
%         B6=[2.1427 2.1260];
%         dB6=[+2.12e-2 +5.18e-4 +2.42e-3];
%         
%         % Calculate the nuisance factors
%         B61 = B6(1) + (dE1 * sqrt(dB6(1)));
%         B62 = B6(2) + (dE2 * sqrt(dB6(2))) + (dE1 * (dB6(3) / sqrt(dB6(1))));
%         
%         % Convert ground motion to logarithmic scale
%         GM6 = log10(GM);
% 
%         % Calculate the probabilities of nuisance
%         Mu = (exp(B61 + (B62 .* GM6))) ./ (1 + (exp(B61 + (B62 .* GM6))));
        

    else % if motion is neither 'PGA' nor 'PGV'
        error("Nuisance function requires motion to be specified as 'PGA' or 'PGV'")
    end
    

end 