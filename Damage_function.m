% Calculation of damage from ground motion

% This function calculates the probability of damage caused by a given
% ground motion. This code generates the probability of damage for DSI1,
% note that code to generate nuisance for other DSI levels is included but
% commented out.

% 'GM' is an array of ground motions with form of [ground motion,
% realization].

% This function perturbs the damage, by perturbing mean (mu) and standard
% deviation (sigma) of a distribution function used to define the damage.
% The perturbation is consistent along the first dimension of 'GM', but
% varies along the second dimension. Thus, a different perturbation is
% applied to each realization.

% *** Function expects ground motion in PGA ***


% Code modified from:
% Schultz, R., Beroza, G. C., & Ellsworth, W. L. (2021). A strategy for
% choosing red-light thresholds to manage hydraulic fracturing induced
% seismicity in North America. Journal of Geophysical Research: Solid
% Earth, 126

% Based on equations from:
% Federal Emergency Management Agency (FEMA). (2015). “Hazus-MH 2.1” in
% Multi-Hazard Loss Estimation Methodology Technical and User Manuals
% (Federal Emergency Management Agency, 2015).


function P = Damage_function(GM2)


    % Convert ***PGA*** values of ground motion into a fraction of
    % gravitational acceleration
    GM=GM2/9.8;
    
    % Initialize output array
    P = zeros(size(GM));
    
    % Create variables for random perturbations
    dE1 = normrnd(0,1,1,size(GM,2));
    dE2 = normrnd(0,1,1,size(GM,2));
    
    
    %%%%%%% For DS1 % Slight/minor
    
    % Define variables for DS1
    u1= [0.20, 0.653];
    du1=[0.00, 0.014];

    % Calculate interim values
    u1_1 = u1(1) + (dE1 .* du1(1));
    u1_2 = u1(2) + (dE2 .* du1(2));
    
    % Populate output array
    for i = 1:size(GM,2)
        P(:,i,:) = normcdf(log(GM(:,i,:)), log(u1_1(i)), u1_2(i));
    end


    %%%%%%% For DS2  % Moderate

%     % Define variables for DS2
%     u2= [0.40, 0.672];
%     du2=[0.00, 0.009];
%     
%     % Calculate interim values
%     u2_1 = u2(1) + (dE1 .* du2(1));
%     u2_2 = u2(2) + (dE2 .* du2(2));
%     
%     % Populate output array
%     for i=1:length(GM)
%         P(:,i,:) = normcdf(log(GM(:,i,:)), log(u2_1(i)), u2_2(i));
%     end 

    %%%%%%% For DS3  % Extensive

%     % Define variables for DS3
%     u3= [0.80, 0.668];
%     du3=[0.00, 0.014];
%     
%     % Calculate interim values
%     u3_1 = u3(1) + (dE1 .* du3(1));
%     u3_2 = u3(2) + (dE2 .* du3(2));
%     
%     % Populate output array
%     for i=1:length(GM)
%         P(:,i,:) = normcdf(log(GM(:,i,:)), log(u3_1(i)), u3_2(i));
%     end
    

    %%%%%%% For DS4  % Complete
    
%     % Define variables for DS4
%     u4= [1.60, 0.668];
%     du4=[0.00, 0.014];
%     
%     % Calculate interim values
%     u4_1 = u4(1) + (dE1 .* du4(1));
%     u4_2 = u4(2) + (dE2 .* du4(2));
%     
%     % Populate output array
%     for i=1:length(GM)
%       P(:,i,:) = normcdf(log(GM(:,i,:)), log(u4_1(i)), u4_2(i));
%     end 

end 