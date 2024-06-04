% Calculating Ground Motion

% This function has been vectorized to allow calculation of ground motion
% at multiple locations given earthquakes of multiple magnitudes.

% 'Mag' can be a single value, or an array of magnitudes. If using an
% array, the first dimension must match the length of distances and Vs30.
% The variable 'distances' is an array containing the distances for GM
% calculation. The first dimension contains the distances, and second
% dimension contains different realizations.
% 'Vs30' is a value, vector, or array of near-surface velocity
% corresponding to the locations described by 'distances'. If defined as a
% vector, it's length should match the length of the first dimension of
% 'distances'. If defined as an array, it's size should match 'distances'.
% The variable 'motion' selects the units in which the results are
% returned: 'PGA' returns ground motion as Peak Ground Acceleration
% in units of cm/s2 in a logarithmic scale; 'PGV' returns ground motion as
% Peak Ground Velocity in units of cm/s in a logarithmic scale.

% This function perturbs the GM. The perturbation is consistent along
% the first dimension of 'distances', but varies along the second
% dimension. Thus, a different perturbation is applied to each realization.


% Code modified from:
% Schultz, R., Beroza, G. C., & Ellsworth, W. L. (2021). A strategy for
% choosing red-light thresholds to manage hydraulic fracturing induced
% seismicity in North America. Journal of Geophysical Research: Solid
% Earth, 126

% Based on the GMPEs from:
% Schultz, R. and Nanometrics Inc. (2019): Initial seismic hazard
% assessment for the induced earthquakes near Fox Creek, Alberta (between
% January 2013 and January 2016); Alberta Energy Regulator, Alberta
% Geological Survey, AER/AGS Special Report 104, 115 p.

function GM_cat = GM_calc(Mag, distances, Vs30, motion)

% 'motion' should either be 'PGA' or 'PGV'

M=Mag;
Rh=distances;
dEg = normrnd(0,1,1,size(Rh,2)); % For perturbing the GM results, one value per realization

if strcmp(motion, 'PGA') % For PGA
      c0=-2.376;
      dc0=-0.212;
      c1=1.818;
      c2=-0.1153;
      c3=-1.752;
      dc3=1.992;
      ea=0.28;
      ei=0.24;
      et=0.37;
      
      c=-0.22;
      Vc=638.08;
      Vr=760;
      
elseif strcmp(motion, 'PGV') % For PGV
      c0=-4.151;
      dc0=0.00;
      c1=1.762;
      c2=-0.09509;
      c3=-1.669;
      dc3=1.582;
      ea=0.27;
      ei=0.19;
      et=0.33;
      
      c=0.18;
      Vc=515.9;
      Vr=760;
      
else % if motion is neither 'PGA' nor 'PGV'
    error("GMPE function requires motion to be specified as 'PGA' or 'PGV'")
end


% Compute effective hypocentre
heff=max(ones(size(M)),10.^(-1.72+0.43.*M));
R=sqrt(Rh.^2+heff.^2);

% Compute GMPE for close distances
Y=(c0+dc0)+(c1.*M)+(c2.*M.^2)+c3.*log10(R);

% Modify GMPE values for mid range distances
I=(R>=70)&(R<140);
Y(I)=Y(I)+dc3.*log10(R(I)./70);

% Modify GMPE values for long range distances
I=(R>=140);
Y(I)=Y(I)+dc3.*log10(140./70);


%%% Add site amplificaiton effects

% Calculate effect
if size(Vs30, 2) == 1 % In case a single vector of Vs30 is provided
    Vs30 = repmat(Vs30, [1,size(Y, 2)]);
end
Fs=zeros(size(Vs30));
I=Vs30 < Vc;
Fs(I) = c * log(Vs30(I)/Vc) - c * log(Vr/Vc);

% Apply site amplification
Y=Y+Fs;

% Apply perturbation.
GM_cat = Y + et * dEg;
  
end

