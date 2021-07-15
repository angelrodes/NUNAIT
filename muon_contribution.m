function [ mucont ] = muon_contribution(lat,elv,nuclide)
% This function generates the muon contribution and its ucertainty based on
% latitude and elevation for a given nuclide. If either lat or elv is not a
% number, a global average is given.
% A single value of latitude and elevation is used to calculate the
% contribution of muons to the total surface production of Be-10. All other
% productions are scaled accordingly.
% Angel Rodes, 2020

%% version
mucont.ver='1.4';

%% Load constats
load('consts.mat')

mucont.consts=consts;

% selct nuclide
nucidx=(nuclide==consts.nuclides);

if sum(nucidx)==0
    warning(['Nuclide ' num2str(nuclide) ' not included in consts'])
end

% These approximations are based on the production at 1678 sites
% equally distributed on land areas according to
% ETOPO1_Bed_g_geotiff.tif (Eakins et al., 2012)
% and calculated using
% P_mu_total_alpha1.m and stone2000.m from CRONUS calculators v2.3
% (Balco, 2017 and Balco, 2008)

if sum(isnan([lat,elv]))>0
    % Surface muon contribution to Be-10 production (global)
    C_mu_P_10=0.0126; %1.26%
    dC_mu_P_10=0.00044; % Scatter over land areas
    % The scatter of this value over the Earth surface introduces a large
    % uncertainty (35%) on the muon production. This could be reduced by
    % considering a surface muon contribution based on altitude and latitude.
else
    % Surface muon contribution to Be-10 production (site)
    C_mu_P_10 = (1.29+lat/900+1.056*exp(-((lat+1)/30.31).^2)) .* (0.1+0.9*exp(-elv/2000)) / 100 ;
    % This formula fits the original data within a 5%
    dC_mu_P_10=0.05*C_mu_P_10;
end

% Sources of uncertainties on the muon contributions:
%    - Total surface production rate uncertainties: 8 to 12% using LSD scaling
%        (see Phillips et al., 2016 https://doi.org/10.1016/j.quageo.2015.09.006)
%    - Uncertainty on the surface P26/P10 ratio: ~10% according to the
%        calibration dataset from Borchers et al., 2016
%        https://doi.org/10.1016/j.quageo.2015.01.009
%    - Scatter of the muon production prediceted by P_mu_total_alpha1.m:
%          5-13%
%         (see  Balco 2017 https://doi.org/10.1016/j.quageo.2017.02.001)
%    - Scatter of the muon P26/P10: ~20% according to data scatter in Balco (2017)
% Considering all these unceratinties, a minimum uncertainty of 20% on muon produced
% cosmonuclides is considered in the models:
dC_mu_P_10=max(dC_mu_P_10,0.2*C_mu_P_10);


% Ratio between Al-26 and Be-10 muon contributions to total proction rates
% (unitless, %/%)
C_mu_P2610=1.4587;

% Ratio between Cl-36 and Be-10 muon contributions to total proction rates
% (unitless, %/%). Consistent with Heisinger & Nolte (2000)
C_mu_P3610=3.2720;

% Ratio between Ne-21 and Be-10 muon contributions to total proction rates
% (unitless, %/%). Consistent with Balco & Shuster (2009)
C_mu_P2110=4.086;

% Ratio between He-3 and Be-10 muon contributions to total proction rates
% (unitless, %/%). Consistent with Blard et al. (2013)
C_mu_P310=1;

% Ratio between C-14 and Be-10 muon contributions to total proction rates
% (unitless, %/%). Consistent with Heisinger & Nolte (2000)
C_mu_P1410=8.2767;

if nuclide==10
    mucont.value=C_mu_P_10;
    mucont.uncert=dC_mu_P_10;
elseif nuclide==26
    mucont.value=C_mu_P_10*C_mu_P2610;
    mucont.uncert=dC_mu_P_10*C_mu_P2610;
elseif nuclide==36
    mucont.value=C_mu_P_10*C_mu_P3610;
    mucont.uncert=dC_mu_P_10*C_mu_P3610;
elseif nuclide==14
    mucont.value=C_mu_P_10*C_mu_P1410;
    mucont.uncert=dC_mu_P_10*C_mu_P1410;
elseif nuclide==21
    mucont.value=C_mu_P_10*C_mu_P2110;
    mucont.uncert=dC_mu_P_10*C_mu_P2110;
elseif nuclide==3
    mucont.value=C_mu_P_10*C_mu_P310;
    mucont.uncert=dC_mu_P_10*C_mu_P310;
else
    warning(['Nuclide <' num2str(nuclide) '> not recognized by ' mfilename])
    mucont.value= NaN;
    mucont.uncert=NaN;
end

end