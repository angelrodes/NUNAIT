function [ out ] = nuna_model(Production_rates,Attenuation_lengths,decay_constant,...
    weathering,glacial_erosion_rate,uplift_rate,delta_h,...
    site_elv,base_level,mucont,...
    D18z,consts)
% NUNATAK ACCUMULATION MODEL
% Production_rates and Attenuation_lengths are expected to be matrices
% density is expected to be a single value
% weathering,glacial_erosion_rate,uplift_rate,delta_h and D18z are expected to be single values
% decay_constant, zsamples,zice and mucont are expected to be columns with the same lengths 
% This model CONSIDERS production under ice
%% Angel Rodes, 2021


%% Load climate curves
load('climatecurves.mat')
% remove data not needed to speed up computing
if uplift_rate==0
    minDneeded=D18z*min(site_elv-base_level)+climatecurves.d18O(1);
    needed=(climatecurves.d18O>=minDneeded);
    needed(1:2)=1;
    climatecurves.d18O(~needed)=minDneeded;
    remove=(diff([climatecurves.d18O,1])==0);
    climatecurves.d18O=climatecurves.d18O(~remove);
    climatecurves.age=climatecurves.age(~remove);
end

%% Simplify names
l=decay_constant+log(2)/(100*4543e6); % add a small decay (100 times earth age) to simplify math with stable isotopes
% P=Production_rates;
% L=Attenuation_lengths;

%% Check inputs 
if numel(Production_rates)~=numel(Attenuation_lengths) || ...
        numel(site_elv)~=numel(base_level) || ...
        numel(weathering)>1 || numel(glacial_erosion_rate)>1
    warning('wrong size of input parameters')
end

%% load constants
% load('consts.mat')
density=consts.rho;
densityice=consts.rhoice;


%% Make input matrices 
[sample_index,p_index,climate_index]=ndgrid(1:numel(site_elv),1:size(Production_rates,2),1:numel(climatecurves.age));


% Dimensions:
% 1 samples: zsamples,zice and mucont
% 2 production rates
% 3 climate curves

T=climatecurves.age(climate_index); 
HS=site_elv(sample_index);
BL=base_level(sample_index)+delta_h;
MC=mucont.value(sample_index);
dMC=mucont.uncert(sample_index);
lmatrix=l(sample_index);
D=D18z*(HS-BL)+climatecurves.d18O(1);% thresholds
E=glacial_erosion_rate;
W=weathering;
P=repmat(Production_rates,[1,1,numel(climatecurves.age)]);
L=repmat(Attenuation_lengths,[1,1,numel(climatecurves.age)]);


% Define step lengths (years)
climatecurves.dage=diff([climatecurves.age,4543e6]);
dT=climatecurves.dage(climate_index);

% glaciated times
G=(climatecurves.d18O(climate_index)+...
    climatecurves.age(climate_index)/1e6*uplift_rate*D18z...
    >D);

% erosion during each step
Zi=E.*dT.*G+W.*dT.*~G;

% depths
Z=cumsum(Zi,numel(size(Zi)))-Zi;

% depths under ice
ZICE=max(0*D,...
    (climatecurves.d18O(climate_index)+...
    climatecurves.age(climate_index)/1e6*uplift_rate*D18z...
    -D)/D18z)*100; % cm


% 
% size_P=size(P) % TESTING ONLY
% size_L=size(L)
% size_lmatrix=size(lmatrix)
% size_W=size(W)
% size_Z=size(Z)
% P(:,:,35)
% caca

% Accumulation model
Cii=...
    P./(lmatrix+W.*density./L).*...
    exp(-(Z.*density+ZICE.*densityice)./L).*...
    (1-exp(-(lmatrix+W.*density./L).*dT)).*...
    exp(-lmatrix.*T);

%% Make output arrays

% Sum for all time steps
Ci=sum(Cii,numel(size(Cii)));

% Concentration produced by spallation (assuming that first production is spallation)
Csp=Ci(:,1);

% Concentration produced by muons (assuming that first production is spallation)
Cmu=(sum(Ci,numel(size(Ci)))-Csp).*MC(:,1,1);
dCmu=(sum(Ci,numel(size(Ci)))-Csp).*dMC(:,1,1);

% final concentrations and uncertainties
out.C=Csp+Cmu;
out.dC=dCmu;

% altitude over ice
out.alt_over_ice=(HS(:,1,1)-BL(:,1,1)); 

end

