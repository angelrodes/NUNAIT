function [  ] = fit_nuna_model(selectedfile)
% This funtion fit the nunatak model to sample data by randomizing
% parameters and using convergence criteria.
% 
%% Written by Ángel Rodés -- SUERC
% angelrodes@gmail.com
% 2021

fitter.ver='2.2';
disp(' ')

%% load data
load(selectedfile)

% clear any previous model data
clearvars -except samples selectedfile fitter
% close all hidden


% samples with data
sam=~isnan(samples.C);

% load constants
if exist('consts.mat', 'file') ~= 2 % create if needed
    constants
end
load('consts.mat')

% Allow user to change model parameters
reset_default_values=1;
while reset_default_values==1
    Questions = {...
        'Periglacial weathering minimum (cm/a):',...
        'Periglacial weathering maximum (cm/a):'...
        'Glacial erosion minimum (cm/a):',...
        'Glacial erosion maximum (cm/a):',...
        'Ice thinning minimum (m):',...
        'Ice thinning maximum (m):',...
        'Uplift rate minimum (m/Ma):',...
        'Uplift rate maximum (m/Ma):',...
        'Current ice surface min. deviation (m):',...
        'Current ice surface max. deviation (m):',...
        'Fit type (0=normal, 1=minimum ages, 2=maximum ages, 3=no fit)',...
        'Reset default values? (0=false, 1=true)'...
        };
    numinput=[consts.ice_free_erosion,...
        consts.glacial_erosion,...
        consts.ice_thinning,...
        consts.uplift_rate,...
        consts.present_ice_surface_uncert,...
        consts.minimum_ages,...
        0];
    dlgtitle = 'Model parameters';
    dims = 1;
    definput =cellstr(num2str(numinput(:)));
    answers = inputdlg(Questions,dlgtitle,dims,definput);
    if ~isempty(answers)
        % convert to numbers
        numanswers=cellfun(@str2num,answers)';
        % minimum erosion rates = 0.1 mm/Ma (avoid 0 in log space)
        numanswers(1:4)=max(numanswers(1:4),0.0001/1e4);
        consts.ice_free_erosion=numanswers(1:2);
        consts.glacial_erosion=numanswers(3:4);
        consts.ice_thinning=numanswers(5:6);
        consts.uplift_rate=numanswers(7:8);
        consts.present_ice_surface_uncert=numanswers(9:10);
        consts.minimum_ages=numanswers(11);
        reset_default_values=numanswers(12);
        if reset_default_values~=0
            constants
            load('consts.mat')
        end
        % check minimum and maximum
        tmpdiffs=diff(numanswers(1:10));
        if sum(tmpdiffs(1:2:9)<0)>0
            questdlg('Minimum cant be higher than maximum...', ...
                'Error', ...
                'OK','OK');
            reset_default_values=1;
        end
    end
end
save('consts.mat','consts','-v7')


% Remember that Psp=1 and sum(Pmu)=1 (muon contributions not defined yet)
% Build Production matrix
P=zeros(numel(samples.lat),numel(consts.Lmu)+1);
L=zeros(numel(samples.lat),numel(consts.Lmu)+1);
l=zeros(numel(samples.lat),1);
for n=1:numel(samples.lat)
    P(n,:)=[consts.Psp(consts.nuclides==samples.isotope(n)),consts.Pmu(consts.nuclides==samples.isotope(n),:)];
    L(n,:)=[consts.Lsp,consts.Lmu];
    l(n,:)=consts.l(consts.nuclides==samples.isotope(n));
end
density=consts.rho;
maxnmodels=consts.maxnmodels;
maxnmodelsonesigma=consts.targetnmodelsonesigma;
minmodelstoconverge=consts.minmodelstoconverge;
convergencestep=consts.convergencestep;
convergencemodels=consts.targetnmodelsonesigma;

% load climate curves
if mean(samples.lat>-55)
    make_climatecurves
else % use antarctic curves
    make_climatecurves_ant
end
load('climatecurves.mat')

%% Init variables
% wathering: erosion rate when ice-free
W=consts.ice_free_erosion(1)*(consts.ice_free_erosion(2)/consts.ice_free_erosion(1)).^rand(maxnmodels,1);

% glacial erosion
E=consts.glacial_erosion(1)*(consts.glacial_erosion(2)/consts.glacial_erosion(1)).^rand(maxnmodels,1);

% ice-thinning since maximum extension
ic=consts.ice_thinning(1)*(consts.ice_thinning(2)/consts.ice_thinning(1)).^rand(maxnmodels,1);
% D18z: ratio between d18 threshold and elevation over the ice surface
D18z=(max(climatecurves.d18O)-climatecurves.d18O(1))./ic;

% uplift rate
UR=consts.uplift_rate(1)+(consts.uplift_rate(2)-consts.uplift_rate(1)).*rand(maxnmodels,1);

% Current-ice-reference uncertainty
DELTAH=consts.present_ice_surface_uncert(1)+(consts.present_ice_surface_uncert(2)-consts.present_ice_surface_uncert(1)).*rand(maxnmodels,1);

% Chisquare
CHISQ=0.*W+NaN;
onesigma=~isnan(CHISQ);
bestmodels=~isnan(CHISQ);
models2s=0;

% Degrees of freedom (at least 1)
DOF=max(1,sum(~isnan(samples.C)) -(range(W)>0) -(range(E)>0) -(range(ic)>0) -(range(UR)>0) );


%% Run models
convergence=0;
n=0;
h = waitbar(0,'Get some coffee...');
tic
disp('Running models...')

while convergence==0
    n=n+1;
    
    % waitbar
    if mod(n-1,convergencestep)==0
        nconv=max(sum(onesigma),convergencemodels*(1-(n-minmodelstoconverge)/(maxnmodels/3))^3)*(n>=minmodelstoconverge)+sum(~isnan(CHISQ))*(n<minmodelstoconverge); % minimum convergence models
        waitbar(min(1,max(n/maxnmodels,sum(onesigma)/maxnmodelsonesigma)),h,...
            ['N_{tot.}=' num2str(n-1)...
            ' \chi^{2}_{\nu}=' num2str(min(CHISQ(~isnan(CHISQ)))/DOF,3)...
            ' N_{1\sigma}=' num2str(sum(onesigma))...
            ' N_{2\sigma}=' num2str(models2s)...
            ' N_{conv.}=' num2str(max(sum(bestmodels),round(nconv)))...
            ]);
    end
    
    % convergence methods
    if n>minmodelstoconverge && mod(n-1,convergencestep)==0 && sum(~isnan(CHISQ))>nconv
        minchi=min(CHISQ(~isnan(CHISQ)));
        minchi=minchi(1);
        onesigma=(CHISQ<minchi+DOF); % probably minchi>>DOF
        bestmodels=(CHISQ<minchi+2*DOF); % converge to 2 sigma initially
        models2s=sum(bestmodels);
        if n>maxnmodels*2/3
            bestmodels=(CHISQ<minchi+DOF); % converge to 1 sigma for the last 1/3 of the models
        end
        k=0;
        if minchi<1000*DOF
            while sum(bestmodels)<max(1,nconv) % or to the best nconv models
                k=k+1;
                bestmodels=(CHISQ<minchi+2*DOF*2^k);
%                 disp(['k=' num2str(k) ' bestmodels=' num2str(sum(bestmodels))]) % test
            end
        elseif sum(bestmodels)<max(1,nconv) % in case of terrible fitting
            sortedCHISQ=sort(CHISQ(~isnan(CHISQ))); % sort chisq to find the best models
            bestmodels=(CHISQ<sortedCHISQ(min(numel(sortedCHISQ),round(nconv+1)))); 
        end
        % converge to best models and expand 10% not to miss fitting models
        W(n:n+convergencestep)=min(W(bestmodels))*0.95*(max(W(bestmodels))/min(W(bestmodels))*1.05).^rand(convergencestep+1,1);
        W(n:n+convergencestep)=min(W(n:n+convergencestep),max(W(1:n-1))); % do not extend more than the original limits
        W(n:n+convergencestep)=max(W(n:n+convergencestep),min(W(1:n-1))); % do not extend more than the original limits
        E(n:n+convergencestep)=min(E(bestmodels))*0.95*(max(E(bestmodels))/min(E(bestmodels))*1.05).^rand(convergencestep+1,1);
        E(n:n+convergencestep)=min(E(n:n+convergencestep),max(E(1:n-1))); % do not extend more than the original limits
        E(n:n+convergencestep)=max(E(n:n+convergencestep),min(E(1:n-1))); % do not extend more than the original limits
        ic=(max(climatecurves.d18O)-climatecurves.d18O(1))./D18z;
        ic(n:n+convergencestep)=min(ic(bestmodels))*0.95*(max(ic(bestmodels))/min(ic(bestmodels))*1.05).^rand(convergencestep+1,1);
        ic(n:n+convergencestep)=min(ic(n:n+convergencestep),max(ic(1:n-1))); % do not extend more than the original limits
        ic(n:n+convergencestep)=max(ic(n:n+convergencestep),min(ic(1:n-1))); % do not extend more than the original limits
        D18z(n:n+convergencestep)=(max(climatecurves.d18O)-climatecurves.d18O(1))./ic(n:n+convergencestep);
        UR(n:n+convergencestep)=min(UR(bestmodels))*0.95+(max(UR(bestmodels))-min(UR(bestmodels)))*1.1.*rand(convergencestep+1,1);
        UR(n:n+convergencestep)=min(UR(n:n+convergencestep),max(UR(1:n-1))); % do not extend more than the original limits
        UR(n:n+convergencestep)=max(UR(n:n+convergencestep),min(UR(1:n-1))); % do not extend more than the original limits
        DELTAH(n:n+convergencestep)=min(DELTAH(bestmodels))*0.95+(max(DELTAH(bestmodels))-min(DELTAH(bestmodels)))*1.1.*rand(convergencestep+1,1);
        DELTAH(n:n+convergencestep)=min(DELTAH(n:n+convergencestep),max(DELTAH(1:n-1))); % do not extend more than the original limits
        DELTAH(n:n+convergencestep)=max(DELTAH(n:n+convergencestep),min(DELTAH(1:n-1))); % do not extend more than the original limits
    end
    
    % model concentrations 
    model=nuna_model(P,L,l,...
        W(n),E(n),UR(n),DELTAH(n),...
        samples.site_elv,samples.base_level,samples.mucont,...
        D18z(n),consts);
    
    % add uncertainty of the recent climate data (~10a)
    model.dC=model.dC+min(diff(climatecurves.age));
    
    % chi-squared
    if consts.minimum_ages==1 && max(model.C(sam)<samples.C(sam))==1
        % if ages are minimum: discard models below sample concentrations
        CHISQ(n)=NaN;
    elseif consts.minimum_ages==2 && max(model.C(sam)>samples.C(sam))==1
        % if ages are maximum: discard models above sample concentrations
        CHISQ(n)=NaN;
    elseif consts.minimum_ages==3
        % if no fitting
        CHISQ(n)=DOF;
        if n>maxnmodelsonesigma
            convergence=DOF;
            close(h)
        end
    else
        CHISQ(n)=sum((...
            (model.C(sam)-samples.C(sam))...
            ./(model.dC(sam).^2+samples.dC(sam).^2).^0.5...
            ).^2);
    end
        
    % stopping criteria
    if n>=maxnmodels-100 ||...
        sum(onesigma)>maxnmodelsonesigma &&...
            n>minmodelstoconverge
        convergence=1;
        close(h)
    end
end
disp(['...done after ' num2str(n) ' models'])
disp(['Red.Chisq=' num2str(min(CHISQ(~isnan(CHISQ)))/DOF)])
disp(['DOF=' num2str(DOF)])
disp(['N(one-sigma)=' num2str(sum(onesigma))])
toc

% recalculate ice thinning for last glacial cycle
ic=(max(climatecurves.d18O(climatecurves.age<1.2e5))-climatecurves.d18O(1))./D18z;

%% Probabilities
minchi=min(CHISQ(~isnan(CHISQ)));
minchi=minchi(1); % avoid multiple values
onesigma=(CHISQ<minchi+DOF); % because usually minchi>>DOF
best=find(CHISQ==minchi,1);
chisqpdf=@(x,dof)1./(2.^(dof/2)*gamma(dof/2)).*x.^(dof/2-1).*exp(-x/2); % chisq prob
% PROBS=chisqpdf(CHISQ,DOF); % this does not work for bad fitting values
PROBS=chisqpdf(CHISQ/DOF,1); % this is > 0 even with very low GOF


%% One sigma results
results.nmodels_onesigma=sum(onesigma);
results.totalmodels=n;
results.minchi=minchi;
results.DOF=DOF;
results.D18z=[min(D18z(onesigma)) D18z(best) max(D18z(onesigma))];
results.ic=[min(ic(onesigma)) ic(best) max(ic(onesigma))];
results.W=[min(W(onesigma)) W(best) max(W(onesigma))];
results.E=[min(E(onesigma)) E(best) max(E(onesigma))];
results.UR=[min(UR(onesigma)) UR(best) max(UR(onesigma))];
results.DELTAH=[min(DELTAH(onesigma)) DELTAH(best) max(DELTAH(onesigma))];
results.onesigma.D18z=D18z(onesigma);
results.onesigma.ic=ic(onesigma);
results.onesigma.W=W(onesigma);
results.onesigma.E=E(onesigma);
results.onesigma.UR=UR(onesigma);
results.onesigma.DELTAH=DELTAH(onesigma);
results.onesigma.PROBS=PROBS(onesigma);


%% Format results
% disp(' ')
% disp('One sigma results [min - max]:')
if range(log(results.ic))<0.90*range(log(ic)) % display only if relevant
    precision=max(2,-floor(log10(range(results.ic)/results.ic(2)))+1);
    string_ic=['[ ' num2str(round(results.ic(1))) ' - ' num2str(round(results.ic(3))) ' ] m'];
%     disp(['Ice thinning: ' string_ic])
else
    string_ic=['m'];
end
if range(log(results.W))<0.90*range(log(W))
    precision=max(2,-floor(log10(range(results.W)/results.W(2)))+1);
    string_W=['[ ' num2str(results.W(1)*1e4,precision) ' - ' num2str(results.W(3)*1e4,precision) ' ] m/Ma'];
%     disp(['Weathering: ' string_W])
else
    string_W=['m/Ma'];
end
if range(log(results.E))<0.90*range(log(E))
    precision=max(2,-floor(log10(range(results.E)/results.E(2)))+1);
    string_E=['[ ' num2str(results.E(1)*1e4/1e3,precision) ' - ' num2str(results.E(3)*1e4/1e3,precision) ' ] mm/a'];
%     disp(['Glacial erosion: ' string_E])
else
    string_E=['mm/a'];
end

if range((results.UR))<0.90*range(UR)
    precision=max(2,-floor(log10(range(results.UR)/results.UR(2)))+1);
    string_UR=['[ ' num2str(results.UR(1),precision) ' - ' num2str(results.UR(3),precision) ' ] m/Ma'];
%     disp(['Uplift rate: ' string_UR])
else
    string_UR=['m/Ma'];
end

string_DELTAH=['m'];

% disp(' ')

%% Generate fake samples to plot model results 
tic

uniquenuclides=unique(samples.isotope)'; % nuclides
nelvs=300; % number of elevations
nfs=nelvs*numel(uniquenuclides); % number of fake samples

fakesamples.lat=mean(samples.lat)*ones(nfs,1);
fakesamples.base_level=mean(samples.base_level)*ones(nfs,1);
minelv=mean(samples.base_level)-(climatecurves.d18O(1)-min(climatecurves.d18O))./min(results.D18z);
maxelv=mean(samples.base_level)+max(max(samples.elv_above_base),max(climatecurves.d18O-climatecurves.d18O(1))./min(results.D18z))*1.2; 
fakesamples.site_elv=repmat(linspace(max(0,minelv),maxelv,nelvs)',[numel(uniquenuclides) 1]);
nuctmp=repmat(uniquenuclides,[nelvs 1]);
fakesamples.isotope=nuctmp(:);
fakesamples.basin_elv=fakesamples.site_elv;

% zeros
fakesamples.mucont.value=0.*fakesamples.lat;
fakesamples.mucont.uncert=0.*fakesamples.lat;
fakesamples.C=0.*fakesamples.lat;
fakesamples.Cmin=0.*fakesamples.lat;
fakesamples.Cmax=0.*fakesamples.lat;
fakesamples.elv_above_base=0.*fakesamples.lat;

% calculate fakesample-specific data
for n=1:numel(fakesamples.lat)
    mucont = muon_contribution(fakesamples.lat(n),fakesamples.basin_elv(n),fakesamples.isotope(n));
    fakesamples.mucont.value(n)=mucont.value;
    fakesamples.mucont.uncert(n)=mucont.uncert;
    fakesamples.elv_above_base(n)=fakesamples.site_elv(n)-fakesamples.base_level(n);
end


% best model
P=zeros(numel(fakesamples.lat),numel(consts.Lmu)+1);
L=zeros(numel(fakesamples.lat),numel(consts.Lmu)+1);
l=zeros(numel(fakesamples.lat),1);
for n=1:numel(fakesamples.lat)
    P(n,:)=[consts.Psp(consts.nuclides==fakesamples.isotope(n)),consts.Pmu(consts.nuclides==fakesamples.isotope(n),:)];
    L(n,:)=[consts.Lsp,consts.Lmu];
    l(n,:)=consts.l(consts.nuclides==fakesamples.isotope(n));
end

model=nuna_model(P,L,l,...
    results.W(2),results.E(2),results.UR(2),results.DELTAH(2),...
    fakesamples.site_elv,fakesamples.base_level,fakesamples.mucont,...
    results.D18z(2),consts);

fakesamples.C=model.C;
fakesamples.dC=model.dC;
fakesamples.Cmin=model.C;
fakesamples.Cmax=model.C;

% model uncert.
h = waitbar(0,'Calculating model scatter...');
disp('Calculating model scatter...')
tic
for n=1:numel(results.onesigma.D18z)
    if mod(n,10)==0
        waitbar(n/numel(results.onesigma.D18z),h);
    end
    model=nuna_model(P,L,l,...
        results.onesigma.W(n),results.onesigma.E(n),results.onesigma.UR(n),results.onesigma.DELTAH(n),...
        fakesamples.site_elv,fakesamples.base_level,fakesamples.mucont,...
        results.onesigma.D18z(n),consts);
    fakesamples.Cmin=min(fakesamples.Cmin,model.C);
    fakesamples.Cmax=max(fakesamples.Cmax,model.C);
end
fakesamples.Cmin=fakesamples.Cmin-fakesamples.dC;
fakesamples.Cmax=fakesamples.Cmax+fakesamples.dC;
toc
close(h)

% save(samples.filename);
[filedir, filename, fileext] = fileparts (selectedfile);
newfile=fullfile(filedir,[filename '_model.mat']);
save(newfile,'-v7')



end % end of function

