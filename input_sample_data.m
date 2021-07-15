function [ ] = input_sample_data(selectedfile)
% This function is used to ask the user for input data

fid = fopen(selectedfile);
%% load data
% Bug: In octave textscan is only taking the first sample!!!
mydata = textscan(fid, '%s %f %f %f %f %f %f %f %f ',...
    'HeaderLines', 1,'Delimiter',',');
fclose(fid);
samples.name=mydata{1};
samples.lat=mydata{2};
samples.site_elv=mydata{3};
samples.basin_elv=samples.site_elv; % this is included for sediments
samples.isotope=mydata{4};
samples.base_level=mydata{5};
samples.apparent_years=mydata{6};
samples.dapparent_years=mydata{7};

% Check common errors
if median(samples.apparent_years(~isnan(samples.apparent_years)))<1e3
    warning('Make sure that your apparent ages are in years, not in ka or Ma!')
end
if max(abs(samples.lat))>90
    warning('Latitude must be between -90 and 90!')
end


%% load constants
if exist('consts.mat', 'file') ~= 2 % create if needed
    constants
end
load('consts.mat')


%% init variables
samples.l=0.*samples.lat;
samples.mucont.value=0.*samples.lat;
samples.mucont.uncert=0.*samples.lat;
samples.C=0.*samples.lat;
samples.dC=0.*samples.lat;
samples.elv_above_base=0.*samples.lat;




%% calculate
for n=1:numel(samples.lat)
    if sum(consts.nuclides==samples.isotope(n))>0
        samples.l(n)=consts.l(consts.nuclides==samples.isotope(n));
    else
        warning(['Nuclide ' num2str(samples.isotope(n)) ' not recognized by ' mfilename])
    end
    mucont = muon_contribution(samples.lat(n),samples.basin_elv(n),samples.isotope(n));
    samples.mucont.value(n)=mucont.value;
    samples.mucont.uncert(n)=mucont.uncert;
    
    if samples.apparent_years(n)>0 % igonere ages with values = 0
        samples.C(n)=round(1/samples.l(n)*(1-exp(-samples.l(n)*samples.apparent_years(n))));
        samples.dC(n)=round(samples.dapparent_years(n)*exp(-samples.l(n)*samples.apparent_years(n)));
    else
        samples.C(n)=NaN;
        samples.dC(n)=NaN;
    end
    
    samples.elv_above_base(n)=samples.site_elv(n)-samples.base_level(n);
    % Check errors
    if samples.site_elv(n)<samples.base_level(n)
        warning([samples.name{n} ' is under the ice surface!'])
    end
    if isnan(samples.C(n)) 
        warning([samples.name{n} ' has no cosmo data!'])
    end
end

display_sample_data(samples)


[filedir, filename, fileext] = fileparts (selectedfile);
samples.filename=filename;

newfile=fullfile(filedir,[filename '_sampledata.mat']); 

save(newfile,'samples','-v7')

end % end of function

