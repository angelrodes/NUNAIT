 function [ ] = START(  )
%% Nunatak ice-thinning model
%% Written by Ángel Rodés -- SUERC
% angelrodes@gmail.com
% 2021

%% TODO list
% Solve textscan in input_sample_data in octave (only first sample taken)
%    - explore dlmread instead of textscan
%    - maybe ignore sample names (we don't really use them beyond listing)
%        dlmread(filename,',',1,1)

%% Ask what to do
btn = questdlg ('Select', 'NUNAIT calculator', 'Run simulation', 'Display results', 'Run simulation');

%% remove previous plots
close all hidden

if strcmp(btn,'Run simulation')
    %% Input data
    %     [file,path] = uigetfile({'*.csv','*.CSV','*.mat','*.MAT'}, 'Select input file or production rates');
    [file,path] = uigetfile({'*.*'}, 'Select input file (csv) or production rates (mat)');
    selectedfile = fullfile(path,file);
    [filedir, filename, fileext] = fileparts (selectedfile);
    
    if ~strcmp(fileext,'.mat') && ~strcmp(fileext,'.MAT')
        input_sample_data(selectedfile); % load data and calculate production rates
        newfile=fullfile(filedir,[filename '_sampledata.mat']); % name of the file generated by input_sample_data
        selectedfile=newfile;
    else
        load(selectedfile) % load data and production rates
    end
    
    %% Fit model
    disp(['Fitting: ' selectedfile '']) 
    plot_samples(selectedfile); % run model
    fit_nuna_model(selectedfile); % run model
    [filedir, filename, fileext] = fileparts (selectedfile);
    selectedfile=fullfile(filedir,[filename '_model.mat']); 
    
else
    [file,path] = uigetfile({'*model.mat','*model.MAT'}, 'Select output file'); 
    selectedfile = fullfile(path,file);
    disp(['File: ' selectedfile])
end

%% Display results
display_results(selectedfile);
plot_results(selectedfile);


%% todo
% allow several files at a time by taking file salection out of input
% samples
% display credits and VERSION CONTROL (consts (general),climate data, fitting script, muon contribution, nuna_model)
% display result text

end

