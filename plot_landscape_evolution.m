function [  ] = plot_landscape_evolution(filename)
%% Calculates the evolution of the elevation of the sampling surfaces 
%% according to the models showing minimum and maximum glacial erosion 
%% rates (extreme scenarios).
%
%% Angel Rodes, 2022
%% www.angelrodes.com

%% Load result file
% filename='/home/angel/Desktop/Angel-work_2022/Roberts_Lane_2022/NUNAIT-main_NEG_20220816/Examples/NEG_samples_all_20220816a_sampledata_model.mat'
load(filename)

%% Select models with minimum and maximum glacial erosion rates
selected_models=[...
    find(onesigma & E==results.E(1),1,'first'),...
    find(onesigma & E==results.E(3),1,'first')];

for base_elv=unique(samples.base_level)'
    selected_samples=base_elv==samples.base_level;
    
    for model_index=selected_models
        figure('units','normalized','outerposition',[0 0 1 1],'Name','Surfaces')
        hold on
        %         title(['base\_level=' num2str(base_elv) 'm ; \epsilon=' num2str(E(model_index)*1e4,4) 'mm/ka'])
        title(['Glacial erosion rate=' num2str(E(model_index)*1e4/1000,4) 'mm/a'])
        
        T=climatecurves.age;
        dT=diff([climatecurves.age,4543e6]);
        ICE_ELV=base_elv+DELTAH(model_index)+(climatecurves.d18O-climatecurves.d18O(1))/D18z(model_index);
        plot(T,ICE_ELV,'-b')
        plot([0,0],[0,0],'-k')
        legend('Ice surface','Sampled surfaces','AutoUpdate','off')
        for sample_index=find(selected_samples)'
            ELV=0.*T;
            ELV(1)=samples.site_elv(sample_index);
            for n=2:numel(T)
                BL=samples.base_level(sample_index)+DELTAH(model_index);
                D=D18z(model_index)*(ELV(n-1)-BL)+climatecurves.d18O(1);% threshold
                % glaciated?
                G=(climatecurves.d18O(n)+climatecurves.age(n)/1e6*UR(model_index)*D18z(model_index)...
                    >D);
                % erosion during this step
                Zi=E(model_index).*dT(n).*G+W(model_index).*dT(n).*~G;
                % elevations
                ELV(n)=Zi/100+ELV(n-1);
            end
            % Plot sample trajectory
            plot(T,ELV,'-k')
            text(T(1),ELV(1),samples.name{sample_index})
            text(1e3,ELV(1),samples.name{sample_index},...
                'HorizontalAlignment', 'right',...
                'VerticalAlignment', 'bottom')
            
        end
%         plot(T,ICE_ELV,'-b')
        set(gca, 'Xdir', 'reverse')
        set(gca,'YAxisLocation','right')
        set(gca, 'XScale', 'log')
        xlim([1e3 1e7])
        ylim([0 max(samples.site_elv)*1.5])
        ylabel('Elevation (m)')
        xlabel('Age (a)')
        box on
        grid on
        
    end
end

end