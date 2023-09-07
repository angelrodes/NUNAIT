clear
close all hidden

%% Load result file
filename=[pwd '/Examples/NEG_samples_all_20220816a_sampledata_model.mat']
load(filename)

%% Select models with minimum and maximum glacial erosion rates
selected_models=[...
    find(onesigma & E==results.E(1),1,'first'),...
    find(onesigma & E==results.E(3),1,'first')];

for base_elv=unique(samples.base_level)'
    selected_samples=base_elv==samples.base_level;
    % remove repeated samples
    names=[];
    for sample_index=find(selected_samples)'
        if ~isempty(names)
            if sum(strcmp(names,samples.name(sample_index)))>0
                selected_samples(sample_index)=0;
            end
        end
        names=[names;samples.name(sample_index)];
    end

    
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

        % calculate label_elv
        label_elv=0.*samples.site_elv(selected_samples);
        [~,ind]=sort(samples.site_elv(selected_samples));
        unsorted = 1:max(ind);
        order=0.*unsorted;
        order(ind) = unsorted;
        label_elevations=interp1([min(order),max(order)],...
            [min(samples.site_elv(selected_samples)),max(samples.site_elv(selected_samples))],...
            order);
        label_elv(selected_samples)=label_elevations+20;
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
            % plot label
%             text(T(1),ELV(1),samples.name{sample_index})
%             text(1e3,ELV(1),samples.name{sample_index},...
%                 'HorizontalAlignment', 'right',...
%                 'VerticalAlignment', 'bottom')
            text(1.3e3,label_elv(sample_index),samples.name{sample_index},...
                'HorizontalAlignment', 'right',...
                'VerticalAlignment', 'middle',...
                'Color','k')
            plot([1.3e3,1e3],[label_elv(sample_index),ELV(1)],':','Color','k')

                        % Plot sample trajectory
            plot(T,ELV,'-k')
            
        end
%         plot(T,ICE_ELV,'-b')
        set(gca, 'Xdir', 'reverse')
        set(gca,'YAxisLocation','right')
        %         xlim([0 5e5])
        set(gca, 'XScale', 'log')
        xlim([1e3 1e7])
        ylim([0 max(samples.site_elv)*1.5])
        ylabel('Elevation (m)')
        xlabel('Age (a)')
        box on
        grid on
        
    end
end
