function [  ] = plot_individual_profiles(filename)
% This function plots results generated by fit_nuna_model
% Same output as the last plots generated by plot_results.m but samples
% sharing the same base_level are plotted in separated graphs.
%% Angel Rodes, 2023


%% load data
load(filename)

% load climate curves
if mean(samples.lat>-55)
    make_climatecurves
else % use antarctic curves
    make_climatecurves_ant
end
load('climatecurves.mat')

% load constants
if exist('consts.mat', 'file') ~= 2 % create if needed
    constants
end
load('consts.mat')


%% Plot apparent_age(concentration) individual profiles
uniquenuclides=unique(samples.isotope)';


for nuc=uniquenuclides
    figure('units','normalized','outerposition',[0 0 0.75 0.75],'Name',[consts.nuclidesstring{consts.nuclides==nuc} ' apparent age profiles'])
    nplot=0;
    for this_elevation=sort(unique(samples.base_level),'descend')'
        nplot=nplot+1;
        subplot(1,numel(sort(unique(samples.base_level),'descend')'),nplot)
        hold on
        selfake=(fakesamples.isotope==nuc);
        selsamp=(samples.isotope==nuc & this_elevation==samples.base_level);
        lambda=mean(samples.l(selsamp));
        apparent_age=@(concentration) log(1-min(max(concentration,0)*lambda,1))/(-lambda); % a
        base=this_elevation;
        minx=100; % plot minimum 100 year
        %     minx=max(1,min(samples.apparent_years)/2);
        x1=fakesamples.Cmin(selfake);
        x2=fakesamples.C(selfake);
        x3=fakesamples.Cmax(selfake);
        y=base+fakesamples.elv_above_base(selfake);
        miny=min(samples.base_level);
        maxy=max(1,max(samples.site_elv)*1.2);
        %     maxx=max(max(samples.apparent_years),apparent_age(max(x2)))*10;
        maxx=1e6;
        %     for n=1:numel(x2)
        %         % horizontla lines
        %         plot(max(minx,apparent_age([x1(n),x3(n)])),y(n)*[1,1],'-b')
        %     end
        for plotline=linspace(0,1,100)
            % vertical lines
            xplotline=apparent_age(x1).*plotline+apparent_age(x2).*(1-plotline);
            plot(xplotline,y,'-b')
            xplotline=apparent_age(x3).*plotline+apparent_age(x2).*(1-plotline);
            plot(xplotline,y,'-b')
        end

        plot(apparent_age(x2),y,'-k')
        %     plot(apparent_age(x1),y,':k') % min
        %     plot(apparent_age(x3),y,':k') % max

        count_samples=0;
        selsamp_index=find(selsamp)';
        [~,order_elv]=sort(samples.elv_above_base(selsamp_index),'descend');
        for n=selsamp_index(order_elv)
            count_samples=count_samples+1;
            y=base+samples.elv_above_base(n);
            %         text(apparent_age(samples.C(n)+samples.dC(n)),y,samples.name{n},'Color','k','Clipping', 'on')
            final_x=apparent_age(samples.C(n)+samples.dC(n));
            final_y=y;
            initial_x=(apparent_age(samples.C(n)+samples.dC(n)))/5;
            initial_y=max(base+samples.elv_above_base(selsamp))-...
                (count_samples-1)/(numel(find(selsamp)')-1)*...
                (max(base+samples.elv_above_base(selsamp))-min(base+samples.elv_above_base(selsamp)));
            text(initial_x,initial_y,samples.name{n},'Color','k','Clipping', 'on','HorizontalAlignment','right')
            plot([initial_x,final_x],[initial_y,final_y],'-k')
            plot(max(minx,apparent_age(samples.C(n)+[-1,1]*samples.dC(n))),y*[1,1],'-r')
            plot(apparent_age(samples.C(n)),y,'or', 'MarkerFaceColor', 'r')
        end
        ylim([miny maxy])
        xlim([minx  maxx])
        % set(gca,'YTickLabel',[]);
        set(gca, 'XScale', 'log')
        %     title([consts.nuclidesstring{consts.nuclides==nuc} ' apparent age'])
        % xlabel('a')
        if nplot==1
            if base==-1e-6 % NODATA value
                ylabel('\Delta Elevation (m)')
            else
                ylabel('Elevation (m)')
            end
        else
            set(gca,'Yticklabel',[])
        end
        xlabel(['(' consts.nuclidesstring{consts.nuclides==nuc} ' a)'])
        box on
        grid on

    end
end

end % end function

