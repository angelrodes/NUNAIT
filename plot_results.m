function [  ] = plot_results(filename)
% This function plots results generated by fit_nuna_model
%% Angel Rodes, 2020


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

%% Plot full probabilities (variables)
figure('units','normalized','outerposition',[0 0 0.5 1],'Name','Parameters')

subplot(5,1,1)
hold on
x=ic;
res=results.ic;
plot(x,PROBS,'.k')
title('Ice-thinning since maximum')
xlabel(string_ic)
ylabel('P(\chi^2)')
% xlim([max(min(x),res(2)-range(res)*3) min(max(x),res(2)+range(res)*3)])
if range(x)>0
    xlim([min(x) max(x)])
end
% ylim([0 max(PROBS(:))*1.2])
set(gca,'YTickLabel',[]);
box on
grid on

subplot(5,1,2)
hold on
x=W*1e4;
res=results.W*1e4;
plot(x,PROBS,'.k')
title('Weathering')
xlabel(string_W)
ylabel('P(\chi^2)')
% xlim([max(min(x),res(2)-range(res)*3) min(max(x),res(2)+range(res)*3)])
if range(x)>0
    xlim([min(x) max(x)])
end
set(gca, 'XScale', 'log')
% ylim([0 max(PROBS(:))*1.2])
set(gca,'YTickLabel',[]);
box on
grid on


subplot(5,1,3)
hold on
x=E*1e4/1e3;
res=results.E*1e4/1e3;
plot(x,PROBS,'.k')
title('Glacial erosion')
xlabel(string_E)
ylabel('P(\chi^2)')
% xlim([max(min(x),res(2)-range(res)*3) min(max(x),res(2)+range(res)*3)])
if range(x)>0
    xlim([min(x) max(x)])
end
set(gca, 'XScale', 'log')
% ylim([0 max(PROBS(:))*1.2])
set(gca,'YTickLabel',[]);
box on
grid on

subplot(5,1,4)
hold on
x=UR;
res=results.UR;
plot(x,PROBS,'.k')
title('Uplift rate')
xlabel(string_UR)
ylabel('P(\chi^2)')
% xlim([max(min(x),res(2)-range(res)*3) min(max(x),res(2)+range(res)*3)])
if range(x)>0
    xlim([min(x) max(x)])
end
% set(gca, 'XScale', 'log')
% ylim([0 max(PROBS(:))*1.2])
set(gca,'YTickLabel',[]);
box on
grid on

subplot(5,1,5)
hold on
x=DELTAH;
res=results.DELTAH;
plot(x,PROBS,'.k')
title('Z dev.')
xlabel(string_DELTAH)
ylabel('P(\chi^2)')
% xlim([max(min(x),res(2)-range(res)*3) min(max(x),res(2)+range(res)*3)])
if range(x)>0
    xlim([min(x) max(x)])
end
% set(gca, 'XScale', 'log')
% ylim([0 max(PROBS(:))*1.2])
set(gca,'YTickLabel',[]);
box on
grid on


%% Plot glacial model
figure('units','normalized','outerposition',[1 0.5 0.5 0.5],'Name','Glacial model')
% subplot(3,3,[1 3])
hold on
agescaling=1;
x=climatecurves.age/agescaling;
% rerpesent the elevation respect the current average ice surface in the area
base=mean(unique(samples.base_level)); 
y2=base+results.DELTAH(2)+(climatecurves.d18O-climatecurves.d18O(1))/results.D18z(2); % best
% find min and max curves
y1=y2; % min curve
y3=y2; % max curve
% h = waitbar(0,'Calculating glacial surface...');
disp('Calculating glacial surface...')
tic
for n=find(onesigma)'
%    if mod(n,10)==0
%        waitbar(n/sum(onesigma),h);
%    end
    yi=base+DELTAH(n)+(climatecurves.d18O-climatecurves.d18O(1))/D18z(n);
    y1=min(y1,yi);
    y3=max(y3,yi);
end
toc
%close(h)

% plot uncert
for beta=linspace(0,1,100)
    plot(x,y1*beta+y3*(1-beta),'-b')
end
% plot samples 
for n=1:numel(samples.site_elv)
    ys=samples.site_elv(n);
    plot([x(end),0],[ys-results.UR(2)*x(end)*agescaling/1e6,ys],'-r')
    plot([x(end),0],[ys-results.UR(1)*x(end)*agescaling/1e6,ys],':r')
    plot([x(end),0],[ys-results.UR(3)*x(end)*agescaling/1e6,ys],':r')
    plot(0,ys,'or', 'MarkerFaceColor', 'r') % samples
    
end
% plot ice elevation
plot(x,y2,'-k')
% beautify
set(gca, 'Xdir', 'reverse')
set(gca,'YAxisLocation','right')
miny=min(min(y1(climatecurves.age<100e3)),min(samples.site_elv));
maxy=max(max(y3),max(samples.site_elv))+100;
xlim([0 max(30e3,max(samples.C))*1.2]/agescaling)
ylim([miny maxy])
ylabel('Elevation (m)')
xlabel('Age (a)')
title(['Elevation of the ice-surface (' climatecurves.ver ')'])
box on
grid on

%% Plot apparent_age(concentration) profiles
uniquenuclides=unique(samples.isotope)';

figure('units','normalized','outerposition',[0.5 0 0.5 0.5],'Name','Profile models')
nplot=0;
for nuc=uniquenuclides
    nplot=nplot+1;
    subplot(1,numel(uniquenuclides),nplot)
    hold on
    selfake=(fakesamples.isotope==nuc);
    selsamp=(samples.isotope==nuc);
    lambda=mean(samples.l(selsamp));
    apparent_age=@(concentration) log(1-min(max(concentration,0)*lambda,1))/(-lambda); % a
    base=-1e-6; % set base to NODATA value
    minx=1; % plot minimum 1 year
    %     minx=max(1,min(samples.apparent_years)/2);
    x1=fakesamples.Cmin(selfake);
    x2=fakesamples.C(selfake);
    x3=fakesamples.Cmax(selfake);
    y=base+fakesamples.elv_above_base(selfake);
    miny=min(y);
    maxy=max(max(y),max(samples.elv_above_base)*1.2);
    maxx=max(max(samples.apparent_years),apparent_age(max(x2)))*10;
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

    for n=find(selsamp)'
        y=base+samples.elv_above_base(n);
        text(apparent_age(samples.C(n)+samples.dC(n)),y,samples.name{n},'Color','k','Clipping', 'on')
        plot(max(minx,apparent_age(samples.C(n)+[-1,1]*samples.dC(n))),y*[1,1],'-r')
        plot(apparent_age(samples.C(n)),y,'or', 'MarkerFaceColor', 'r')
    end
    ylim([miny maxy])
    xlim([minx  maxx])
    % set(gca,'YTickLabel',[]);
    set(gca, 'XScale', 'log')
    title([consts.nuclidesstring{consts.nuclides==nuc} ' apparent age'])
    % xlabel('a')
    if base==-1e-6 % NODATA value
        ylabel('\Delta Elevation (m)')
    else
        ylabel('Elevation (m)')
    end
    xlabel('(a)')
    box on
    grid on
    
end

end % end function

