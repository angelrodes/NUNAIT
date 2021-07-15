function [  ] = plot_samples(selectedfile)
% Plot site elevations

%% load data
load(selectedfile)

figure('units','normalized','outerposition',[0.5 0 0.5 0.5],'Name','Profile')
subplot(1,2,1)
hold on
title(samples.filename)

x=samples.apparent_years;
dx=samples.dapparent_years;
y=samples.site_elv;

for n=1:numel(x)
    plot(x(n)+[-1,1]*dx(n),y(n)*[1,1],'-r')
    plot(x(n),y(n),'or')
    text(x(n)+dx(n),y(n),[samples.name{n} ' (' num2str(samples.isotope(n)) ')'])
end

ylabel('Elevation (m)')
xlabel('Age (a)')

subplot(1,2,2)
hold on

x=samples.apparent_years;
dx=samples.dapparent_years;
y=samples.elv_above_base;

for n=1:numel(x)
    plot(x(n)+[-1,1]*dx(n),y(n)*[1,1],'-r')
    plot(x(n),y(n),'or')
    text(x(n)+dx(n),y(n),[samples.name{n} ' (' num2str(samples.isotope(n)) ')'])
end

ylabel('\Delta Elevation (m)')
xlabel('Age (a)')


end

