function [ ] = make_climatecurves_ant()

%% Written by Ángel Rodés -- SUERC
% angelrodes@gmail.com
% 2018

% This script generates paleo-climate curves:
%  -Ages
%  -d18O: combined curves from the Antarctic FIVE-CORE average, LR04, Zachos et al. (2008) 
%    and Hansen et al. (2013), scaled to LR04 (Lisiecki & Raymo, 2005)
%  -MIS names and MIS start and end ages

% Version
climatecurves.ver='FIVE-CORE+ v.1.3';

%% load LR04 data
% Lisiecki, L. E., and M. E. Raymo (2005), A Pliocene-
% Pleistocene stack of 57 globally distributed benthic d18O records, 
% Paleoceanography,20, PA1003, doi:10.1029/2004PA001071.
data=csvread('LR04.csv');
LR04.age=1000*data(:,1)';
LR04.d18O=data(:,2)';

%% load and scale FIVE_CORE data
% Buizert, C., Sigl, M., Severi, M. et al. (2018)
% Abrupt ice-age shifts in southern westerly winds and Antarctic climate
% forced from the north. Nature 563, 681–685. 
% https://doi.org/10.1038/s41586-018-0727-5
data=csvread('five_core.csv');
FIVE_CORE.age=data(:,1)';
FIVE_CORE.d8Otemp=data(:,2)';
p = polyfit(FIVE_CORE.d8Otemp,interp1(LR04.age,LR04.d18O,FIVE_CORE.age),1);
FIVE_CORE.d18O=polyval(p,FIVE_CORE.d8Otemp); % scale to LR04

%% Load an scale Zachos et al. (2008) and Hansen et al. (2013) data
% SOURCE: http://www.columbia.edu/~mhs119/Sensitivity+SL+CO2/Table.txt
% NAME OF DATA SET: Cenozoic Global Deep-Sea Stable Isotope Data
% LAST UPDATE: 10/2008 (Original receipt by WDC Paleo)
% CONTRIBUTOR: Jim Zachos, University of California, Santa Cruz
% IGBP PAGES/WDCA CONTRIBUTION SERIES NUMBER: 2008-098
% WDC PALEO CONTRIBUTION SERIES CITATION: 
% Zachos, J., et al. 2008.
% Cenozoic Global Deep-Sea Stable Isotope Data. 
% IGBP PAGES/World Data Center for Paleoclimatology 
% Data Contribution Series # 2008-098. 
% NOAA/NCDC Paleoclimatology Program, Boulder CO, USA.
data=csvread('Zachos2008Hansen2013.csv');
ZACHOS.age=1e6*data(:,1)';
ZACHOS.d8Otemp=data(:,2)';
p = polyfit(interp1(ZACHOS.age,ZACHOS.d8Otemp,LR04.age(LR04.age>min(ZACHOS.age))),LR04.d18O(LR04.age>min(ZACHOS.age)),1);
ZACHOS.d18O=polyval(p,ZACHOS.d8Otemp); % scale to LR04

%% combine curves
climatecurves.age=[FIVE_CORE.age,LR04.age(LR04.age>max(FIVE_CORE.age)),ZACHOS.age(ZACHOS.age>max(LR04.age))];
climatecurves.d18O=[FIVE_CORE.d18O,LR04.d18O(LR04.age>max(FIVE_CORE.age)),ZACHOS.d18O(ZACHOS.age>max(LR04.age))];

%% remove duplicated ages
[~,ia,~] = unique(climatecurves.age,'first');
climatecurves.age=climatecurves.age(ia);
climatecurves.d18O=climatecurves.d18O(ia);

%% simplify curve for cosmonuclide accumulation
% newages=[PALEO.age(PALEO.age<100e3),10.^(5:log10(1.005):log10(max(PALEO.age)))]; % 0.5% increase after 100ka
newages=[0:10:90,...
    100:100:20e3-100,...
    20e3:200:50e3-200,...
    50e3:500:1e5-200,...
    10.^(5:log10(1.01):log10(max(climatecurves.age)))]; % 1% increase after 100ka
newd18O=interp1(climatecurves.age,climatecurves.d18O,newages);
climatecurves.age=newages;
climatecurves.d18O=newd18O;

%% MIS names and ages
climatecurves.MISname=[{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'10'},{'11'},{'12'},{'13'},{'14'},{'15'},{'16'},{'17'},{'18'},{'19'},{'20'},{'21'},{'22'},{'23'},{'24'},{'25'},{'26'},{'27'},{'28'},{'29'},{'30'},{'31'},{'32'},{'33'},{'34'},{'35'},{'36'},{'37'},{'38'},{'39'},{'40'},{'41'},{'42'},{'43'},{'44'},{'45'},{'46'},{'47'},{'48'},{'49'},{'50'},{'51'},{'52'},{'53'},{'54'},{'55'},{'56'},{'57'},{'58'},{'59'},{'60'},{'61'},{'62'},{'63'},{'64'},{'65'},{'66'},{'67'},{'68'},{'69'},{'70'},{'71'},{'72'},{'73'},{'74'},{'75'},{'76'},{'77'},{'78'},{'79'},{'80'},{'81'},{'82'},{'83'},{'84'},{'85'},{'86'},{'87'},{'88'},{'89'},{'90'},{'91'},{'92'},{'93'},{'94'},{'95'},{'96'},{'97'},{'98'},{'99'},{'100'},{'101'},{'102'},{'103'},{'104'},{'G1'},{'G2'},{'G3'},{'G4'},{'G5'},{'G6'},{'G7'},{'G8'},{'G9'},{'G10'},{'G11'},{'G12'},{'G13'},{'G14'},{'G15'},{'G16'},{'G17'},{'G18'},{'G19'},{'G20'},{'G21'},{'G22'},{'K1'},{'K2'},{'KM1'},{'KM2'},{'KM3'},{'KM4'},{'KM5'},{'KM6'},{'M1'},{'M2'},{'MG1'},{'MG2'},{'MG3'},{'MG4'},{'MG5'},{'MG6'},{'MG7'},{'MG8'},{'MG9'},{'MG10'},{'MG11'},{'MG12'},{'Gi1'},{'Gi2'},{'Gi3'},{'Gi4'},{'Gi5'},{'Gi6'},{'Gi7'},{'Gi8'},{'Gi9'},{'Gi10'},{'Gi11'},{'Gi12'},{'Gi13'},{'Gi14'},{'Gi15'},{'Gi16'},{'Gi17'},{'Gi18'},{'Gi19'},{'Gi20'},{'Gi21'},{'Gi22'},{'Gi23'},{'Gi24'},{'Gi25'},{'Gi26'},{'Gi27'},{'Gi28'},{'Co1'},{'Co2'},{'Co3'},{'Co4'},{'CN1'},{'CN2'},{'CN3'},{'CN4'},{'CN5'},{'CN6'},{'CN7'},{'CN8'},{'N1'},{'N2'},{'N3'},{'N4'},{'N5'},{'N6'},{'N7'},{'N8'},{'N9'},{'N10'},{'NS1'},{'NS2'},{'NS3'},{'NS4'},{'NS5'},{'NS6'},{'Si1'},{'Si2'},{'Si3'},{'Si4'},{'Si5'},{'Si6'},{'ST1'},{'ST2'},{'ST3'},{'ST4'},{'T1'},{'T2'},{'T3'},{'T4'},{'T5'},{'T6'},{'T7'},{'T8'},{'TG1'},{'TG2'},{'TG3'},{'TG4'},{'TG5'},{'TG6'}];
climatecurves.MISend=1e3*[0,14,29,57,71,130,191,243,300,337,374,424,478,533,563,621,676,712,761,790,814,866,900,917,936,959,970,982,1014,1031,1062,1081,1104,1114,1141,1190,1215,1244,1264,1286,1304,1320,1344,1362,1383,1405,1424,1452,1469,1492,1510,1530,1547.5,1570,1585,1608,1628.5,1642.5,1670,1697.5,1715,1743,1758,1782,1802.5,1816,1826,1832.5,1849,1859.5,1875,1898,1915,1941,1965,1990,2017,2043,2088,2103,2125,2146,2168,2192,2207.5,2236,2250,2273,2291,2309,2333,2350,2373,2387,2407,2427,2452,2477,2494,2510,2540,2554,2575,2595,2614,2638,2652,2681,2690,2704,2730,2759,2777,2798,2820,2838,2858,2876,2893,2913,2937,2966,2982.5,2999,3025,3039,3055,3087,3097,3119,3150,3167,3184,3212,3238,3264,3312,3332,3347,3372,3387,3444,3471,3517,3532,3546,3566,3578,3592,3619,3637,3660,3676,3705,3719,3742,3752,3768,3798,3822,3835,3862,3879,3912,3923,3939,3952,3978,4007,4029,4048,4085,4098,4146,4175,4192,4211,4232,4259,4286,4303,4327,4335,4356,4371,4395,4420,4446,4457,4487,4508,4523,4538,4570,4587,4603,4622,4648,4658,4684,4702.5,4722.5,4737,4766,4778,4807,4821,4840,4860,4883,4904,4931,4952.5,4976,4985,5002,5016,5038,5070,5094,5116,5165,5188,5241,5266,5289,5301,5315];
climatecurves.MISstart=1e3*[14,29,57,71,130,191,243,300,337,374,424,478,533,563,621,676,712,761,790,814,866,900,917,936,959,970,982,1014,1031,1062,1081,1104,1114,1141,1190,1215,1244,1264,1286,1304,1320,1344,1362,1383,1405,1424,1452,1469,1492,1510,1530,1547.5,1570,1585,1608,1628.5,1642.5,1670,1697.5,1715,1743,1758,1782,1802.5,1816,1826,1832.5,1849,1859.5,1875,1898,1915,1941,1965,1990,2017,2043,2088,2103,2125,2146,2168,2192,2207.5,2236,2250,2273,2291,2309,2333,2350,2373,2387,2407,2427,2452,2477,2494,2510,2540,2554,2575,2595,2614,2638,2652,2681,2690,2704,2730,2759,2777,2798,2820,2838,2858,2876,2893,2913,2937,2966,2982.5,2999,3025,3039,3055,3087,3097,3119,3150,3167,3184,3212,3238,3264,3312,3332,3347,3372,3387,3444,3471,3517,3532,3546,3566,3578,3592,3619,3637,3660,3676,3705,3719,3742,3752,3768,3798,3822,3835,3862,3879,3912,3923,3939,3952,3978,4007,4029,4048,4085,4098,4146,4175,4192,4211,4232,4259,4286,4303,4327,4335,4356,4371,4395,4420,4446,4457,4487,4508,4523,4538,4570,4587,4603,4622,4648,4658,4684,4702.5,4722.5,4737,4766,4778,4807,4821,4840,4860,4883,4904,4931,4952.5,4976,4985,5002,5016,5038,5070,5094,5116,5165,5188,5241,5266,5289,5301,5315,5320];

% %% generate a paleotemperature curve (not used)
% temperatures=[-6,0]; % references for maximum d18O <115 ka and mean dO18 <100 a
% % see https://en.wikipedia.org/wiki/File:All_palaeotemps.png#Summary
% deltas=[max(PALEO.d18O(PALEO.age<115000)),mean(PALEO.d18O(PALEO.age<100))];
% p = polyfit(deltas,temperatures,1);
% PALEO.temp=polyval(p,PALEO.d18O);

% save climatecurves climatecurves
save('climatecurves.mat','climatecurves','-v7')

%% plot curves
figure('units','normalized','outerposition',[0 0.5 0.5 0.5],'Name','Climate proxy')
subplot(2,1,1)
hold on
% plot(climatecurves.age,climatecurves.d18O,'-','Color',[0.7 0.7 0.7],'LineWidth',3)
plot(ZACHOS.age,ZACHOS.d18O,'-g','LineWidth',3)
plot(LR04.age,LR04.d18O,'-b','LineWidth',3)
plot(FIVE_CORE.age,FIVE_CORE.d18O,'-r','LineWidth',3)
plot(climatecurves.age,climatecurves.d18O,'-','Color','k','LineWidth',1)
legend('Zachos+Hansen','LR04','FIVE-CORE',climatecurves.ver,...
    'Location','southwest')
ylabel('\delta^{18}O')
set(gca, 'Xdir', 'reverse')
xlabel('Age (a)')
set(gca, 'XScale', 'log')
xlim([1e2 1e8])
box on 
grid on
title('Scaled records')

subplot(2,1,2)
hold on
plot(climatecurves.age/1e3,climatecurves.d18O,'-','Color','k','LineWidth',1)
for n=1:numel(climatecurves.MISend)
    text(climatecurves.MISend(n)/1e3,max(climatecurves.d18O),...
        '|',...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Clipping', 'on','Color','b')
    text((climatecurves.MISend(n)/2+climatecurves.MISstart(n)/2)/1e3,max(climatecurves.d18O),...
        climatecurves.MISname(n),...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Clipping', 'on','Color','b')
end
ylabel('\delta^{18}O')
xlabel('Age (ka)')
set(gca, 'Xdir', 'reverse')
% set(gca, 'XScale', 'log')
xlim([0 1.2e5/1e3])
ylim([min(climatecurves.d18O(climatecurves.age<1.2e5)) max(climatecurves.d18O)*1.1])
set (gca, 'Clipping', 'on');
box on
grid on
title(climatecurves.ver)

%figure(2000)
%subplot(2,1,2)
%hold on
%% plot(climatecurves.age,climatecurves.d18O,'-','Color',[0.7 0.7 0.7],'LineWidth',3)
%plot(ZACHOS.age,ZACHOS.d18O,'-g','LineWidth',2)
%plot(LR04.age,LR04.d18O,'-b','LineWidth',2)
%plot(FIVE_CORE.age,FIVE_CORE.d18O,'-r','LineWidth',2)
%plot(climatecurves.age,climatecurves.d18O,'-','Color','k','LineWidth',1)
%legend('Zachos+Hansen','LR04','FIVE-CORE',climatecurves.ver,...
%    'Location','southeast')
%ylabel('\delta^{18}O')
%set(gca, 'Xdir', 'reverse')
%xlabel('Age (a)')
%set(gca, 'XScale', 'log')
%xlim([1e2 1e8])
%box on 
%grid on
%title('B')

end
