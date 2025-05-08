%% DSC4: Loop over max_satellites of DSC3

clear all

%% Import Kernels and Constants from SPICE

path_to_generic_kernels = '/Users/jpenot/Documents/MATLAB/SPICE/naif.jpl.nasa.gov/pub/naif/generic_kernels';
path_to_mice            = '/Users/jpenot/Documents/MATLAB/SPICE/mice';
addpath(strcat(path_to_mice,'/src/mice'))
addpath(strcat(path_to_mice,'/lib'))

% Load the datafiles (kernels)
% (1) leap-seconds
cspice_furnsh([path_to_generic_kernels,'/lsk/naif0012.tls.pc']);
% (2) planets
cspice_furnsh([path_to_generic_kernels,'/spk/planets/de430.bsp']);
% (3) gravity constants
cspice_furnsh([path_to_generic_kernels,'/pck/gm_de431.tpc']);
% (4) planetary constant - you can also open this file with a text editor
% to read its contentmex -setup
cspice_furnsh([path_to_generic_kernels,'/pck/pck00010.tpc']);

%% Loading Data

global muSu etR AU
AU = 1.496e8; % 1 AU in km
muSu = cspice_bodvcd(10, 'GM', 10); % GM of the Sun with 10 significant digits
interval = 3600*24*5; %5days
ref = 'ECLIPJ2000';

% Dates
date0 = '2035 Jan 1 12:00:00 UTC'; % Simulation start
date1 = '2040 Jan 1 12:00:00 UTC'; % Simulation ends
et0 = cspice_str2et(date0);
et1 = cspice_str2et(date1);
etR = et0:interval:et1; % Row vector of times between start and end date

% Ephemeris Data (km, km/s)
XEa = cspice_spkezr('3', etR, ref, 'NONE', '10'); % Earth ephemeris wrt the Sun in ECLIPJ2000
XMa = cspice_spkezr('4', etR, ref, 'NONE', '10'); % Mars ephemeris wrt the Sun in ECLIPJ2000
XVe = cspice_spkezr('2', etR, ref, 'NONE', '10'); % Venus ephemeris wrt the Sun in ECLIPJ2000

% Num timesteps
[~,nT] = size(XEa);

%% Optimization

[aEai, ~, ~, ~, wEai, ~] = Cartesian2Keplerian(XEa(:,1),muSu);
Ki = [aEai; 0.25; 0.0; 0.0; 1.7412; 0.0];

% Genetic Algorithm

% Constraints!

MAX_MAXSATS = 500; % Has to be an integer
MIN_MAXSATS = 500; % Has to be an integer
SPACING = 1; % Has to be an integer
SATNUMS = MIN_MAXSATS:SPACING:MAX_MAXSATS;
[~,N] = size(SATNUMS);

A = [1 -1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
b = [0 0 0 0 0];
Aeq = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
beq = [0 0 0 0 0];
intcon = [3,4,5];

% Initialize data fields
numsats_optN(N) = 0;
emin_optN(N) = 0;
emax_optN(N) = 0;
numshells_optN(N) = 0;
numbranches_optN(N) = 0;
numcarts_optN(N) = 0;

%% Running

i = 1;
for MAXSATS = SATNUMS
    MAXSATS % For progress tracking
    lb = [0.0 0.0  1       1       1];
    ub = [0.6 0.6  MAXSATS MAXSATS MAXSATS];

    options = optimoptions('ga','Display','iter','FunctionTolerance', 1e-3, ...
            'MaxStallGenerations', 7, 'PopulationSize', min(100*MAXSATS,1000));
    
    % Running the GA:
    
    [X_opt, fval, ~, ~] = ga(@(X) wrapperFunc3(XEa,XMa,Ki,X),5,A,b,[],[],lb,ub,@(X)nonlcon(X,MAXSATS),intcon,options);
    
    % Analysis
    
    % Unpacking Data Output
    [es_opt,shells_opt,cartwheels_opt] = unpackVars(X_opt) % Print for progress tracking
    numsats_optN(i) = sum(shells_opt)*cartwheels_opt;
    emin_optN(i) = min(es_opt);
    emax_optN(i) = max(es_opt);
    numshells_optN(i) = length(shells_opt);
    numbranches_optN(i) = max(shells_opt);
    numcarts_optN(i) = cartwheels_opt;

    i = i+1;
end

%% Analysis

DequivMean(N) = 0.0;
DequivMin(N) = 0.0;
DequivMax(N) = 0.0;
DrealMean(N) = 0.0;
DrealMin(N) = 0.0;
DrealMax(N) = 0.0;

for i = 1:N
    % Recalculate obtained orbits
    es_opt = linspace(emin_optN(i),emax_optN(i),numshells_optN(i));
    shells_opt = numbranches_optN(i)*ones(numshells_optN(i),1);
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells_opt,es_opt,numcarts_optN(i));
    
    % Extract Performance Parameters
    [graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats);
    DequivMean(i) = mean(Dequivs)/AU;
    DequivMin(i) = min(Dequivs)/AU;
    DequivMax(i) = max(Dequivs)/AU;
    DrealMean(i) = mean(Dpaths)/AU;
    DrealMin(i) = min(Dpaths)/AU;
    DrealMax(i) = max(Dpaths)/AU;
end

DEaMa = distanceBetweenXs(XEa,XMa);


%% Plotting

c = 3*10^5; % Speed of light, km/s
fitequiv = fit(SATNUMS',DequivMean', 'power1');

figure(1)
plot(SATNUMS,DequivMean, 'Color', 'blue');
hold on
plot(fitequiv);
plot(SATNUMS,DequivMax, 'Color', 'blue', 'LineStyle', '--');
plot(SATNUMS,DequivMin, 'Color', 'blue', 'LineStyle', '--');
hold off
legend('Mean', 'Power Function Fit', 'Bounds at each # of satellites')
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Distance [AU]')
title('Equivalent Transmission Distance')

figure(2)
plot(SATNUMS,DrealMean*AU/c/60, 'Color', 'red');
hold on
plot(SATNUMS,DrealMax*AU/c/60, 'Color', 'red', 'LineStyle', '--');
yline(mean(DEaMa)/c/60, 'green');
yline(max(DEaMa)/c/60, 'Color', 'green', 'LineStyle', '--');
yline(min(DEaMa)/c/60, 'Color', 'green', 'LineStyle', '--');
plot(SATNUMS,DrealMin*AU/c/60, 'Color', 'red', 'LineStyle', '--');
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Transmission Time [min]')
title('Real Transmission Time')
legend('Constellation Mean', 'Bounds at each # of satellites', ...
    'Earth-Mars direct transmission time (mean)', 'Earth-Mars direct transmission time (max/min throughout orbits)')

figure(3)
plot(SATNUMS, numshells_optN);
hold on
plot(SATNUMS, numbranches_optN);
plot(SATNUMS, numcarts_optN);
plot(SATNUMS, numsats_optN);
fplot(@(x) x, [0 max(SATNUMS)]);
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('#')
title('Cartwheel Constellation Configuration Numbers as a Function of Number of Satellites')
legend('Number of Shells per CFFC', 'Number of Branches per CFFC', 'Number of CFFCs', 'Number of satellites in optimal constellation', 'y=x')

figure(4)
plot(SATNUMS, emin_optN);
hold on
plot(SATNUMS, emax_optN);
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Eccentricity')
title('Cartwheel Constellation Eccentricities as a Function of Number of Satellites')
legend('emin', 'emax')


%% Saving Data

filename = "DSC4_Data/DSC4_min"+ num2str(MIN_MAXSATS) + "max" + num2str(MAX_MAXSATS) + "spacing" + num2str(SPACING);
save(filename)

%% Functions

function Out = wrapperFunc3(X1,X2,Ki,X)
    global etR muSu
    
    [es,shells,cartwheels] = unpackVars(X);
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells,es,cartwheels);

    Out = bestLinkBudget(X1,X2,XSats);
end

function [es,shells,cartwheels] = unpackVars(X)
    % Unpack variables
    emin = X(1);
    emax = X(2);
    numshells = X(3);
    numbranches = X(4);
    cartwheels = X(5);
    es = linspace(emin,emax,numshells);
    shells = numbranches*ones(numshells,1);
end

function metric = bestLinkBudget(X1,X2,XSats)
    global AU
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

function [c,ceq] = nonlcon(X,MAXSATS)
    [~,shells,cartwheels] = unpackVars(X);
    numsats = sum(shells)*cartwheels;
    ceq = []; % No equality constraints for MINLPs
    c = numsats - MAXSATS;
end
    