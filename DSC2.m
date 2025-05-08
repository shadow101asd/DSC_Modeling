%% DSC2

clear all

%% Import Kernels and Constants from SPICE

path_to_generic_kernels = '/Users/jpenot/Documents/MATLAB/SPICE/naif.jpl.nasa.gov/pub/naif/generic_kernels';
path_to_mice            = '/Users/jpenot/Documents/MATLAB/SPICE/mice';
addpath(strcat(path_to_mice,'/src/mice'))
addpath(strcat(path_to_mice,'/lib'))

% Load the datafiles (kernels)
% (1) leap-seconds
cspice_furnsh( [path_to_generic_kernels,'/lsk/naif0012.tls.pc']);
% (2) planets
cspice_furnsh( [path_to_generic_kernels,'/spk/planets/de430.bsp']);
% (3) gravity constants
cspice_furnsh( [path_to_generic_kernels,'/pck/gm_de431.tpc']);
% (4) planetary constant - you can also open this file with a text editor
% to read its contentmex -setup
cspice_furnsh( [path_to_generic_kernels,'/pck/pck00010.tpc']);

%% Loading Data

global muSu etR AU
AU = 1.496e8; % 1 AU in km
muSu = cspice_bodvcd(10, 'GM', 10); % GM of the Sun with 10 significant digits
interval = 3600*24*1; %24h
ref = 'ECLIPJ2000';
MAX_ECCENTRICITY = 0.8;
FIXED_ECCENTRICITY = 0.5;
boolFixedE = 0;

% Shells
numshells = 8;
numbranches = 4;
symmetrical = 1; % Bool
emin = 0.1;
emax = 0.7;
shells = numbranches*ones(numshells,1);
es = linspace(emin,emax,numshells);
numsats = sum(shells)*(1+symmetrical);
%numsats = 15; % Has to be an integer

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


%% Optimization

[aEai, ~, ~, ~, wEai, ~] = Cartesian2Keplerian(XEa(:,1),muSu);
K0 = [aEai; 0.25; 0.0; 0.0; 1.7412; 0.0]; % Initial guess


% Optimization conditions
% Eccentricity must be positive
A = [0 0 0 0 0 0; 0 -1 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
b = [0;0;MAX_ECCENTRICITY;0;0;0];

% First attempt - keep satellites in J2000 ecliptic plane
Aeq = [0 0 0 0 0 0; 0 boolFixedE 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
beq = [0;FIXED_ECCENTRICITY;0;0;0;0];

KiOpt = K0; % TESTING

[~, nT] = size(etR);
% OPTIONS = optimoptions('fmincon', 'Display','iter', ...
%                                 'StepTolerance', 1e-3, ...
%                                 "Algorithm","interior-point",...
%                                 "EnableFeasibilityMode",true,...
%                                 "SubproblemAlgorithm","cg");
% [KiOpt, fval] = fmincon(@(Ki) wrapperFunc(XEa,XMa,Ki,numsats,shells,es,symmetrical), K0, A, b, [], [], [], [] ,[], OPTIONS);

%% Analysis

% Recalculate obtained orbits
XSats = NSATSpropagateFromKepleriansSHELLS(KiOpt,muSu,etR,shells,es,symmetrical);

% Distances & Effective Distances
 
DEaMa = distanceBetweenXs(XEa,XMa);
DEaSats = distanceBetweenXs(XEa,XSats);
DMaSats = distanceBetweenXs(XMa,XSats);

[graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats);
[graphslistVe,numedgesVe,pathsVe,DpathsVe,DequivsVe] = networkAnalysis(XEa,XVe,"Earth","Venus",XSats);
Xs = collateXs(XEa,XMa,XSats);
Xs2 = collateXs(XEa,XVe,XSats);

%% Plotting
clf(1);clf(2);clf(3);clf(4);clf(5);clf(10);

figure(1)
scatter3(0,0,0,'filled', 'O', 'yellow')
hold on
plotSPICEPositions(XEa, '.', 'blue');
plotSPICEPositions(XMa, '.', 'red');
plotSPICEPositions(XVe, '.', 'green');
for n = 1:numsats
    plotSPICEPositions(XSats(:,:,n), '.', 'black');
end
hold off
set(gca, 'DataAspectRatio', [1 1 1]) % Fix axes' aspect ratio
legend('Sun', 'Earth', 'Mars', 'Venus', 'Satellites');
title('Full orbits');


figure(2)
plot(etR,DEaMa/AU, 'LineWidth', 4);
title('Distances to Earth [AU]')
LEGEND(1) = "Earth-Mars distance";
hold on
for n = 1:numsats
    plot(etR,DEaSats(:,:,n)/AU);
    LEGEND(n+1) = strcat("Earth-Satellite ", num2str(n), " distance");
end
hold off
legend(LEGEND);
xlabel('t [s]')
ylabel('Distances [AU]')

figure(3)
plot(etR,DEaMa/AU, 'LineWidth', 4);
title('Distances to Mars [AU]')
LEGEND(1) = "Earth-Mars distance";
hold on
for n = 1:numsats
    plot(etR,DMaSats(:,:,n)/AU);
    LEGEND(n+1) = strcat("Mars-Satellite ", num2str(n), " distance");
end
hold off
legend(LEGEND);
xlabel('t [s]')
ylabel('Distances [AU]')

figure(4)
plot(etR,numedges)
xlabel('t [s]');
ylabel('Number of edges in graph at timestep t');
title('Number of internode connections of shorter distance than Earth<->Mars, over time')

figure(5)
plot(etR, DEaMa/AU)
hold on
plot(etR, Dequivs/AU)
plot(etR, Dpaths/AU)
hold off
legend('Earth-Mars Distance', 'Equivalent Transmission Distance', 'Total Actual Transmission Distance');
xlabel('time t [s]');
ylabel('Distance [AU]');

figure(10)
for t=1:nT
    scatter(0,0,'filled', 'O', 'yellow')
    viscircles([0,0], 2.2*AU, Color="black")
    set(gca, "XLim", [-2.5e8, 2.5e8], "YLim", [-2.5e8, 2.5e8]);
    hold on
    p = plot(graphslist{t}, 'XData', reshape(Xs(1,t,:),numsats+2,1,1), ...
        'YData', reshape(Xs(2,t,:),numsats+2,1,1), "EdgeColor", "blue");
    set(gca, 'DataAspectRatio', [1 1 1]) % Fix axes' aspect ratio
    p2 = plot(graphslistVe{t}, 'XData', reshape(Xs2(1,t,:),numsats+2,1,1), ...
        'YData', reshape(Xs2(2,t,:),numsats+2,1,1), "EdgeColor", "blue");
    % Highlight paths
    highlight(p, paths{t}, 'LineWidth', 4, 'EdgeColor','red');
    highlight(p2, pathsVe{t}, 'LineWidth', 4, 'EdgeColor','green');
    pause(0.1);
    hold off
end

%% Functions

function Out = wrapperFunc(X1,X2,Ki,numsats,shells,es,symmetrical)
    global etR muSu
    if isempty(shells) || isempty(es) || isempty(symmetrical)
        XSats = NSATSpropagateFromKeplerians(Ki,muSu,etR,numsats);
    else
        XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells,es,symmetrical);
    end

    Out = bestLinkBudget(X1,X2,XSats);
end

function metric = bestLinkBudget(X1,X2,XSats)
    global AU
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end
    