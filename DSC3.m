%% DSC3: Carthweel FF

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


%% Optimization

[aEai, ~, ~, ~, wEai, ~] = Cartesian2Keplerian(XEa(:,1),muSu);
Ki = [aEai; 0.25; 0.0; 0.0; 1.7412; 0.0];

% Genetic Algorithm

% Constraints!

MAXSATS = 50; % Has to be an integer

A = [1 -1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
b = [0 0 0 0 0];
Aeq = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
beq = [0 0 0 0 0];
lb = [0.0 0.0  1       1       1];
ub = [0.6 0.6  MAXSATS MAXSATS MAXSATS];
intcon = [3,4,5];
options = optimoptions('ga','Display','iter','FunctionTolerance', 1e-3, ...
    'MaxStallGenerations', 7, 'PopulationSize', 1000);

% Running the GA:
% X_opt = [0.15,0.3,2,30,5]; % For testing/plotting

[X_opt, fval, ~, ~] = ga(@(X) wrapperFunc3(XEa,XMa,Ki,X),5,A,b,[],[],lb,ub,@(X)nonlcon(X,MAXSATS),intcon,options);

%% Analysis

% Unpacking Data Output
[es_opt,shells_opt,cartwheels_opt] = unpackVars(X_opt);
numsats_opt = sum(shells_opt)*cartwheels_opt;
emin_opt = min(es_opt);
emax_opt = max(es_opt);
numshells_opt = length(shells_opt);
numbranches_opt = max(shells_opt);

% Recalculate obtained orbits
XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells_opt,es_opt,cartwheels_opt);

% Distances & Effective Distances
 
DEaMa = distanceBetweenXs(XEa,XMa);
DEaSats = distanceBetweenXs(XEa,XSats);
DMaSats = distanceBetweenXs(XMa,XSats);

[graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats);
[graphslistVe,numedgesVe,pathsVe,DpathsVe,DequivsVe] = networkAnalysis(XEa,XVe,"Earth","Venus",XSats);
Xs = collateXs(XEa,XMa,XSats);
Xs2 = collateXs(XEa,XVe,XSats);
[~,nT] = size(XEa);

%% Plotting
%clf(1);clf(2);clf(3);clf(4);clf(5);clf(10);

figure(1)
scatter3(0,0,0,'filled', 'O', 'yellow')
hold on
plotSPICEPositions(XEa, '.', 'blue');
plotSPICEPositions(XMa, '.', 'red');
plotSPICEPositions(XVe, '.', 'green');
for n = 1:numsats_opt
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
for n = 1:numsats_opt
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
for n = 1:numsats_opt
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
    clf(10)
    viscircles([0,0], 2.2, Color="black");
    set(gca, "XLim", [-1.75, 1.75], "YLim", [-1.75, 1.75], "DataAspectRatio", [1, 1, 1]);
    hold on
    % Plot graphs
    p2 = plot(graphslistVe{t}, 'XData', reshape(Xs2(1,t,:)/AU,numsats_opt+2,1,1), ...
        'YData', reshape(Xs2(2,t,:)/AU,numsats_opt+2,1,1), "EdgeColor", [0.3010 0.7450 0.9330]);
    p1 = plot(graphslist{t}, 'XData', reshape(Xs(1,t,:)/AU,numsats_opt+2,1,1), ...
        'YData', reshape(Xs(2,t,:)/AU,numsats_opt+2,1,1), "EdgeColor", [0.3010 0.7450 0.9330]);
    % Highlight Sun + Planets
    scatter(0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black');
    scatter(XEa(1,t)/AU,XEa(2,t)/AU,2e2, 'o', 'blue','filled','MarkerEdgeColor', 'black');
    scatter(XMa(1,t)/AU,XMa(2,t)/AU,2e2, 'o', 'red','filled','MarkerEdgeColor', 'black');
    scatter(XVe(1,t)/AU,XVe(2,t)/AU,2e2, 'o', 'green','filled','MarkerEdgeColor', 'black');
    % Highlight paths
    highlight(p1, paths{t}, 'LineWidth', 6, 'EdgeColor','red');
    highlight(p2, pathsVe{t}, 'LineWidth', 6, 'EdgeColor','green');
    hold off
    % Add axis labels and title
    xlabel('x [AU]');
    ylabel('y [AU]');
    title(sprintf(['Graph Plots Over Time, Mars-Earth MBSP & Venus-Earth MBSP. \n ' ...
        '%d Satellites: %d cartwheel formations, each with %d branches and %d shells. \n ' ...
        'Eccentricity range: [%f - %f]. Resulting ED: %f AU'], ...
        numsats_opt, cartwheels_opt, numbranches_opt, numshells_opt, emin_opt, emax_opt,mean(Dequivs)/AU));
    pause(0.1);
end

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
    