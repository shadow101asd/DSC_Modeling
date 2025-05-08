%% Zhan et al. Sim - Compute Performance

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

%% Replicating Zhan et. al. Circular Orbits

% Minimum hop topology

K1 = [0.53338*AU, 0, 0, 0, 0, 0];
K2 = [1.0*AU, 0, 0, 0, 0, 0];
K3 = [1.0582*AU, 0, 0, 0, 0, 0];

n1 = 5;
n2 = 9;
n3 = 10;

loop1 = Circular_Kns(K1, n1, muSu);
loop2 = Circular_Kns(K2, n2, muSu);
loop3 = Circular_Kns(K3, n3, muSu);
Kns = [loop1, loop2, loop3];

numsats = length(Kns);
for n = 1:numsats
    Kn = Kns(:,n);
    XSats(:,:,n) = propagateFromKeplerians(Kn,muSu,etR);
end

%% Analysis

% Extract Performance Parameters
[graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats);
DequivMean = mean(Dequivs)/AU;
DequivMin = min(Dequivs)/AU;
DequivMax = max(Dequivs)/AU;
DrealMean = mean(Dpaths)/AU;
DrealMin = min(Dpaths)/AU;
DrealMax = max(Dpaths)/AU;

DEaMa = distanceBetweenXs(XEa,XMa);
DEaSats = distanceBetweenXs(XEa,XSats);
DMaSats = distanceBetweenXs(XMa,XSats);
Xs = collateXs(XEa,XMa,XSats);


%% Plotting
%clf(1);clf(2);clf(3);clf(4);clf(5);clf(10);

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
    clf(10)
    viscircles([0,0], 2.2, Color="black");
    set(gca, "XLim", [-1.75, 1.75], "YLim", [-1.75, 1.75], "DataAspectRatio", [1, 1, 1]);
    hold on
    % Plot graphs
    p1 = plot(graphslist{t}, 'XData', reshape(Xs(1,t,:)/AU,numsats+2,1,1), ...
        'YData', reshape(Xs(2,t,:)/AU,numsats+2,1,1), "EdgeColor", [0.3010 0.7450 0.9330]);
    % Highlight Sun + Planets
    scatter(0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black');
    scatter(XEa(1,t)/AU,XEa(2,t)/AU,2e2, 'o', 'blue','filled','MarkerEdgeColor', 'black');
    scatter(XMa(1,t)/AU,XMa(2,t)/AU,2e2, 'o', 'red','filled','MarkerEdgeColor', 'black');
    % Highlight paths
    highlight(p1, paths{t}, 'LineWidth', 6, 'EdgeColor','red');
    hold off
    % Add axis labels and title
    xlabel('x [AU]');
    ylabel('y [AU]');
    title(sprintf(['Graph Plots Over Time, Mars-Earth MBSP. \n ' ...
        'Zhan et al. Optimal Cicular Loops Constellation for DCutoff = 0.7 AU \n' ...
        '%d Satellites: Resulting ED: %f AU'], ...
        numsats, mean(Dequivs)/AU));
    pause(0.1);
end