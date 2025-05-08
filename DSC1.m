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

global muSu etR OPTSCALING
AU = 1.496e8; % 1 AU in km
muSu = cspice_bodvcd(10, 'GM', 10); % GM of the Sun with 10 significant digits
interval = 3600*24; %24h
ref = 'ECLIPJ2000';
ANIMATION_PAUSE = 0.02; % s
OPTSCALING = 1e-16;

% Dates
date0 = '2035 Jan 1 12:00:00 UTC'; % Simulation start
date1 = '2040 Jan 1 12:00:00 UTC'; % Simulation ends
et0 = cspice_str2et(date0);
et1 = cspice_str2et(date1);
etR = et0:interval:et1; % Row vector of times between start and end date

% Ephemeris Data (km, km/s)
XEa = cspice_spkezr('3', etR, ref, 'NONE', '10'); % Earth ephemeris wrt the Sun in ECLIPJ2000
XMa = cspice_spkezr('4', etR, ref, 'NONE', '10'); % Mars ephemeris wrt the Sun in ECLIPJ2000



%% Simulation

y0 = [-2.0e8; 0.0; 0.0; 0.0; -20.0; 0.0]; % Initial guess

% Optimization conditions
% First attempt - keep satellites in J2000 ecliptic plane
A = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1];
b = [0;0;0;0;0;0];

[~, nT] = size(etR);
OPTIONS = optimoptions('fmincon', 'Display','iter', 'StepTolerance', 1e-10);
[yiOpt, fval] = fmincon(@(yi) wrapperFunc(XEa,XMa,yi), y0, A, b, [], [], [], [] ,[], OPTIONS);
%[yiOpt, fval] = fminsearch(@(yi) wrapperFunc(XEa,XMa,yi), y0);
[~,XSat] = ode113(@(t,y) twobodyode(t,y,muSu), etR, yiOpt);
fval = fval/OPTSCALING;
XSat1 = XSat';
XSat2 = -XSat1;

% Distances

DEaMa = sqrt((XEa(1,:)-XMa(1,:)).^2 + (XEa(2,:)-XMa(2,:)).^2 + (XEa(3,:)-XMa(3,:)).^2);
DEaSat1 = sqrt((XEa(1,:)-XSat1(1,:)).^2 + (XEa(2,:)-XSat1(2,:)).^2 + (XEa(3,:)-XSat1(3,:)).^2);
DEaSat2 = sqrt((XEa(1,:)-XSat2(1,:)).^2 + (XEa(2,:)-XSat2(2,:)).^2 + (XEa(3,:)-XSat2(3,:)).^2);
DMaSat1 = sqrt((XMa(1,:)-XSat1(1,:)).^2 + (XMa(2,:)-XSat1(2,:)).^2 + (XMa(3,:)-XSat1(3,:)).^2);
DMaSat2 = sqrt((XMa(1,:)-XSat2(1,:)).^2 + (XMa(2,:)-XSat2(2,:)).^2 + (XMa(3,:)-XSat2(3,:)).^2);


%% Plotting
clf

figure(1)
scatter3(0,0,0,'filled', 'O', 'yellow')
hold on
plotSPICEPositions(XEa, '.', 'blue');
plotSPICEPositions(XMa, '.', 'red');
plotSPICEPositions(XSat1, '.', 'black');
plotSPICEPositions(XSat2, '.', 'green');
hold off
set(gca, 'DataAspectRatio', [1 1 1]) % Fix axes' aspect ratio
legend('Sun', 'Earth', 'Mars', 'Satellite 1', 'Satellite 2');
title('Full orbits');

% figure(2)
% scatter3(0,0,0,'filled', 'O', 'yellow')
% title('Animation');
% ax = gca;
% ax.XLimMode = 'manual'; ax.YLimMode = 'manual'; ax.ZLimMode = 'manual';
% ax.XLim = 1e8*[-2.5 2.5]; ax.YLim = 1e8*[-2.5 2.5]; ax.ZLim = 1e7*[-1.0 1.0];
% set(gca, 'DataAspectRatio', [1 1 1]); % Fix axes limits
% hold on
% plotSPICEPositions(XEa(:,1), '.', 'blue');
% plotSPICEPositions(XMa(:,1), '.', 'red');
% plotSPICEPositions(XSat(:,1), '.', 'black');
% legend('Sun', 'Earth', 'Mars', 'Satellite','AutoUpdate','off');
% for i = 2:4:nT
%     plotSPICEPositions(XEa(:,i), '.', 'blue');
%     plotSPICEPositions(XMa(:,i), '.', 'red');
%     plotSPICEPositions(XSat(:,i), '.', 'black');
%     pause(ANIMATION_PAUSE)
% end  
% hold off

figure(2)
title('Distances [AU]')
plot(etR,DEaMa/AU);
hold on
plot(etR,max(DEaSat1,DMaSat1)/AU);
plot(etR,max(DEaSat2,DMaSat2)/AU);
hold off
xlabel('t [s]')
ylabel('Distances [AU]')
legend('Earth-Mars distance', 'Equivalent distance through Satellite 1', 'Equivalent distance through Satellite 2');

%% Functions

% function Out = wrapperFunc(X1,X2,yi)
% global etR muSu OPTSCALING
%     [~,XSat] = ode15s(@(t,y) twobodyode(t,y,muSu), etR, yi);
%     XSats = XSat';
%     XSats(:,:,2) = -XSat';
%     Out = bestLinkBudget(X1,X2,XSats,@computeLinkBudget)*OPTSCALING;
% end


function Out = wrapperFunc(X1,X2,yi)
    global etR muSu OPTSCALING
    [~, nT] = size(etR);
    interval = (etR(nT)-etR(1))/(nT-1);
    XSat(:,1) = yi;
    XSat(:,nT) = yi;
    x0_C = yi;
    for i = 1:nT
        [x1_C, ~] = ds_2bpSTM(x0_C,interval,muSu);
        XSat(:,i) = x1_C;
        x0_C = x1_C;
    end
    XSats(:,:,2) = -XSat;
    Out = bestLinkBudget(X1,X2,XSats,@computeLinkBudget)*OPTSCALING;
end

function Out = bestLinkBudget(X1,X2,XSats,PairwiseFunc)
    minout = 10^25;
    [~, nT, numSats] = size(XSats);
    for i = 1:numSats
        LB1 = PairwiseFunc(X1,XSats(:,:,i));
        LB2 = PairwiseFunc(X2,XSats(:,:,i));
        minout = min(minout, max(LB1,LB2));
    end
    Out = minout;
end
    