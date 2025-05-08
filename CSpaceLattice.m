%% Cspace Lattice Performance Analysis
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

%% Create Lattice

a = 1.25*AU;
i = 0;
Om = 0;
maxE = 0.2;

numEs = 10;
numws = 10;
numf0s = 10;
numSats = numEs*numws*numf0s;

XSats = [];

for e = linspace(0,maxE,numEs)
    for w = linspace(0,2*pi,numws)
        for f0 = linspace(0,2*pi,numf0s)
            Ki = [a, e, i, Om, w, f0];
            X = propagateFromKeplerians(Ki,muSu,etR);
            XSats = cat(3, XSats, X);
        end
    end
end

%% Plotting

DSC_TVG_Plotting_from_Xs(20, XSats, XEa, XMa, "Earth", "Mars", 'blue', 'red', etR)