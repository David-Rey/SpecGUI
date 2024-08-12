
clear; clc; close all;
rng(1);

% Path to raman spectra .spe file
ramanDataPath = '../data/2017 février 21 17_02_49.spe';

% Instrument function data path
instrumentPath = '../data/instrument/Fct_instrument_1BIN_2400g.csv';

% Define and create file paths from file number for thomson signal
fileNum = 24;
baseDir = fullfile('..', 'data');

% Construct the full paths to the signal and background files
signalName = fullfile(baseDir, ['2020-07-21  ', num2str(fileNum), '.spe']);
backgroundName = fullfile(baseDir, ['2020-07-21  ', num2str(fileNum + 1), '.spe']);

% Energy in eV
Bev = 2.48e-4;

% Incident light wavelength in nm
centerWavelength = 532;

% Temperature in K
environment.temperature = 288.15;

% Pressure in Pa
environment.pressure = 1000;

% Volume in m^3
environment.volume = 0.90478;

% Laser power in W
environment.power = 4;

psoOptinos.numParticles = 500;
psoOptinos.numIterations = 50;

% Create Optimzation object
opt = Optimize(ramanDataPath, instrumentPath, Bev, environment, centerWavelength, psoOptinos);

% perform curve fitting operation
opt.optimize();

% create Thomson object to calculate 
tom = Thomson(signalName, backgroundName, psoOptinos);

fig1 = figure();
ax1 = gca;
opt.drawLandscape(ax1);

fig2 = figure();
ax2 = gca;
opt.drawOverlay(ax2);

fig3 = figure();
ax3 = gca;
tom.drawThomson(ax3);

electronDensity = tom.area_SI * opt.bestScale;
fprintf('The electron density is: %.3e\n', electronDensity)



