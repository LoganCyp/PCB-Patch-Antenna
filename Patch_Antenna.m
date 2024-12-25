%% OshPark 915 MHz Patch Antenna Design
% Logan Cypser
clear;

%% Physical constant
% Center Frequency
fc = 915e6; 
% Dielectric Constant
er = 4.35;
% Speed of Light 
c = physconst('LightSpeed');
% Free Space Wave Number 
k = (2 * pi * fc) / c;
% Intrinsic Impedance of Free Space
eta = 377;

Z = 50;
%% Design Parameters

% Substrate Height
h = 1.524e-3;

F = (Z / 60) * sqrt((er + 1) / 2) + ((er - 1) / (er + 1)) * (0.23 + (0.11 / er));

if F < 1
    W_h = (8 * exp(F)) / (exp(2 * F) - 2);
else 
    B = (eta * pi) / (2 * Z * sqrt(er));
    W_h = (2 / pi) * (B - 1 * log(2*B - 1) + ((er - 1) / (2*er)) * (log(B -1) + 0.39 - (0.61 / er)));
end

W0 = W_h * h;

% Stripeline Width @ 50 Ohms
mw50 = W0; %m

% Patch Width
W = c / (2 * fc * sqrt((er + 1) / 2));

% Effective Dielectric Constant
e_eff = ((er + 1) / 2) + ((er - 1) / 2) * (1 + (12 * h) / W) ^ (-1/2);

% Effective length
L_eff = c / (2 * fc * sqrt(e_eff));

% Length Extension
d_L = (0.412 * h) * (((e_eff + 0.3) * ((W / h) + 0.264)) / ((e_eff - 0.258) * ((W / h) + 0.8)));

% Adjusted Length
L = L_eff - d_L;

% Input Impedance
%Za = 90 * (er^2 / (er - 1)) * (L / W)^2;

% Conductance
X = k * W;
G1 = (1 / (pi * eta)) * (-2 + cos(X) + X*sinint(X) + (sin(X) / X));

% Mutual Conductance
integrand = @(theta) (sin(((k * W) / 2) .* cos(theta)) ./ cos(theta)).^2 .* besselj(0, k * L .* sin(theta)) .* sin(theta).^3;
G12 = 1/(120 * pi^2) * integral(integrand, 0, pi);

% Input Resistance
Rin = 1 / (2 * (G1 + G12));

% Inset Feed Point Distance
y0 = (L / pi) * acos(sqrt(Z / Rin));

% Ground Plane Length
Lg = 2 * L; 

% Ground Plane Width
Wg = 2 * W; 

%Inset notch width
n = 0.3 * W0;
Nw = (n*2) + W0;

%Inset notch length
Nl = y0;

%% Create antenna primatives
patch = antenna.Rectangle('Length', L, 'Width', W);
groundplane = antenna.Rectangle('Length', Lg, 'Width', Wg);
notch = antenna.Rectangle('Length', Nl, 'Width', Nw, 'Center', [(L/2)-(Nl/2),0]);
microstrip_feed = antenna.Rectangle('Length', Lg/2, 'Width', mw50, 'Center', [(Lg/4),0]);
substrate_material = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h); 

build_patch = (patch-notch) + microstrip_feed;

%% Define the properties of the PCB stack.
basicPatch = pcbStack;
basicPatch.Name = 'CACI 915 MHz Antenna';
basicPatch.BoardThickness = h;
basicPatch.BoardShape = groundplane;
basicPatch.Layers = {build_patch,substrate_material,groundplane};
basicPatch.FeedLocations = [Lg/2 0 1 3];
basicPatch.FeedDiameter = mw50/2;
figure
show(basicPatch)
figure;
%mesh(basicPatch, 'MaxEdgeLength',2.5e-3,'MinEdgeLength',0.8e-3);
%mesh(basicPatch, 'MaxEdgeLength',1e-2, 'MinEdgeLength',2e-3);
%mesh(basicPatch, 'MaxEdgeLength',1e-2, 'MinEdgeLength',2e-3);
mesh(basicPatch, 'MaxEdgeLength', 7e-3);

%% Plot the radiation pattern of the basic patch antenna.
figure
pattern(basicPatch, fc)

%% Plot the impedance of the basic patch antenna.
%enumerate frequencies
%freqs = linspace(fc-0.05*fc,fc + 0.1*fc,10);
freqs = linspace(875e6, 975e6, 100);


%plot complex impedance
figure
impedance(basicPatch,freqs)

%plot RF s-parameters
S = sparameters(basicPatch, freqs);
figure; 
rfplot(S);

%connector
connector = PCBConnectors.SMAEdge;
connector.SignalLineWidth = mw50;
connector.EdgeLocation = 'east';
connector.ExtendBoardProfile = false;

%pcb service
service = PCBServices.OSHParkWriter;
service.Filename = 'OshPark_915MHz_Patch.zip';

%write gerber
PW = PCBWriter(basicPatch,service,connector);
gerberWrite(PW);