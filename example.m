% The example code. 
clear;
% close all;
%%
% Model coefficents
n    = 100;     % accuracy parameter
m    = 30;      % number of spanwise modes

% flow parameters
c0   = 343;     % the speed of sound
rho0 = 1.25;    % the density of the fluid
M0   = 0.1;     % Mach number of the flow

% aerofoil parameters
c       = 1;                        % the chord length
r_lam_h = 8;                        % ratio of lambda to h
r_h_c   = 0.025;                    % ratio of h to c
kc      = 10.^linspace(-1,2,100);   % dimensionless frequencies required
omega   = kc*c0/c;                  % dimensional angular frequencies
freq    = omega/2/pi;               % dimensional frequencies

obs_location = [0, 0, c];           % obsever location vector

%% Use the function Pred_Green.m to obtain the code for predicting the PSD of the TE noise under the specified sawtooth profile and flow parameters.
h      = c*r_h_c;                   % half root-to-tip amplitude
lambda = r_lam_h*h;                 % serration wavelength

spec_saw  = Pred_Green(lambda, h*2, 1, M0, rho0, obs_location, freq,c0,n,m); % Sawtooth with lambda and h
spec_base = Pred_Green(0.04,  1e-6, 1, M0, rho0, obs_location, freq,c0,n,m); % baseline


%% Figure
% Code for drawing the figure
Line_color = [248 118 109;  190	190	190]/256; % Line color RGB used in figures
fig_phi = figure;
pic_phi = tiledlayout("flow",'TileSpacing','compact','Padding','compact');

nexttile(pic_phi);
semilogx(kc,10*log10(spec_base/4e-10), 'Color',Line_color(1,:),'LineStyle','-','LineWidth',2); hold on;
semilogx(kc,10*log10(spec_saw/4e-10),'Color',Line_color(2,:),'LineStyle','-','LineWidth',2); hold on;
ylim([-60 30]);
yticks([-60 -30 0]);
xlim([0.1 100]);
legend('Baseline','Serrated','Location','southeast');
xlabel('$kc$','Interpreter','latex','FontName','Time New Roman');
ylabel('PSD (dB ref 4\times 10^{-10}Pa^2/Hz)');
fontsize(pic_phi,22,"pixels");
fontname(pic_phi,"Times New Roman");