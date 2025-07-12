% Measurements obtained from
% R. E. Wells and E. W. Merrill, "Influence of Flow Properties of Blood 
% Upon Viscosity-Hematocrit Relationships," J Clin Invest, vol. 41, no. 8, 
% pp. 1591–1598, Aug. 1962, doi: 10.1172/JCI104617.
clear; clc; close all


% Add function to path
addpath("src/");

% Conversion factors
cp_to_p   = 0.01;
cp_to_PaS = 0.001;

% Colors for plots
colors = hsv(5);

% Import power-law parameters obtained from Wells measures (obtained by
% running fit_powerlaw_to_wells.m
load("fitted_parameters.mat");


%% Fitting
% Shear rate range for fitting
x = linspace(10,130,100)'; % vector de shear rate de 10 a 130 sec^-1
% x = linspace(1,2000,100)'; % vector de shear rate de 10 a 130 sec^-1

% HCT values for testing
% HCTvalues = [28.2, 35.2, 40.2, 46.6, 50.1];
HCTvalues = [20.0, 32.5, 43.0, 57.5, 70.0];

% HCT values and viscosities of synthetic measurements
HCTvalues_meas = [12.2892, 24.4337, 61.8313, 123.6627];
mu = zeros([numel(HCTvalues), numel(HCTvalues_meas)]);

% Model fitting for each HCT
m = zeros(size(HCTvalues));
n = zeros(size(HCTvalues));
legends = cell(1,numel(HCTvalues));

for i=1:numel(HCTvalues)
    hct = sprintf('Hct = %.1f',HCTvalues(i));
    legends{i} = hct;

    % Fit
    [m_fitted, n_fitted] = powerLawParams(HCTvalues(i), hct_fit, m_fit, n_fit);
    m(i) = m_fitted;
    n(i) = n_fitted;

    % Evaluate synthetic measurements
    mu(i,:) = m(i)*HCTvalues_meas.^(n(i)-1);

end

% Plots
figure(1);
% mu_newtonian = mean(mu, 2);
mu_newtonian = zeros([1, numel(m)]);
range = [0, 2800];
for i=1:numel(m)    
    area = m(i)/n(i)*(range(2)^n(i) - range(1)^n(i));
    mu_newtonian(i) = area/(range(2)-range(1));
end

for i=1:numel(HCTvalues)
    % % Newtonian viscosities
    % plot(x, mu_newtonian(i)*ones(size(x)),':','LineWidth',4,'Color',colors(i,:)); hold on

    % Evaluate and store fitted model
    plot(HCTvalues_meas, mu(i,:), 'x', 'LineWidth', 4, 'MarkerSize', 24, 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', colors(i,:), 'Color', colors(i,:)); hold on
    plot(x, m(i)*x.^(n(i)-1),'-','LineWidth',4,'Color',colors(i,:)); hold on
end

% Dummy plot for legend
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(5,:), 'MarkerFaceColor', colors(5,:)); hold on
plot(NaN, NaN, 'x', 'LineWidth', 4, 'MarkerSize', 24, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k'); hold on
% plot(NaN, NaN,':','LineWidth', 4,'Color','k'); hold on
plot(NaN, NaN,'-','LineWidth', 4,'Color','k'); hold off

% legend('', '', '', '', '', '', '', '', '', '', '', '', '', '','', 'Hct 20', 'Hct 32.5', 'Hct 45', 'Hct 57.5', 'Hct 70', 'Interp. meas.', '$\mu_{\mathrm{NF}}$', '$\mu_{\mathrm{PL}}$', 'Interpreter', 'latex', 'fontsize', 35, 'numcolumns',1)
legend('', '', '', '', '', '', '', '', '','', 'Hct 20', 'Hct 32.5', 'Hct 45', 'Hct 57.5', 'Hct 70', 'Interp. meas.', '$\mu_{\mathrm{PL}}$', 'Interpreter', 'latex', 'fontsize', 35, 'numcolumns',1)

xlabel("$\dot{\gamma}$ (1/s)", "Interpreter", "latex")
ylabel("$\mu$ (centipoise)", "Interpreter", "latex")

ax = gca;
% ax.XAxis.TickValues = [0.5 1 1.5 2.0 2.5];
% ax.XAxis.TickLabels = [20, 32.5, 45, 57.5, 70];
% ax.XAxis.TickLabels = ["Hct 28.2", "Hct 35.2", "Hct 40.2", "Hct 46.6", "Hct 50.1"]; % 28.2, 35.2, 40.2, 46.6, 50.1
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.06;
ax.Box = 'on';
ax.FontWeight = 'bold';
ax.FontSmoothing = 'on';
ax.FontSize = 35;
ax.LineWidth = 1.5;
ax.XAxis.MinorTick = 'off';
ax.TickLength = [0.025, 0.25];
ax.XAxis.TickLength = [0.0, 0.25];
ax.TickDir = 'in';
ax.TickLabelInterpreter = 'latex';
xlim([5, 135]);

set(gca, 'Position', [0.130000000000000,0.123058425254726,0.775000000000000,0.801941574745274])
set(gcf, 'Position', [19.222222222222220,5.552222222222222e+02,6.902222222222222e+02,5.177777777777778e+02])
% set(gca, 'Position', [0.130000000000000,0.173260341285209,0.775000000000000,0.751739658714791])
% set(gcf, 'Position', [409,229,836,638])