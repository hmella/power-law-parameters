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


%% READ DATA
% Import digitalized measurements
data = readtable('wells_measures/wpd_datasets.csv','decimal', '.');
mu = struct('h16', struct('x', data.Hct16, 'y', data.Var2),...
             'h33', struct('x', data.Hct33, 'y', data.Var4),...
             'h43', struct('x', data.Hct43, 'y', data.Var6),...
             'h57', struct('x', data.Hct57, 'y', data.Var8),...
             'h70', struct('x', data.Hct70, 'y', data.Var10));
HCTvalues = [16, 33, 43, 57, 70];


%% LEAST-SQUARES FITTING
% Weighted fit?
weighted = true;

% Shear rate range for fitting
x = linspace(10,130,100)';

% Adjusted viscosities
mu_adj = struct('h16', struct('x', x, 'y', zeros(size(x))),...
             'h33', struct('x', x, 'y', zeros(size(x))),...
             'h43', struct('x', x, 'y', zeros(size(x))),...
             'h57', struct('x', x, 'y', zeros(size(x))),...
             'h70', struct('x', x, 'y', zeros(size(x))));
RMSE = zeros([1, numel(HCTvalues)]);
Rsqr = zeros([1, numel(HCTvalues)]);

% Model fitting for each HCT
m = zeros(size(HCTvalues));
n = zeros(size(HCTvalues));
for i=1:numel(HCTvalues)
    hct = sprintf('h%d',HCTvalues(i));

    % Fit
    [m_fitted, n_fitted, Rsquare, rmse] = powerLawFit(mu.(hct).x, mu.(hct).y, weighted);
    fprintf('\nFit for Hct %d. R^2 = %d', HCTvalues(i), Rsquare)
    m(i) = m_fitted;
    n(i) = n_fitted;

    % Evaluate and store fitted model
    mu_adj.(hct).y = m(i)*x.^(n(i)-1);
    RMSE(i) = rmse;
    Rsqr(i) = Rsquare;
end
fprintf('\n')

% Export power-law parameters to mat file
m_fit = m;
n_fit = n;
hct_fit = HCTvalues;
save('fitted_parameters.mat', 'm_fit', 'n_fit', 'hct_fit');

% Plot
figure(1);
for i=1:numel(HCTvalues)
    % Current Hct
    hct = sprintf('h%d',HCTvalues(i));

    % Calculate error
    error = sqrt( sum( ( mu.(hct).y - m(i)*mu.(hct).x.^(n(i) - 1) ).^2 )/( numel(mu.(hct).y) - 2 ) );

    % Plot measurements
    plot(mu.(hct).x, mu.(hct).y, 'x', ...
        'LineWidth', 4, ...
        'MarkerSize', 24, ...
        'MarkerEdgeColor', colors(i,:), ...
        'MarkerFaceColor', colors(i,:)); 
    hold on

    % plot fitted model
    plot(x, m(i)*x.^(n(i)-1), '-', ...
        'LineWidth', 4, ...
        'Color',colors(i,:), ...
        'DisplayName', sprintf('Hct %d', HCTvalues(i)));
    hold on
end

% Dummy plots for legends
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(4,:), 'MarkerFaceColor', colors(4,:)); hold on
plot(NaN, NaN, 's', 'LineWidth', 2, 'MarkerSize', 24, 'MarkerEdgeColor', colors(5,:), 'MarkerFaceColor', colors(5,:)); hold on
plot(NaN, NaN, 'x', 'LineWidth', 4, 'MarkerSize', 24, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k'); hold on
plot(NaN, NaN, '-', 'LineWidth', 4, 'Color', 'k'); hold off

legend('', '', '', '', '', '', '', '', '', '',...
    sprintf('Hct 16 ($R^2$ = %.3f, RMSE = %.2f)', Rsqr(1), RMSE(1)),...
    sprintf('Hct 33 ($R^2$ = %.3f, RMSE = %.2f)', Rsqr(2), RMSE(2)),...
    sprintf('Hct 43 ($R^2$ = %.3f, RMSE = %.2f)', Rsqr(3), RMSE(3)),...
    sprintf('Hct 57 ($R^2$ = %.3f, RMSE = %.2f)', Rsqr(4), RMSE(4)),...
    sprintf('Hct 70 ($R^2$ = %.3f, RMSE = %.2f)', Rsqr(5), RMSE(5)),...
    'Wells et al.', '$\mu_{\mathrm{PL}}$', ...
    'Interpreter', 'latex', 'fontsize', 35, 'numcolumns', 1, 'location', 'northeast')

xlabel("$\dot{\gamma}$ (1/s)", "Interpreter", "latex")
ylabel("$\mu$ (centipoise)", "Interpreter", "latex")

ax = gca;
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

% set(gca, 'Position', [0.130000000000000,0.123058425254726,0.775000000000000,0.801941574745274])
% set(gcf, 'Position', [19.222222222222220,5.552222222222222e+02,6.902222222222222e+02,5.177777777777778e+02])