function [m, n, Rsqr, rmse] = powerLawFit(gamma, mu, weighted)
% powerLawFit Fits a power‐law model to data
%
% Syntax:
%   [m, n, Rsqr, rmse] = powerLawFit(gamma, mu)
%   [m, n, Rsqr, rmse] = powerLawFit(gamma, mu, weighted)
%
% Description:
%   Estimates parameters of the power‐law relationship
%       mu = m * gamma^n
%   by performing linear regression on log‐transformed data. Optionally,
%   a weighted least squares fit can be applied.
%
% Inputs:
%   gamma    - Numeric vector of independent variable observations
%   mu       - Numeric vector of dependent variable observations (same size as gamma)
%   weighted - (optional) Logical flag to enable weighted fitting (default: false).
%              When true, weights are set to the inverse variance of log‐residuals.
%
% Outputs:
%   m    - Estimated scaling coefficient of the power law
%   n    - Estimated exponent of the power law
%   Rsqr - Coefficient of determination (goodness‐of‐fit)
%   rmse - Root‐mean‐square error of residuals on the original data scale
%
% Example:
%   gamma = logspace(0, 2, 50);
%   mu = 3.5 * gamma.^(-1.2) + 0.2 * randn(size(gamma));
%   [m, n, R2, err] = powerLawFit(gamma, mu, true);

    % Check arguments
    if nargin < 3
        weighted = false;
    end  

    % Build interest quantities
    if weighted
        mu_sqr = mu.^(2);
    else
        mu_sqr = ones(size(mu));
    end
    ln_gamma = log(gamma);
    ln_mu    = log(mu);

    % Build linear system
    a11 = sum(mu_sqr);
    a12 = sum(mu_sqr.*ln_gamma);
    a22 = sum(mu_sqr.*(ln_gamma.^2));
    A = [a11, a12;
         a12, a22];
    b = [sum(mu_sqr.*ln_mu);
        sum(mu_sqr.*ln_mu.*ln_gamma)];

    % Solve linear system
    params = A\b;
    M = params(1);
    N = params(2);

    % Get parameters for the exponential model
    m = exp(M);
    n = N + 1;

    % Evaluate errors (nonlinear regression)
    SSres = sum((mu - m*gamma.^(n - 1)).^2);
    SStot = sum((mu - mean(mu)).^2);
    Rsqr = 1.0 - SSres/SStot;
    rmse = sqrt(sum((mu - m*gamma.^(n - 1)).^2)/numel(mu));

end

