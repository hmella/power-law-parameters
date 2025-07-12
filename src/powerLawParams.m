function [m_target, n_target] = powerLawParams(HCT_target, hct_fit, m_fit, n_fit)
%% powerLawParams
% powerLawParams  Compute interpolated power‐law parameters for a given hematocrit
%
% Syntax:
%   [m_target, n_target] = powerLawParams(HCT_target, hct_fit, m_fit, n_fit)
%
% Inputs:
%   HCT_target  – double scalar, desired hematocrit value for which to obtain
%                 the power‐law coefficients
%   hct_fit     – double vector, hematocrit values at which the model was fitted
%   m_fit       – double vector, fitted m coefficients corresponding to hct_fit
%   n_fit       – double vector, fitted n coefficients corresponding to hct_fit
%
% Outputs:
%   m_target    – double scalar, interpolated m coefficient at HCT_target
%   n_target    – double scalar, interpolated n coefficient at HCT_target
%
% Description:
%   This function takes a set of power‐law model coefficients (m_fit, n_fit)
%   obtained at known hematocrit levels (hct_fit) and interpolates them to
%   compute the coefficients (m_target, n_target) corresponding to the specified
%   HCT_target. Internally, a 1-D interpolation (e.g., interp1) is used.
%
% Example:
%   % Interpolate parameters for HCT = 0.45
%   hct_vals = [0.30, 0.35, 0.40, 0.50];
%   m_vals   = [1.10, 1.15, 1.20, 1.30];
%   n_vals   = [0.90, 0.92, 0.94, 0.98];
%   [m_t, n_t] = powerLawParams(0.45, hct_vals, m_vals, n_vals);

  % Check Hct range
  if HCT_target < 16 || HCT_target > 70
      error('Invalid Hct value: %d is out of the range [16, 70]', HCT_target);
  end


  % Power-law parameters fitted from measurements
  hct_fit = [16, 33, 43, 57, 70];

  % Shear rate range for fitting
  x = linspace(10,130,100)'; % vectror de shear rate de 10 a 130 sec^-1

  % Adjusted viscosities
  mu_adj = struct('h16', struct('x', x, 'y', zeros(size(x))),...
                'h33', struct('x', x, 'y', zeros(size(x))),...
                'h43', struct('x', x, 'y', zeros(size(x))),...
                'h57', struct('x', x, 'y', zeros(size(x))),...
                'h70', struct('x', x, 'y', zeros(size(x))));

  % Target viscosities
  mu_target = zeros(size(x));
  for i=1:numel(hct_fit)
      % Evaluate and store fitted model
      hct = sprintf('h%d',hct_fit(i));
      mu_adj.(hct).y = m_fit(i)*x.^(n_fit(i)-1);
  end

  % Build data points to adjust viscosity model to target HCT
  for i=1:numel(x)
      % Extract values corresponding to shear rate xi
      tmp = zeros(size(hct_fit));
      for j=1:numel(hct_fit)
          hct = sprintf('h%d',hct_fit(j));
          tmp(j) = mu_adj.(hct).y(i);
      end
  
      % Interpolate value to obtain target Hct
      mu_target(i) = interp1(hct_fit, tmp, HCT_target);
  end
  
  % Adjust viscosity model to obtain m and n value
  weighted = true;
  [m_target, n_target, Rsqr, rmse] = powerLawFit(x, mu_target, weighted);

end

