% HPW_fit.m
% 
% Function to calculate viscosity (eta) and inclusion radius (a)
% from translational (DT) and rotational (DR) diffusion coefficients, 
% by fitting to the exact equations of
%    B. D. Hughes, B. A. Pailthorpe, White L. R.,
%    The translational and rotational drag on a cylinder moving in a membrane. 
%    J Fluid Mech. 110, 349-372 (1981)
% (see HPW_drag.m)
%
% Finds the eta, a that minimize the squared distance between calculated
% and experimental Dr, Dt.
% Distance^2: scaleD*scaleD*(calc_DT-Dt)^2 + (calc_DR-Dr)^2;
% scaleD = 1e12 (hard-coded) to account for scale difference between
% typical DT and DR values.
%
% Inputs:
%    Dr = rotational diffusion coefficient (rad^2/s)
%    Dt = translational diffusion coefficient (m^2/s)
%    eta_w = external (water) viscosity, default 8.9e-4 Pa.s
%    T = temperature (K), default 295;
%    nTerms =  number of terms to calculate in "infinite" series; default 30
%    num_lterms = number of "l" terms to calculate in ; default 20
%    num_mterms = number of "l" terms to calculate in ; default 15;
%    params0 =  initial values for log10_eta, log10_a (/Pa.s, m); default [-9, -6];
%    LB     : [optional] lower bounds of search parameters(log10_eta,
%             log10_a); default [-10 -9]
%    UB     : [optional] upper bounds of search parameters(log10_eta,
%             log10_a); default [-7 -4]
%
% Outputs:
%   etafit: best-fit viscosity, Pa.s
%   afit: best-fit inclusion radius, m
%   chi2 : squared distance with best-fit parameters
%
%%
% Copyright 2018, Raghuveer Parthasarathy, The University of Oregon
%
%%
% Disclaimer / License  
%   This program is free software: you can redistribute it and/or 
%     modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.
%   This set of programs is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
%   General Public License for more details: <http://www.gnu.org/licenses/>.
%
% Raghuveer Parthasarathy
% Aug. 29, 2018
% last modified Sept. 2, 2018

function [etafit, afit, chi2fit] = HPW_fit(Dr, Dt, eta_w, Temperature, nTerms, num_lterms, num_mterms, params0, LB, UB)

if ~exist('eta_w', 'var') || isempty(eta_w)
    eta_w = 8.9e-4; % water viscosity, Pa.s
end
if ~exist('Temperature', 'var') || isempty(Temperature)
    Temperature = 295; % K
end
if ~exist('nTerms', 'var') || isempty(nTerms)
    nTerms = 30;
end
if ~exist('num_lterms', 'var') || isempty(num_lterms)
    num_lterms = 20;
end
if ~exist('num_mterms', 'var') || isempty(num_mterms)
    num_mterms = 15;
end
if ~exist('params0', 'var') || isempty(params0)
    params0 = [-8.5 -7];
end
if ~exist('LB', 'var') || isempty(LB)
    LB = [-10 -9];
end
if ~exist('UB', 'var') || isempty(UB)
    UB = [-7 -4];
end

% find best-fit eta, a; defined by minimal distance between calc. and
% experimental (Dt,Dr)

% Because Dt << Dr (in SI units), multiply the Dt chi^2 by this factor:
scaleD = 1e12;

% optimization options
lsqoptions.TolFun = 1e-8;  %  % MATLAB default is 1e-6
lsqoptions.TolX = 1e-8;  % default is 1e-6
lsqoptions.Display = 'final'; % 'off' or 'final'; 'iter' for display at each iteration
params = lsqnonlin(@(P) objfun(P, scaleD, Dr, Dt, eta_w, Temperature, nTerms, num_lterms, num_mterms),params0,LB,UB,lsqoptions);

etafit = 10^params(1);
afit = 10^params(2);
chi2fit = sum(objfun(params, scaleD, Dr, Dt, eta_w, Temperature, nTerms, num_lterms, num_mterms));

plotopt = false;
if plotopt
    figure; plot(Dr, Dt, 'o', 'markersize', 14)
    hold on
    [calc_D_T, calc_D_R] = HPW_drag(10^params0(1), 10^params0(2), eta_w, Temperature, nTerms, num_lterms, num_mterms);
    plot(calc_D_R, calc_D_T, 'x', 'markersize', 14)
    [calc_D_T, calc_D_R] = HPW_drag(10^params(1), 10^params(2), eta_w, Temperature, nTerms, num_lterms, num_mterms);
    plot(calc_D_R, calc_D_T, 's', 'markersize', 14)
    xlabel('D_R (rad^2/s)')
    ylabel('D_T (m^2/s)')
    legend('Measured', 'Initial params', 'Best fit')
end

end

    function chi2 = objfun(params, scaleD, Dr, Dt, eta_w, Temperature, nTerms, num_lterms, num_mterms)
    %  distance between calc. and experimental (Dt,Dr)
    [calc_DT, calc_DR] = HPW_drag(10^params(1), 10^params(2), eta_w, Temperature, nTerms, num_lterms, num_mterms);
    % chi2 = scaleD*scaleD*(calc_DT-Dt)^2 + (calc_DR-Dr)^2;
    chi2 = [scaleD*scaleD*(calc_DT-Dt)^2  (calc_DR-Dr)^2];
    % param_dist = scaleD*abs(calc_DT-Dt);
    % param_dist = abs(calc_DR-Dr)
    end

