%% LLM_fit_viscosity_L.m
%
% Calulates membrane viscosity and rod length for a rod-like membrane
% inclusion, using the theory of 
%    A. J. Levine, T. B. Liverpool, F. C. MacKintosh, 
%    Mobility of extended bodies in viscous films and membranes. 
%    Phys. Rev. E. 69, 021503 (2004)
% Use only D_parallel and D_R, not D_perpendicular (redundant; noisy)
%
% Note: extracted points from LLM 2004 graphs, and make simple polynomial
% interpolations. See notes, 2021, and Jahl and Parthasarathy 2021.
%
% For uncertainties, repeated sampling from a Gaussian distribution with
% the mean, std. deviation of the D values.
%
% Inputs:
%    D_par : Parallel diffusion coefficient (m^2/s)
%    unc_D_par : uncertainty in Parallel diffusion coefficient (m^2/s)
%    D_R : Rotational diffusion coefficient (rad^2/s)
%    unc_D_R : uncertainty in Rotational diffusion coefficient (rad^2/s)
%    N_trials : number of trials for repeated sampling, for uncertainty.
%       (Default 1, so no uncertainty calculation)
%    eta_w = external (water) viscosity, default 8.9e-4 Pa.s
%    Temperature = temperature (K), default 295;
%    plotopt: if True, make plots (intersecting Dr and Dparl curves on eta vs.
%             L axes; histogram of eta values if sampling multiple values
%             for uncertainties); default false
%
% Outputs:
%    eta_m : Membrane viscosity (Pa.s.m)
%    stdDev_eta_m : standard deviation of membrane viscosity from 
%            repeated sampling (Pa.s.m) -- note that 1/2 of the 68%
%            confidence interval (0.5*(eta_m_84perc-eta_m_16perc)) is a
%            better measure of uncertainty -- see July 22, 2021 notes.
%    L : Bead length (m)
%    unc_L : uncertainty in Bead length (m)
%    eta_m_16perc, eta_m_84perc : 16th and 84th percentile eta_m values,
%        from sampling for uncertainty estimates; empty if N_trials = 1.
%        Note that for a Gaussian distribution these are +/- 1 sigma, but
%        the eta_m distributions are non-Gaussian (see 22 July 2021 notes)
%    eta_m_each : all the eta_m values from repeated sampling
%    L_each : all the eta_m values from repeated sampling
%
%%
% Copyright 2021, Raghuveer Parthasarathy, The University of Oregon
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
%%
% Based on Philip Jahl: Levine_solve_parl.m
% Raghuveer Parthasarathy
% June 30, 2021
% last modified July 22, 2021


function [eta_m, stdDev_eta_m, L, unc_L, eta_m_16perc, eta_m_84perc, eta_m_each, L_each] = ...
    LLM_fit_viscosity_L(D_par, unc_D_par, ...
    D_R, unc_D_R, N_trials, Temperature, eta_w, plotopt)

if ~exist('N_trials', 'var') || isempty(N_trials)
    N_trials = 1;
end
if ~exist('eta_w', 'var') || isempty(eta_w)
    eta_w = 8.9e-4; % water viscosity, Pa.s
    % Or use:     eta_w = 1.01e-3; % 0.1M sucrose viscosity, Pa.s
end
if ~exist('Temperature', 'var') || isempty(Temperature)
    Temperature = 295; % K
end
if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = 295; % K
end

k_B = 1.380649e-23; % Boltzmann's constant
eta_seed = 5e-9;
L_seed = 1e-6;

% Make sure that D_parl is in m^2/s, not um^2/s
if D_par > 1e-10
    fprintf('Is D_parl in m^2/s? \n')
end

syms eta_m_eq L_eq positive % "_eq" added for symbolic variables
assumeAlso(eta_m_eq, 'real')
assumeAlso(L_eq, 'real')

%% Diffusion coefficient equations

Dr_equation = D_R*fig11_fit(eta_w, L_eq, eta_m_eq) == k_B*Temperature/(4*pi*eta_m_eq*L_eq^2);
D_parl_equation = D_par*fig3_fit(eta_w, L_eq, eta_m_eq) == k_B*Temperature/(4*pi*eta_m_eq);

[eta_m, L] = vpasolve([Dr_equation, D_parl_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
eta_m = double(eta_m); % Necessary because it's possible to have no solution
L = double(L);
% Note that these will be overwritten if we're sampling a distribution.

if isempty(eta_m)
    eta_m = NaN;
    L = NaN;
end


%% Repeated calculations for uncertainties
if N_trials>1
    eta_m_each = zeros(N_trials,1);
    L_each = zeros(N_trials,1);
    rand_Dparl = D_par + unc_D_par*randn(N_trials,1);
    rand_Dr = D_R + unc_D_R*randn(N_trials,1);
    for j=1:N_trials
        if mod(j,100)==0
            fprintf('trial %d of %d.\n', j, N_trials);
        end
        this_D_parl = rand_Dparl(j);
        this_Dr = rand_Dr(j);
        % I think I need to re-define the equations...
        Dr_equation = this_Dr*fig11_fit(eta_w, L_eq, eta_m_eq) == k_B*Temperature/(4*pi*eta_m_eq*L_eq^2);
        D_parl_equation = this_D_parl*fig3_fit(eta_w, L_eq, eta_m_eq) == k_B*Temperature/(4*pi*eta_m_eq);
        [sampled_eta_m_out, sampled_L_out] = vpasolve([Dr_equation, D_parl_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
        if isempty(sampled_eta_m_out)
            eta_m_each(j) = NaN;
            L_each(j) = NaN;
        else
            eta_m_each(j) = double(sampled_eta_m_out);
            L_each(j) = double(sampled_L_out);
        end
    end
    nanInd = isnan(eta_m_each);
    eta_m_each(nanInd)=[];
    L_each(nanInd)=[];
    
    % Take the mean of the sampled eta and L as the output values, rather
    % than the value calculated from the most probable D's, above.
    % Omit NaNs (failed solutions)
    eta_m = mean(eta_m_each, 'omitnan');
    L = mean(L_each, 'omitnan');
    stdDev_eta_m = std(eta_m_each, 'omitnan');
    unc_L = std(L_each, 'omitnan');
    % percentiles, ignoring NaNs
    sort_eta = sort(eta_m_each(~isnan(eta_m_each)));
    N_trials_notNaN = length(sort_eta);
    index_16 = round(0.16*N_trials_notNaN);
    index_84 = round(0.84*N_trials_notNaN);
    if index_16 > 1 && index_84 < N_trials_notNaN
        ok_percentiles = true;
        eta_m_16perc = sort_eta(round(0.16*N_trials_notNaN));
        eta_m_84perc = sort_eta(round(0.84*N_trials_notNaN));
    else
        ok_percentiles = false;
        eta_m_16perc = NaN;
        eta_m_84perc = NaN;
    end
    
    if plotopt
        figure; histogram(eta_m_each(~isnan(eta_m_each))); xlabel('\eta_m (Pa.s.m)')
        hold on
        ax = axis;
        plot(median(eta_m_each, 'omitnan')*[1 1], [ax(3) ax(4)], '--', 'color', [0.7 1.0 0.4], 'linewidth', 2.0) 
        plot(mean(eta_m_each, 'omitnan')*[1 1], [ax(3) ax(4)], '-', 'color', [0.9 0.6 0.2], 'linewidth', 2.0) 
        % plot(eta_m*[1 1], [ax(3) ax(4)], '--', 'color', [0.9 0.6 0.2], 'linewidth', 2.0) 
        plot(mean(eta_m_each, 'omitnan') + stdDev_eta_m*[1 1], [ax(3) ax(4)], ':', 'color', [0.9 0.6 0.2]) 
        plot(mean(eta_m_each, 'omitnan') - stdDev_eta_m*[1 1], [ax(3) ax(4)], ':', 'color', [0.9 0.6 0.2]) 
        legend('histogram', 'median', 'mean', 'mean + 1 std.', 'mean - 1 std.')
        
        % Ranked plot of eta values
        figure; plot(sort(eta_m_each), 'o')
        hold on
        plot([1 N_trials], mean(eta_m_each, 'omitnan')*[1 1], '-', 'color', [0.9 0.6 0.2], 'linewidth', 2.0)
        plot([1 N_trials], mean(eta_m_each, 'omitnan') + stdDev_eta_m*[1 1], ':', 'color', [0.9 0.6 0.2]) 
        plot([1 N_trials], mean(eta_m_each, 'omitnan') - stdDev_eta_m*[1 1], ':', 'color', [0.9 0.6 0.2])
        if ok_percentiles
            plot([1 N_trials], eta_m_16perc*[1 1], '--', 'color', [0.6 0.9 0.4], 'linewidth', 2.0)
            plot([1 N_trials], eta_m_84perc*[1 1], '--', 'color', [0.6 0.9 0.4], 'linewidth', 2.0)
        end
        xlabel('rank')
        ylabel('\eta_m (Pa.s.m)')
        legend('Viscosity values', 'mean', 'mean + 1 std.', 'mean - 1 std.', '16%', '84%')
        
        figure; histogram(log10(eta_m_each)); xlabel('log_{10}(\eta_m)')
        hold on
        ax = axis;
        plot(log10(eta_m*[1 1]), [ax(3) ax(4)], '-', 'color', [0.9 0.6 0.2], 'linewidth', 2.0) 
        
        figure; histogram(L_each); xlabel('L')
    end
else
    stdDev_eta_m = 0;
    unc_L = 0;
    eta_m_16perc = [];
    eta_m_84perc = [];
    eta_m_each = [];
    L_each = [];
end


%% Plot
if plotopt
    eta_list = (0.1:0.25:14)*1e-9;
    L_list = zeros(1,length(eta_list));
    for j = 1:length(eta_list)
        Dr_equation = D_R*fig11_fit(eta_w, L_eq, eta_list(j)) == k_B*Temperature/(4*pi*eta_list(j)*L_eq^2);
        S = vpasolve(Dr_equation, L_eq, L_seed);
        L_list(j) = S;
    end
    L_list2 = zeros(1,length(eta_list));
    for j = 1:length(eta_list)
        D_parl_equation = D_par*fig3_fit(eta_w, L_eq, eta_list(j)) == k_B*Temperature/(4*pi*eta_list(j));
        S = vpasolve(D_parl_equation, L_eq, L_seed);
        L_list2(j) = S;
    end

    figure('Name', 'L vs. eta_m');
    title('Uses $D_{\parallel}$ value', 'Interpreter', 'latex')
    xlabel('\eta_m (Pa\cdots\cdotm)')
    ylabel('L (m)')
    hold on;
    plot(eta_list, L_list)
    plot(eta_list, L_list2)
    legend('$D_R$', '$D_{\parallel}$', 'Interpreter', 'latex')
end

end


%% Fit functions of the Levine figures

function [c_parl] = fig3_fit(eta_w, L, eta_m)
ell_0 = eta_m / (2*eta_w); 
c_3 = 2.2842e-3;
c_2 = 4.8769e-2;
c_1 = 4.5546e-1;
c_0 = -7.0028e-1;
c_parl = exp(c_0 + c_1*log(L/ell_0) + c_2*log(L/ell_0).^2 + ...
    c_3*log(L/ell_0).^3);
end

function [c_perp] = fig4_fit(eta_w, L, eta_m)
ell_0 = eta_m / (2*eta_w); 
c_5 = -1.0835e-4;
c_4 = -4.9564e-4;
c_3 = 1.0684e-2;
c_2 = 1.1445e-1;
c_1 = 6.4181e-1;
c_0 = -7.9744e-2;
c_perp = exp(c_0 + c_1*log(L/ell_0) + c_2*log(L/ell_0).^2 + ...
    c_3*log(L/ell_0).^3 + c_4*log(L/ell_0).^4 + c_5*log(L/ell_0).^5);
end

function [c_R] = fig11_fit(eta_w, L, eta_m)
ell_0 = eta_m / (2*eta_w); 
c_5 = -4.4659e-5;
c_4 = -2.8110e-4;
c_3 = 5.4812e-3;
c_2 = 6.2408e-2;
c_1 = 2.1796e-1;
c_0 = -1.8594;
c_R = exp(c_0 + c_1*log(L/ell_0) + c_2*log(L/ell_0).^2 + ...
    c_3*log(L/ell_0).^3 + c_4*log(L/ell_0).^4 + c_5*log(L/ell_0).^5);
end
