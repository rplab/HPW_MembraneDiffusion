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
%    varargin : variable input arguments 
%    -------------
%    D_par : Parallel diffusion coefficient (m^2/s), default 0.1e-12 m^2/s
%    unc_D_par : uncertainty in Parallel diffusion coefficient (m^2/s),  default 0.01e-12 m^2/s
%    D_perp : Parallel diffusion coefficient (m^2/s), default 0.1e-12 m^2/s
%    unc_D_perp : uncertainty in Parallel diffusion coefficient (m^2/s),  default 0.01e-12 m^2/s
%    D_R : Rotational diffusion coefficient (rad^2/s)m, default 0.1 
%    unc_D_R : uncertainty in Rotational diffusion coefficient (rad^2/s), default 0.01
%    N_trials : number of trials for repeated sampling, for uncertainty.
%       (Default 1)
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
%            If N_trials = 1, , "standard deviation" is just 0.5*difference
%              between D_par- and D_perp-derived values.)
%    L : Bead length (m)
%    unc_L : uncertainty in Bead length (m)
%            If N_trials = 1, , "standard deviation" is just 0.5*difference
%              between D_par- and D_perp-derived values.)
%    eta_m_16perc, eta_m_84perc : 16th and 84th percentile eta_m values,
%        from sampling for uncertainty estimates; empty if N_trials = 1.
%        Note that for a Gaussian distribution these are +/- 1 sigma, but
%        the eta_m distributions are non-Gaussian (see 22 July 2021 notes)
%    eta_m_each : all the eta_m values from repeated sampling
%    L_each : all the eta_m values from repeated sampling
%
% Based on Philip Jahl: Levine_solve_parl.m
% Raghuveer Parthasarathy:
% Revised syms to work on older MATLAB version
% Write LLM equations using "L/ell_0"
% Uncertainty calculation (sampling)!
%
% Raghuveer Parthasarathy
% June 30, 2021
% October 12, 2021: calculate (eta, L) from (Dperp, D_R) as well as 
%    (D_par, D_R) rather than only the latter; average the two together 
%    and combine uncertainties. See Oct 2021 notes.
%    Also, to avoid errors in inputs, use input parser.
% last modified October 13, 2021



function [eta_m, stdDev_eta_m, L, unc_L, eta_m_16perc, eta_m_84perc, eta_m_each, L_each] = ...
    LLM_fit_viscosity_L(varargin)

p = inputParser;
valid_positive = @(x) isnumeric(x) && isscalar(x) && (x > 0);
valid_D_T = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1e-6); % last check would fail if input as um^2/s, for example
valid_N = @(x) isnumeric(x) && isscalar(x) && (mod(x,1)<eps); % last check enforces integer

% Default values, and check validity
addOptional(p,'D_par', 0.1e-12, valid_D_T); % m^2/s
addOptional(p,'unc_D_par', 0.01e-12, valid_D_T); % m^2/s
addOptional(p,'D_perp', 0.1e-12, valid_D_T); % m^2/s
addOptional(p,'unc_D_perp', 0.01e-12, valid_D_T); % m^2/s
addOptional(p,'D_R', 0.1, valid_positive); % rad^2/s
addOptional(p,'unc_D_R', 0.01, valid_positive); % rad^2/s
addOptional(p,'N_trials', 1, valid_N); 
addOptional(p,'Temperature', 295, valid_positive); % K
addOptional(p,'eta_w', 8.9e-4, valid_positive); % water viscosity, Pa.s
   % Or use:     eta_w = 1.01e-3; % 0.1M sucrose viscosity, Pa.s
addOptional(p,'plotopt', true, @islogical);  
parse(p,varargin{:});

if p.Results.plotopt
    disp('All input parameters: ');
    p.Results
end

% Will assign some to to non-structured variables, to clarify
D_par = p.Results.D_par;
unc_D_par = p.Results.unc_D_par;
D_perp = p.Results.D_perp;
unc_D_perp = p.Results.unc_D_perp;
D_R = p.Results.D_R;
unc_D_R = p.Results.unc_D_R;

%% Other parameters

k_B = 1.380649e-23; % Boltzmann's constant
eta_seed = 5e-9;
L_seed = 1e-6;

syms eta_m_eq L_eq positive % "_eq" added for symbolic variables
assumeAlso(eta_m_eq, 'real')
assumeAlso(L_eq, 'real')

%% Diffusion coefficient equations

Dr_equation = D_R*fig11_fit(p.Results.eta_w, L_eq, eta_m_eq) == k_B*p.Results.Temperature/(4*pi*eta_m_eq*L_eq^2);
D_parl_equation = D_par*fig3_fit(p.Results.eta_w, L_eq, eta_m_eq) == k_B*p.Results.Temperature/(4*pi*eta_m_eq);
D_perp_equation = D_perp*fig4_fit(p.Results.eta_w, L_eq, eta_m_eq) == k_B*p.Results.Temperature/(4*pi*eta_m_eq);

% Calculate eta_m and L for (D_par, D_R), and then (Dperp, D_R)
% Note that these will be overwritten if we're sampling a distribution.
% (Redundant, but ok...)

% from D_parallel
[eta_m_fromDpar, L_fromDpar] = vpasolve([Dr_equation, D_parl_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
eta_m_fromDpar = double(eta_m_fromDpar); % Necessary because it's possible to have no solution
L_fromDpar = double(L_fromDpar);
if isempty(eta_m_fromDpar)
    eta_m_fromDpar = NaN;
    L_fromDpar = NaN;
end

% from D_perpendicular
[eta_m_fromDperp, L_fromDperp] = vpasolve([Dr_equation, D_perp_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
eta_m_fromDperp = double(eta_m_fromDperp); % Necessary because it's possible to have no solution
L_fromDperp = double(L_fromDperp);
if isempty(eta_m_fromDperp)
    eta_m_fromDperp = NaN;
    L_fromDperp = NaN;
end

% Combine (if both not NaN; otherwise just keep one)
[eta_m, stdDev_eta_m, L, unc_L] = combine_eta_L(eta_m_fromDpar, eta_m_fromDperp, ...
    L_fromDpar, L_fromDperp);

%% Repeated calculations for uncertainties
if p.Results.N_trials>1
    eta_m_each = zeros(p.Results.N_trials,1);
    L_each = zeros(p.Results.N_trials,1);
    rand_Dparl = D_par + unc_D_par*randn(p.Results.N_trials,1);
    rand_Dperp = D_perp + unc_D_perp*randn(p.Results.N_trials,1);
    rand_Dr = D_R + unc_D_R*randn(p.Results.N_trials,1);
    for j=1:p.Results.N_trials
        if mod(j,100)==0
            fprintf('trial %d of %d.\n', j, p.Results.N_trials);
        end
        this_D_parl = rand_Dparl(j);
        this_D_perp = rand_Dperp(j);
        this_Dr = rand_Dr(j);
        % I think I need to re-define the equations...
        Dr_equation = this_Dr*fig11_fit(p.Results.eta_w, L_eq, eta_m_eq) == ...
            k_B*p.Results.Temperature/(4*pi*eta_m_eq*L_eq^2);
        D_parl_equation = this_D_parl*fig3_fit(p.Results.eta_w, L_eq, eta_m_eq) == ...
            k_B*p.Results.Temperature/(4*pi*eta_m_eq);
        D_perp_equation = this_D_perp*fig4_fit(p.Results.eta_w, L_eq, eta_m_eq) == ...
            k_B*p.Results.Temperature/(4*pi*eta_m_eq);
        % Solve using D_par and D_perp, separately; then average if not NaN
        [sampled_eta_m_out_Dpar, sampled_L_out_Dpar] = ...
            vpasolve([Dr_equation, D_parl_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
        sampled_eta_m_out_Dpar = double(sampled_eta_m_out_Dpar); % Necessary because it's possible to have no solution
        sampled_L_out_Dpar = double(sampled_L_out_Dpar);
        if isempty(sampled_eta_m_out_Dpar)
            sampled_eta_m_out_Dpar = NaN;
            sampled_L_out_Dpar = NaN;
        end
        [sampled_eta_m_out_Dperp, sampled_L_out_Dperp] = ...
            vpasolve([Dr_equation, D_perp_equation], [eta_m_eq, L_eq], [eta_seed; L_seed]);
        sampled_eta_m_out_Dperp = double(sampled_eta_m_out_Dperp); % Necessary because it's possible to have no solution
        sampled_L_out_Dperp = double(sampled_L_out_Dperp);
        if isempty(sampled_eta_m_out_Dperp)
            sampled_eta_m_out_Dperp = NaN;
            sampled_L_out_Dperp = NaN;
        end
        % combine parallel and perpendicular. Ignore uncertainty since this
        % will be calculated by sampling
        [eta_m_each(j), ~, L_each(j), ~] = ...
            combine_eta_L(sampled_eta_m_out_Dpar, sampled_eta_m_out_Dperp, ...
            sampled_L_out_Dpar, sampled_L_out_Dperp);
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
    
    if p.Results.plotopt
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
        plot([1 p.Results.N_trials], mean(eta_m_each, 'omitnan')*[1 1], '-', 'color', [0.9 0.6 0.2], 'linewidth', 2.0)
        plot([1 p.Results.N_trials], mean(eta_m_each, 'omitnan') + stdDev_eta_m*[1 1], ':', 'color', [0.9 0.6 0.2]) 
        plot([1 p.Results.N_trials], mean(eta_m_each, 'omitnan') - stdDev_eta_m*[1 1], ':', 'color', [0.9 0.6 0.2])
        if ok_percentiles
            plot([1 p.Results.N_trials], eta_m_16perc*[1 1], '--', 'color', [0.6 0.9 0.4], 'linewidth', 2.0)
            plot([1 p.Results.N_trials], eta_m_84perc*[1 1], '--', 'color', [0.6 0.9 0.4], 'linewidth', 2.0)
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
    eta_m_16perc = [];
    eta_m_84perc = [];
    eta_m_each = [];
    L_each = [];
end


%% Plot
if p.Results.plotopt
    eta_list = (0.1:0.25:14)*1e-9;
    L_list = zeros(1,length(eta_list));
    for j = 1:length(eta_list)
        Dr_equation = D_R*fig11_fit(p.Results.eta_w, L_eq, eta_list(j)) == k_B*p.Results.Temperature/(4*pi*eta_list(j)*L_eq^2);
        S = vpasolve(Dr_equation, L_eq, L_seed);
        L_list(j) = S;
    end
    L_list2 = zeros(1,length(eta_list));
    for j = 1:length(eta_list)
        D_parl_equation = D_par*fig3_fit(p.Results.eta_w, L_eq, eta_list(j)) == k_B*p.Results.Temperature/(4*pi*eta_list(j));
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

function [eta, stdDev_eta, L, unc_L] = combine_eta_L(eta_1, eta_2, ...
    L_1, L_2)
% Average of each of the two eta, L values, omitting NaNs
eta = mean([eta_1 eta_2], 'omitnan');
L = mean([L_1 L_2], 'omitnan');
% uncertainties as 1/2 difference, or NaN if any are NaN (no info about
% uncertainty)
if sum(isnan([eta_1 eta_2]))==0
    % Both real numbers
    stdDev_eta = 0.5*abs(eta_1 - eta_2);
    unc_L = 0.5*abs(L_1 - L_2);
else
    stdDev_eta = NaN;
    unc_L = NaN;
end
end
