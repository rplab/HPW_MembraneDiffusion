% HPW_drag.m
% 
% Function to calculate translational and rotational diffusion coefficients 
% as a function of membrane viscosity (eta) and inclusion radius (a)
% using the exact equations of
%    B. D. Hughes, B. A. Pailthorpe, White L. R.,
%    The translational and rotational drag on a cylinder moving in a membrane. 
%    J Fluid Mech. 110, 349-372 (1981)
% 
% Also allows for a polynomial approximation to the Lambda(epsilon) factor
% in the translational drag, to account for numerical divergences around 
% epsilon = 3.45 and larger in the series calculations. (Recommended if 
% considering epsilon in this range.)
%
% Based on earlier notes, programs. (e.g. June 2013, calcDHPW.m) These made
% use of numerical integration for various parts of the calculation, rather
% than the faster and more elegant series expansions. I figured out the
% error in HPW Equation A13, namely ln(2e) should be ln(2/e), which makes
% the series approach work!
%
% Note that the "reduced drag" only depends on epsilon = 2.*a*eta_w./eta, but
% I've written eta and a separately as inputs to make fitting to
% "physical" diffusion coefficients more straightforward.
%
% Also see https://eighteenthelephant.com/2018/09/04/membrane-diffusion-software/
%
% Inputs:
%    eta = 2D membrane viscosity (Pa.s.m); must be a single (scalar) value
%    a   = inclusion radius (m); must be a single (scalar) value
%    eta_w = external (water) viscosity, default 8.9e-4 Pa.s
%    T = temperature (K), default 295;
%    nTerms =  number of terms to calculate in "infinite" series; default 36
%    num_lterms = number of "l" terms to calculate in ; default 36
%    num_mterms = number of "m" terms to calculate in ; default 36;
%    plotopt = true to plot various things including tm, Tm, rm terms (default false)
%    Lambda_T_option: calculation option for translational drag 
%               'full' [default], HPW series calculation of Lambda_T
%                    (Eq. 3.68 of Hughes et al. 1984)
%               'simplified' : Lambda(T) = 1 / (epsilon T) 
%               'polynomial' : Lambda(T) = 1 / (epsilon T) *
%                              polynomial approximation to remaining
%                              factors (a function of epsilon). This avoids
%                              divergences in numerical calculations of
%                              series for epsilon around 3.45, etc.
%
% Outputs:
%    D_T = Translational diffusion coefficient (m^2/s)
%    D_R = Rotational diffusion coefficient (rad^2/s)
%    lambda_T = Translational drag coefficient (N s /m)
%    lambda_R = Rotational drag coefficient (N s m)
%
%%
% Copyright 2018-2021, Raghuveer Parthasarathy, The University of Oregon
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
% Raghuveer Parthasarathy
% August 28, 2018
% Last modified: August 20, 2021

function [D_T, D_R, lambda_T, lambda_R] = HPW_drag(eta, a, eta_w, ...
    Temperature, nTerms, num_lterms, num_mterms, plotopt, lambda_T_option)


%% Defaults
if ~exist('eta_w', 'var') || isempty(eta_w)
    eta_w = 8.9e-4; % water viscosity, Pa.s
end
if ~exist('Temperature', 'var') || isempty(Temperature)
    Temperature = 295; % K
end
if ~exist('nTerms', 'var') || isempty(nTerms)
    nTerms = 36;
end
if ~exist('num_lterms', 'var') || isempty(num_lterms)
    num_lterms = 36;
end
if ~exist('num_mterms', 'var') || isempty(num_mterms)
    num_mterms = 36;
end
if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end


%% Other parameters

k_B = 1.3806e-23; % Boltzmann's constant, J/K
epsilon = 2.*a*eta_w/eta; % dimensionless

%% Figures for tm, rm plots

if plotopt
    % Create figures for "tm" and "rm" plots, and colormaps
    h_tm = figure('name', 'tm', 'position', [50 50 800 400]); 
    h_rm = figure('name', 'rm', 'position', [100 100 800 400]); 
    r_t_colormap = parula(num_mterms+1);
else
    h_tm = [];
    h_rm = [];
    r_t_colormap = [];
end

%% Calculations...

% "t" terms
% Will calculate each tm separately. Inefficient -- I could probably
% vectorize this -- but I'm not going to bother. It would take some
% thought, since the "t1," "t2," and "t3" arrays have different sizes
tm_array = zeros(1, num_mterms);
for m=1:num_mterms
    tm_array(m) = calc_tm(m, epsilon, nTerms, h_tm, r_t_colormap);
end

% "T" terms
T = calc_T(epsilon, nTerms, plotopt);

% "r" terms
rm_array = zeros(1, num_mterms);
for m=1:num_mterms
    rm_array(m) = calc_rm(m, epsilon, nTerms, h_rm, r_t_colormap);
end


%% T'_l, t'_lm, and t''_lm
% Eq. 3.63, 3.60, 3.61
% These depend only on l, m; not on epsilon

[Tp_l, tp_lm, tpp_lm] = calc_Tprime_terms(num_lterms, num_mterms);

%% Xm1 and Xm2

% Solve for Xm1, using Eq. 3.67 (1)
M1 = tp_lm - epsilon*tpp_lm;
% solution to Ax = B is x = A\B, or x = inv(A)*B
%Xm1 = M1 \ Tp_l;
% Solve for Xm2, using Eq. 3.67 (2)
delta_l1 = double((1:num_lterms)'==1);
%Xm2 = M1 \ delta_l1;

% Is this any better? -- I think so...
Xm1 = lsqminnorm(M1,Tp_l);
Xm2 = lsqminnorm(M1,delta_l1);

if plotopt
    % Xm1 and Xm2
    figure('name', 'Xm1 and Xm2');
    plot(1:num_mterms, Xm1, 'd-')
    hold on
    plot(1:num_mterms, Xm2, 'o-')
    xlabel('m');
    ylabel('Xm')
    legend('Xm1', 'Xm2')
    
    % Products, used in later sums
    figure('name', 't and Xm1, Xm2 products');
    subplot(1,2,1)
    plot(1:num_mterms, tm_array'.*Xm1, 'd-')
    hold on
    plot(1:num_mterms, tm_array'.*Xm2, 's-')
    xlabel('m');
    ylabel('Products')
    legend('tm.*Xm1', 'tm.*Xm2')
    subplot(1,2,2)
    plot(1:num_mterms, cumsum(tm_array'.*Xm1), 'd-')
    hold on
    plot(1:num_mterms, cumsum(tm_array'.*Xm2), 's-')
    xlabel('m');
    ylabel('Cumulative Sum')
    legend('\Sigma tm.*Xm1', '\Sigma tm.*Xm2')
end

%% R_lm, R'_lm, Phim1, Phim2

[R_lm, Rp_lm] = calcR_lm(num_lterms, num_mterms);
% Note the sums begin at m==1, not m==0
% Calculate R_lm and Rp_lm including m==0, but don't use those rows
M2 = R_lm(2:end,2:end) - epsilon*Rp_lm(2:end,2:end); % start at l==1, m==1 (Eq. 5.35)
% Solve for Phim1, using Eq. 5.35 (1)
Phim1 = M2 \ (-R_lm(2:end,1) + epsilon*Rp_lm(2:end,1)); % "l" >= 1 , and m==0
% Solve for Phim2, using Eq. 5.35 (2)
delta_l1 = (1:num_lterms)'==1;
Phim2 = M2 \ delta_l1;

if plotopt
    % Phim1 and Phim2
    figure('name', 'Phim1 and Phim2');
    plot(1:num_mterms, Phim1, 'd-')
    hold on
    plot(1:num_mterms, Phim2, 'o-')
    xlabel('m');
    ylabel('Phim')
    legend('Phim1', 'Phim2')
    
    % Products, used in  later sums
    figure('name', 'r and Phim1, Phim2 products');
    subplot(1,2,1)
    plot(1:num_mterms, rm_array'.*Phim1, 'd-')
    hold on
    plot(1:num_mterms, rm_array'.*Phim2, 's-')
    xlabel('m');
    ylabel('Products')
    legend('rm.*Phim1', 'rm.*Phim2')
    subplot(1,2,2)
    plot(1:num_mterms, cumsum(rm_array'.*Phim1), 'd-')
    hold on
    plot(1:num_mterms, cumsum(rm_array'.*Phim2), 's-')
    xlabel('m');
    ylabel('Cumulative Sum')
    legend('\Sigma rm.*Phim1', '\Sigma rm.*Phim2')
end



%% Drag coefficients

%% Drag coefficients

% Translation
% HPW Eq. 3.68:
% Note that the tm can be undefined (numerical issues?) for large m, so 
% restrict the sum to finite values.
% Original:  Lambda_T = (1 + (epsilon^2)/15.*sum(tm_array'.*Xm2))./epsilon./(T + epsilon*sum(tm_array'.*Xm1));
tm_Xm1_array = tm_array'.*Xm1; % Element-wise product of the two arrays
tm_Xm2_array = tm_array'.*Xm2;
if strcmp(lambda_T_option, 'full')
    Lambda_T = (1 + (epsilon^2)/15 .* ...
        sum(tm_Xm2_array(isfinite(tm_Xm2_array))))./epsilon ./ ...
        (T + epsilon*sum(tm_Xm1_array(isfinite(tm_Xm1_array))));
elseif strcmp(lambda_T_option, 'simplified')
    Lambda_T = 1./(epsilon*T);
elseif strcmp(lambda_T_option, 'polynomial')
    Lambda_T_coeffs = [-2.33337e-03  -1.37628e-02  9.77674e-01]; % from fit, Jan. 2021
    Lambda_T = 1./(epsilon*T).*(Lambda_T_coeffs(1)*(log(epsilon)).^2 + ...
        Lambda_T_coeffs(2)*log(epsilon) + Lambda_T_coeffs(3) );
end
    
% Eq. 3.43:
lambda_T = 8*pi*eta_w*a*Lambda_T; % drag coefficient
D_T = k_B * Temperature ./ lambda_T;

% Rotation
r0 = calc_rm(0, epsilon, nTerms, h_rm, r_t_colormap);
% Note that the rm can be undefined (numerical issues?) for large m, so 
% restrict the sum to finite values.
% Original: Lambda_R = (2/9./epsilon).*(1 + (epsilon.^2)/35.*sum(rm_array'.*Phim2))./(r0 + sum(rm_array'.*Phim1));
rm_Phim1_array = rm_array'.*Phim1; % Element-wise product of the two arrays
rm_Phim2_array = rm_array'.*Phim2; 
Lambda_R = (2/9./epsilon).*(1 + (epsilon.^2)/35 .* ...
    sum(rm_Phim2_array(isfinite(rm_Phim2_array)))) ./ ...
    (r0 + sum(rm_Phim1_array(isfinite(rm_Phim1_array))));

% Eq. 5.16
lambda_R = 8*pi*eta_w*a.^3*Lambda_R;
D_R = k_B * Temperature ./ lambda_R;

end

function tm = calc_tm(m, epsilon, nTerms, h_tm, r_t_colormap)
% Calculates "tm" for a given m; HPW Eq. A5-A8
    % each "t" term
    % t1: HPW Eq. A6
    %n1 = 0:nTerms;
    %t1_terms = (-1).^(m-n1).*gamma(2*n1+3).*(epsilon/2).^(2*n1+1)./gamma(n1+2-m)./gamma(n1+2+m)./gamma(n1+1.5-m)./gamma(n1+2.5+m);
    
    % Restrict n to start at m-1, rather than 0, so that the first gamma in
    % the denominator avoids negative integer arguments (undefined!) 
    % Or is Eq. A6 incorrect?
    % -- RP 22Dec 2020
    n1 = (m-1):nTerms;
    allgammas = exp(gammaln(2*n1+3) - gammaln(n1+2-m) - gammaln(n1+2+m) - gammaln(n1+1.5-m) - gammaln(n1+2.5+m));
    t1_terms = (-1).^(m-n1).*allgammas.*(epsilon/2).^(2*n1+1);
    t1 = pi*pi/8*sum(t1_terms);

%   nTerms
%    figure; semilogy(0:nTerms, t1_terms, 'x-')
%    hold on
%    semilogy(0:nTerms, -t1_terms, 'x-')
%figure; plot(n1, t1_terms, 'x-')
%pause

 % *** To do: Note the above, considering positive arguments only for the first denominator gamma! (i.e. start n1 at m-1)

    % t2: HPW Eq. A7
    if m >= 1
        n2=0:(m-1);
        allgammas = exp(gammaln(2*n2+2) + gammaln(m-n2) + gammaln(-0.5+m-n2)- gammaln(1.5+m+n2) - gammaln(2+m+n2));
        t2_terms = allgammas.*sin(pi*(m - n2 - 0.5))./pi.*(epsilon/2).^(2*n2);
        % t2_terms = gamma(2*n2+2).*gamma(m-n2).*(epsilon/2).^(2*n2)./gamma(1.5-m+n2)./gamma(1.5+m+n2)./gamma(2+m+n2);
    else
        t2_terms = 0;
    end
    t2 = pi/8*sum(t2_terms);

   % t3: HPW Eq. A8
    n3 = 0:nTerms;
%    t3_terms = (-1).^n3.*(log(2./epsilon) - psi(2+2*m+2*n3) + 0.5*psi(n3+1) + 0.5*psi(1.5+n3) + 0.5*psi(1.5+2*m+n3) + 0.5*psi(2+2*m+n3)) .* ...
%        gamma(2+2*m+2*n3).*(epsilon/2).^(2*m+2*n3)./gamma(n3+1)./gamma(1.5+n3)./gamma(1.5+2*m+n3)./gamma(2+2*m+n3);
    
    log_prodf = zeros(1,length(n3)); % note that the first element is zero, 
        %as it should be, since it corresponds to 0!
    all_f = (2*m+1) + (2:2*nTerms);
    log_all_f = log(all_f);
    cs_log_all_f = cumsum(log_all_f);
    log_prodf(2) = log_all_f(1);
    for j=3:length(n3)
        log_prodf(j) = cs_log_all_f(2*j-3) - cs_log_all_f(j-2);
    end

% Slower:
%    for j=2:length(n3)
%        % f = 2*m + 2 + n3(j) : 2*m + 1 + 2*n3(j);
%        f = (2*m + 1 + n3(j)) + (1 : n3(j));
%        log_prodf(j) = sum(log(f));
%    end
    log_other_gammas = -gammaln(n3+1) - gammaln(1.5+n3) - gammaln(1.5+2*m+n3);
    psi_factors = -psi(2+2*m+2*n3) + 0.5*psi(n3+1) + 0.5*psi(1.5+n3) + 0.5*psi(1.5+2*m+n3) + 0.5*psi(2+2*m+n3); % not faster; just for clarity
    t3_terms = (-1).^n3.*(log(2./epsilon) + psi_factors) .* ...
               exp(log_prodf + log_other_gammas).*(epsilon/2).^(2*m+2*n3);
    % t3_terms(isinf(t3_terms) | isnan(t3_terms))=0;  

    t3 = pi/4*sum(t3_terms);

    % HPW Eq. A5
    tm = t1 + t2 + t3;
    
    if ~isempty(h_tm)
        % Make plots, t1, t2, t3 as separate subplots
        figure(h_tm); 
        subplot(1,3,1); hold on
        plot(n1, t1_terms, 'o-', 'color', r_t_colormap(m,:))
        xlabel('n')
        ylabel('t terms')
        title('t1')
        subplot(1,3,2); hold on
        plot(n2, t2_terms, 'x-', 'color', r_t_colormap(m,:))
        xlabel('n')
        ylabel('t terms')
        title('t2')
        subplot(1,3,3); hold on
        plot(n3, t3_terms, '+-', 'color', r_t_colormap(m,:))
        xlabel('n')
        ylabel('t terms')
        title('t3')
        % legend('t1', 't2', 't3');
    end
end

function [Tp_l, tp_lm, tpp_lm] = calc_Tprime_terms(num_lterms, num_mterms)
% Calculate T', t'_lm, t''_lm;  Eq. 3.63, 3.60, 3.61
% These depend only on l, m; not on epsilon

    ell = (1:num_lterms)'; % I think l starts at 1, based on Eq. 3.58
    % (Or at least, the l==0 term is irrelevant for Lambda_T)

    % Tp_l: column vector
    % HPW Eq 3.63
    Tp_l = (3*((-1).^(ell+1))/32).*(gamma(ell-0.5)./gamma(ell+2)).^2;

    %tp_lm and tpp_lm
    m = 1:num_mterms;
    ell_array = repmat(ell, 1, num_mterms);
    m_array = repmat(m, num_lterms, 1);
    % HPW Eq 3.60:
    tp_lm = (-1).^(ell_array+m_array+1)/2 ./...
        (4*(ell_array-m_array.^2)-1)./(ell_array+m_array)./(ell_array+m_array+1);
    % HPW Eq 3.61:
    tpp_lm_term1 = ((ell_array-1)==m_array)./(4*ell_array-3)./(4*ell_array-1)./(4*ell_array+1);
    tpp_lm_term2 = 2*(ell_array==m_array)./(4*ell_array-1)./(4*ell_array+1)./(4*ell_array+3);
    tpp_lm_term3 = ((ell_array+1)==m_array)./(4*ell_array+1)./(4*ell_array+3)./(4*ell_array+5);
    tpp_lm = pi/2*(tpp_lm_term1 + tpp_lm_term2 + tpp_lm_term3);
end

function T = calc_T(epsilon, nTerms, plotopt)
% Calculates "T" terms: HPW Eq. A11-A13
% With an imoprtant correction: 2epsilon -> epsilon/2  in the T2 equation!!
    % HPW Eq. A12
    N1 = 0:nTerms;
    %T1_terms = (epsilon/2).^(2*N1+1).*(-1).^N1.*gamma(2.5+2*N1)./(gamma(N1+1.5).^2)./(gamma(N1+2).^2);
    allgammas = exp(gammaln(2.5+2*N1) - 2*gammaln(N1+1.5) - 2*gammaln(N1+2));
    T1_terms = (epsilon/2).^(2*N1+1).*(-1).^N1.*allgammas;
    T1 = (pi^1.5)/4*sum(T1_terms);

    % HPW Eq. A13
    N2 = 0:nTerms;
%     T2_terms = ((-1).^(N2)).*(epsilon/2).^(2*N2).*gamma(2*N2+1.5) .*...
%         (log(2./epsilon)-psi(2*N2+1.5)+psi(N2+1)+psi(N2+1.5)) ./...
%         (gamma(N2+1).^2)./(gamma(N2+1.5).^2);
    allgammas = exp(gammaln(2*N2+1.5) - 2*gammaln(N2+1) - 2*gammaln(N2+1.5));
    T2_terms = ((-1).^(N2)).*(epsilon/2).^(2*N2).*allgammas .*...
        (log(2./epsilon)-psi(2*N2+1.5)+psi(N2+1)+psi(N2+1.5));
    T2 = (pi^0.5)/2*sum(T2_terms);
    
    % Eq. A11
    T = T1 + T2;
    
    if plotopt
        figure; plot(N1, T1_terms, 'o-')
        hold on;
        plot(N2, T2_terms, 'x-')
        plot(N2, cumsum(T1_terms+T2_terms), 's:')
        xlabel('n')
        ylabel('T terms')
        legend('T1', 'T2', 'cumulative sum');
    end
end

function [R_lm, Rp_lm] = calcR_lm(num_lterms, num_mterms)
% HPW Eq. 5.32, 5.33
    ell = (0:num_lterms)'; % I think l starts at 0, based on Eq. 5.32, 5.33
    m = 0:num_mterms;
    ell_array = repmat(ell, 1, num_mterms+1);
    m_array = repmat(m, num_lterms+1, 1);
    R_lm = (-1).^(ell_array+m_array+1)/2./(4*(m_array - ell_array).^2 - 1)./(ell_array+m_array+1)./(ell_array+m_array+2);
    delta_l_mm1 = ell_array==(m_array-1);
    delta_l_m = ell_array==m_array;
    delta_l_mp1 = ell_array==(m_array+1);
    Rp_lm = pi*delta_l_mm1/2./(4*m_array+3)./(4*m_array+1)./(4*m_array-1) + ...
            pi*delta_l_m./(4*m_array+5)./(4*m_array+3)./(4*m_array+1) + ...
            pi*delta_l_mp1/2./(4*m_array+7)./(4*m_array+5)./(4*m_array+3);
end

function rm = calc_rm(m, epsilon, nTerms, h_rm, r_t_colormap)
% Calculates "rm" for a given m; HPW Eq. A14-A16
    % HPW Eq. A14
    % m==0 is irrelevant to the sum in Eq. 5.36 (which starts at m=1), but
    % it's necessary for the r_0 term in the denominator
    if m==0
        nr1 = 1:nTerms;
    else
        nr1 = 0:nTerms;
    end
    r1_terms = (epsilon/2).^(2*m+2*nr1-1).*gamma(1.5-nr1).*gamma(1+2*nr1+2*m)./gamma(1+2*m+nr1)./gamma(1+nr1)./gamma(2*m+nr1+2.5);
    r1 = -pi/8*sum(r1_terms);

    % HPW Eq. A15
    nr2 = 0:m;
    r2_terms = (epsilon/2).^(2*nr2).*gamma(1+m-nr2).*gamma(2+2*nr2)./gamma(1.5+m+nr2)./gamma(1.5-m+nr2)./gamma(3+m+nr2);
    r2 = pi/8*sum(r2_terms);

    % HPW Eq. A16
    nr3 = 0:nTerms;
    r3_terms = (-1).^nr3.*(epsilon/2).^(2*m+2*nr3+2).*gamma(4+2*m+2*nr3).*...
        (log(2./epsilon) - psi(4+2*m+2*nr3) + 0.5*psi(2.5+2*m+nr3) + 0.5*psi(1+nr3) + 0.5*psi(nr3+2.5) + 0.5*psi(4+2*m+nr3)) ./ ...
        gamma(nr3+1)./gamma(2.5+2*m+nr3)./gamma(2.5+nr3)./gamma(4+2*m+nr3);
    r3 = pi/4*sum(r3_terms);

    % HPW Eq. A14
    rm = r1 + r2 + r3;
    
    if ~isempty(h_rm)
        % Make plots, r1, r2, r3 as separate subplots
        % note that m==0 is possible, so colormap is offset by 1
        figure(h_rm); 
        subplot(1,3,1); hold on
        semilogy(nr1, r1_terms, 'x-', 'color', r_t_colormap(m+1,:))
        xlabel('n')
        ylabel('r terms')
        title('r1')
        subplot(1,3,2); hold on
        semilogy(nr2, r2_terms, 'x-', 'color', r_t_colormap(m+1,:))
        xlabel('n')
        ylabel('r terms')
        title('r2')
        subplot(1,3,3); hold on
        semilogy(nr3, r3_terms, 'x-', 'color', r_t_colormap(m+1,:))
        xlabel('n')
        ylabel('r terms')
        title('r3')
        % legend('r1', 'r2', 'r3');
    end
end