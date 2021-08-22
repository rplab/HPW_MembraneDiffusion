% LLM_calcD.m
% 
% Function to calculate parallel, perpendicular, and rotational diffusion
% coefficients as a function of membrane viscosity (eta) and rod length (L)
% using the theory of 
%    A. J. Levine, T. B. Liverpool, F. C. MacKintosh, 
%    Mobility of extended bodies in viscous films and membranes. 
%    Phys. Rev. E. 69, 021503 (2004)
% Closed-form solutions are not available, so this uses polynomial
%    interpolation from the graphs of their paper, for an aspect ratio of 3.
% Write the interpolation functions here, redundantly with
%    LLM_fit_viscosity_L.m
% See notes, 2021
%
% Inputs:
%    eta = 2D membrane viscosity (Pa.s.m)
%    L   = inclusion major axis length (m)
%    eta_w = external (water) viscosity, default 8.9e-4 Pa.s
%    T = temperature (K), default 295;
%
% Outputs:
%    D_par = Translational diffusion coefficient along rod axis (m^2/s)
%    D_perp = Translational diffusion coefficient perpendicular to rod axis (m^2/s)
%    D_R = Rotational diffusion coefficient (rad^2/s)
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
% Raghuveer Parthasarathy
% August 15, 2021
%
% Last modified: August 15, 2021 

function [D_par, D_perp, D_R] = LLM_calcD(eta, L, eta_w, Temperature)

k_B = 1.380648e-23; % Boltzmann's constant;

kTn = k_B*Temperature/(4*pi*eta); % for convenience
D_par = kTn/fig3_fit(eta_w, L, eta);
D_perp = kTn/fig4_fit(eta_w, L, eta);
D_R = kTn/(L^2)/fig11_fit(eta_w, L, eta);

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
        
