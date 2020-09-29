% The ODEs for a simple predator/prey model with 2 equations and 4 parameters.
% This model is from page 6, section 3.1 of "A Methodology For Performing
% Global Uncertainty And Sensitivity Analysis In Systems Biology", Marino, et.
% al., Journal of Theoretical Biology, 2008-09-07, doi:
% 10.1016/j.jtbi.2008.04.011.

% From the above cited paper, with y(1) = Q and y(2) = P.

% It has two state variables (Q, P) and several parameters. Q represents the
% density of prey, P represents the density of predators, α is the intrinsic
% rate of prey population increase, β is the predation rate coefficient, σ is
% the predator mortality rate, and δ is the reproduction rate of predators per
% prey consumed.

% This file can be used as a template for creating a new ODE model. To do so
% copy this file and edit it, replacing the predator/prey parameters and
% equations with those for the new model.
%
% When creating a new model, you will also need to copy and edit file
% lhs_ode_predator_prey_settings_new.m to create a settings file for the new
% model, to define the model parameters values, initial conditions, etc.

function dy = Original_2sp_wi_ode(t, y, params)

neq = size(y,1);
dy = zeros(neq,1);

LB = y(1);
GV = y(2);
MNZext = y(3);
MNZint_lb = y(4);
MNZint_gv = y(5);
MET = y(6);


b1 = params(1);
a1 = params(2);
bin1 = params(3);
km1 = params(4);
K1 = params(5);
EC501 = params(6);
b2 = params(7);
bin2 = params(8);
K2 = params(9);
a2 = params(10);
EC502 = params(11);
sgl = params(12);
slg = params(13);

dy(1) =  b2*LB*(1 - (LB + sgl*GV)/K2) - a2* MNZint_lb/(EC502 + MNZint_lb)*LB;
dy(2) = b1*GV*(1 - (GV+slg*LB)/K1) - a1* MNZint_gv/(EC501 + MNZint_gv)*GV;
dy(3) = - bin2*MNZext*LB - bin1*MNZext*GV;
dy(4) = bin2*MNZext*LB;
dy(5) = (bin1*MNZext - km1*MNZint_gv)*GV;
dy(6) = (km1 * MNZint_gv)*GV;


end


