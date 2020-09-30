% (1) Bacterial-mediated metabolism of internalized MNZ, preventing MNZint
% from directly killing GV

% Experimental Conditions
% (1) MNZ only
% (2) MNZ + GV sup
% (3) MNZ + LI sup
% (4) MNZ + GV
% (5) MNZ + GV

% Response Conditions (collected at 0, 4, 6, 24, 48 hrs)
% (1) MNZ external
% (2) MNZ internal
% (3) Acetamide
% (4) Other Metabolites
% (5) Cell count

function dydt = GvLi_ODE_function(t,y,params)
[beta_GV, alpha_GV, beta_in, kmet, K, EC50, ...
    beta_lb, beta_in_lb,K_lb,alpha_LB,EC50_LB] = params{:};

% Defining parameters
LB = y(1);
GV = y(2);
MNZext = y(3);
MNZint_lb = y(4);
MNZint_gv = y(5);
MET = y(6);

L_li = alpha_LB* MNZint_lb/(EC50_LB + MNZint_lb);
L_gv = alpha_GV* MNZint_gv/(EC50 + MNZint_gv);
cell_vol = 1E-8;
cult_vol = 0.5;
pc_Li = MNZint_lb/(LB*cult_vol);
pc_Gv = MNZint_gv/(GV*cult_vol);

% % Ordinary Differential Equations
dydt=zeros(6,1);
dydt(1) = beta_lb*LB*(1 - LB/K_lb) - L_li*LB; % Li
dydt(2) = beta_GV*GV*(1 - GV/K) - L_gv*GV; % Gv
dydt(3) = - beta_in_lb*MNZext*LB - beta_in*MNZext*GV + ...
    L_li*LB*pc_Li*cell_vol + L_gv*GV*pc_Gv*cell_vol; %MNZ external
dydt(4) = beta_in_lb*MNZext*LB - L_li*LB*pc_Li*cell_vol; % Internal MNZ for LB 
dydt(5) = (beta_in*MNZext - kmet*MNZint_gv)*GV - L_gv*GV*pc_Gv*cell_vol; % Internal MNZ for GV - 
dydt(6) = (kmet * MNZint_gv)*GV; % Metabolites