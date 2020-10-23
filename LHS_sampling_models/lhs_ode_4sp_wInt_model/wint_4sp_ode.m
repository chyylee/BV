% Model structure 4: 4 speices model with interactions
% This model includes interspecies interactions for 4 species that can
% uptake MNZ.

function dy = wint_4sp_ode(t, y, params)

neq = size(y,1);
dy = zeros(neq,1);

% Names of state variables
LB1 = y(1);
GV1 = y(2);
MNZext = y(3);
MNZint_lb1 = y(4);
MNZint_gv1 = y(5);
MET1 = y(6);
LB2 = y(7);
GV2 = y(8);
MNZint_lb2 = y(9);
MNZint_gv2 = y(10);
MET2 = y(11);

% Main Parameters
bgv1 = params(1);
agv1 = params(2);
bingv1 = params(3);
kmgv1 = params(4);
Kgv1 = params(5);
EC50gv1 = params(6);
blb1 = params(7);
binlb1 = params(8);
Klb1 = params(9);
alb1 = params(10);
EC50lb1 = params(11);
bgv2 = params(12);
agv2 = params(13);
bingv2 = params(14);
kmgv2 = params(15);
Kgv2 = params(16);
EC50gv2 = params(17);
blb2 = params(18);
binlb2 = params(19);
Klb2 = params(20);
alb2 = params(21);
EC50lb2 = params(22);

% Additional Interaction Terms
fL1L2 = params(23); % L1 -> L2*
fL2L1 = params(24); % L2 -> L1*
fL1G1 = params(25); % L1 -> G1
fG1L1 = params(26); % G1 -> L1
fL2G1 = params(27); % L2 -> G1*
fG1L2 = params(28); % G1 -> L2
fL1G2 = params(29); % L1 -> G2
fG2L1 = params(30); % G2 -> L1
fL2G2 = params(31); % L2 -> G2*
fG2L2 = params(32); % G2 -> L2
fG1G2 = params(33); % G1 -> G2*
fG2G1 = params(34); % G2 -> G1*

% Relation to carrying capacity:
sL1L2 = (Klb2 - fL1L2*Klb2)/Klb1;
sL2L1 = (Klb1 - fL2L1*Klb1)/Klb2;
sL1G1 = (Kgv1 - fL1G1*Kgv1)/Klb1;
sG1L1 = (Klb1 - fG1L1*Klb1)/Kgv1;
sL2G1 = (Kgv1 - fL2G1*Kgv1)/Klb2;
sG1L2 = (Klb2 - fG1L2*Klb2)/Kgv1;
sL1G2 = (Kgv2 - fL1G2*Kgv2)/Klb1;
sG2L1 = (Klb1 - fG2L1*Klb1)/Kgv2;
sL2G2 = (Kgv2 - fL2G2*Kgv2)/Klb2;
sG2L2 = (Klb2 - fG2L2*Klb2)/Kgv2;
sG1G2 = (Kgv2 - fG1G2*Kgv2)/Kgv1;
sG2G1 = (Kgv1 - fG2G1*Kgv1)/Kgv2;

% Logistic Term for Each Species:
L1_INT = 1/Klb1*(LB1 + sL2L1*LB2 + sG1L1*GV1 + sG2L1*GV2);
G1_INT = 1/Kgv1*(GV1 + sL2G1*LB2 + sG2G1*GV2 + sL1G1*LB1);
L2_INT = 1/Klb2*(LB2 + sL1L2*LB1 + sG1L2*GV1 + sG2L2*GV2);
G2_INT = 1/Kgv2*(GV2 + sL2G2*LB2 + sG1G2*GV1 + sL1G2*LB1);

% Oridanry Differential Equations:
dy(1) =  blb1*LB1*(1 - L1_INT) - alb1* MNZint_lb1/(EC50lb1 + MNZint_lb1)*LB1; %LB1
dy(2) = bgv1*GV1*(1 - G1_INT) - agv1* MNZint_gv1/(EC50gv1 + MNZint_gv1)*GV1; %GV1
dy(3) = - binlb1*MNZext*LB1 - bingv1*MNZext*GV1 - binlb2*MNZext*LB2 - bingv2*MNZext;
dy(4) = binlb1*MNZext*LB1;
dy(5) = (bingv1*MNZext - kmgv1*MNZint_gv1)*GV1;
dy(6) = (kmgv1 * MNZint_gv1)*GV1;
dy(7) =  blb2*LB2*(1 - L2_INT) - alb2* MNZint_lb2/(EC50lb2 + MNZint_lb2)*LB2; %LB2
dy(8) = bgv2*GV2*(1 - G2_INT) - agv2* MNZint_gv2/(EC50gv2 + MNZint_gv2)*GV2; %GV2
dy(9) = binlb2*MNZext*LB2;
dy(10) = (bingv2*MNZext - kmgv2*MNZint_gv2)*GV2;
dy(11) = (kmgv2 * MNZint_gv2)*GV2;

end


