% Model structure 3: 4 species model with no interactions
% This model assumes no microbe-microbe interactions and simply that each
% bacteria is internalizing MNZ. There are 2 BV-associated bacteria, and 2
% LB spp.

function dy = Noint_4sp_ode(t, y, params)

neq = size(y,1);
dy = zeros(neq,1);

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
b3 = params(12);
a3 = params(13);
bin3 = params(14);
km3 = params(15);
K3 = params(16);
EC503 = params(17);
b4 = params(18);
bin4 = params(19);
K4 = params(20);
a4 = params(21);
EC504 = params(22);

dy(1) =  b2*LB1*(1 - LB1/K2) - a2* MNZint_lb1/(EC502 + MNZint_lb1)*LB1; %LB1
dy(2) = b1*GV1*(1 - GV1/K1) - a1* MNZint_gv1/(EC501 + MNZint_gv1)*GV1; %GV1
dy(3) = - bin2*MNZext*LB1 - bin1*MNZext*GV1 - bin4*MNZext*LB2 - bin3*MNZext;
dy(4) = bin2*MNZext*LB1;
dy(5) = (bin1*MNZext - km1*MNZint_gv1)*GV1;
dy(6) = (km1 * MNZint_gv1)*GV1;
dy(7) =  b4*LB2*(1 - LB2/K4) - a4* MNZint_lb2/(EC504 + MNZint_lb2)*LB2; %LB2
dy(8) = b3*GV2*(1 - GV2/K3) - a3* MNZint_gv2/(EC503 + MNZint_gv2)*GV2; %GV2
dy(9) = bin4*MNZext*LB2;
dy(10) = (bin3*MNZext - km3*MNZint_gv2)*GV2;
dy(11) = (km3 * MNZint_gv2)*GV2;

end


