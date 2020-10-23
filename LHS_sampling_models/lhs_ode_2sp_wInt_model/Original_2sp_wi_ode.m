% Model structure 2: 2-species with interactions model
% This model assumes microbe-microbe interactions between Gv and a LB sp.

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


