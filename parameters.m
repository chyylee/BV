function parameters()
clear;
beta_GV = 0.2269;
alpha_GV = 1.004;
beta_in = 0.01393;
kmet = 0.01741;
K = 4.2;
EC50 = 420;
beta_lb = 0.2309;
beta_in_lb = 0.00415;
K_lb = 3.569;
alpha_LB = 1.049;
EC50_LB = 598.87;

params = {beta_GV, alpha_GV, beta_in, kmet, K, EC50, ...
    beta_lb, beta_in_lb,K_lb,alpha_LB,EC50_LB};

paramnames = ["\beta _G_v", "\alpha _G_v", "\beta _{in Gv}", "kmet", ...
    "K_G_v", "EC50 _G_v", "\beta _L_i", "\beta _{in Li}", "K_L_i", ...
    "\alpha _L_i", "EC50 _L_i"];


save('params.mat','params', 'paramnames');
end
 


