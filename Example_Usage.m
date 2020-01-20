%% Examples of Code Usage
% January 20th, 2020

% This code will go through a few examples of how the ODE model was used to
% make figures in the manuscript. All of these functions call the same ODE
% file "GvLi_ODE_function.m", and require the parameter information in
% "parameters.m".

% The functions that can be called are:
% surface_data: Used for Fig. 3a,b and Fig. 4a
% Sense_1D_params: Used for Fig. 2a,b,d and fig. S3
% Sense_1D_ratio: Used for Fig. 2c

% Information on how to call ODE functions can be found on the MathWorks
% website: https://www.mathworks.com/help/matlab/ordinary-differential-equations.html?s_tid=CRUX_lftnav

% Specific questions on code:
% Christina Y. Lee at chyylee@umich.edu

%% 1) Generate a surface data (Fig. 3a,b and Fig. 4a), call surface_data.m
clear;
% inputs to the function
t = 48; % Time in hours
load('params.mat') % Loads parameter information
GvLi_range = logspace(-4.2,4.2,33); % Range of GvLi Ratios
MNZ_range = logspace(1,3.2,33); % Range of MNZ doses

[surface_out, vpGvLi, vpMNZ] = surface_data(t,params,GvLi_range,MNZ_range);

%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LI = squeeze(surface_out(1,:,:));
GV = squeeze(surface_out(2,:,:));

[X Y] = meshgrid(vpGvLi, vpMNZ);

figure(1)
subplot(1,2,1)
s = surf(X,Y,GV,'FaceAlpha',0);
zlabel('% Max Growth G. vaginalis')
xlabel('Gv:Li ratio')
ylabel('[MNZ]')
set(gca, 'YScale', 'log','XScale','log')
set(gca,'Fontsize',14)

subplot(1,2,2)
s = surf(X,Y,LI,'FaceAlpha',0);
zlabel('% Max Growth L. iners')
xlabel('Gv:Li ratio')
ylabel('[MNZ]')
set(gca, 'YScale', 'log','XScale','log')
set(gca,'Fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2) Parameter Sensitivity Analysis (Fig. 2a,b,d and fig. S3)
clear;
% inputs to the function:
t = 48; % time (hrs)
dParams = [1/1000 1/100 1/50 1/10 1/5 1 5 10 50 100 1000]; % Fold Change in baseline parameter
MNZ = 500; % MNZ dose (ug/ul)
load('params.mat')% call parameters
GvLi_ratio = 1; % Fold change of Gv compare to Li (here 1:1)

[gv_out, li_out, output_data] = Sense_1D_Params(t,dParams,MNZ,params,GvLi_ratio);

%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = figure(2);
subplot(2,2,1)
gv_f = [gv_out(1:4,:); gv_out(6:8,:); gv_out(10:11,:)]*100;
li_f = [li_out(1:4,:); li_out(6:8,:); li_out(10:11,:)]*100;
pLabels = {'kgrow_G_V', 'kkill_G_V', 'kint_G_V', 'kmet', 'EC50_G_V'};
cmap1  = colormap(ax2, parula(11));
h1 = bar(gv_f(1:5,:), 'FaceColor', 'flat');
for c = 1:size(squeeze(output_data(2,:,:)),2)
    h1(c).CData = cmap1(c,:);
end
set(gca, 'XTick',1:length(pLabels),'XTickLabel', pLabels)
legend('0.001x', '0.01x', '0.02x', '0.1x', '0.2x', '1x', '5x', '10x', '50x', '100x','1000x')
title('\itG. vaginalis \rm\bfParameters')
ylabel('%Max Growth \itG. vaginalis')
set(gca, 'Fontsize', 14)

subplot(2,2,2)
gv_f = [gv_out(1:4,:); gv_out(6:8,:); gv_out(10:11,:)]*100;
li_f = [li_out(1:4,:); li_out(6:8,:); li_out(10:11,:)]*100;
pLabels = {'kgrow_L_I', 'kint_L_I', 'kkill_L_I', 'EC50_L_I'};
cmap1  = colormap(ax2, parula(11));
h1 = bar(gv_f(6:end,:), 'FaceColor', 'flat');
for c = 1:size(squeeze(output_data(2,:,:)),2)
    h1(c).CData = cmap1(c,:);
end
set(gca, 'XTick',1:length(pLabels),'XTickLabel', pLabels)
legend('0.001x', '0.01x', '0.02x', '0.1x', '0.2x', '1x', '5x', '10x', '50x', '100x','1000x')
title('\itL. iners \rm\bfParameters')
ylabel('%Max Growth \itG. vaginalis')
set(gca, 'Fontsize', 14)

subplot(2,2,3)
gv_f = [gv_out(1:4,:); gv_out(6:8,:); gv_out(10:11,:)]*100;
li_f = [li_out(1:4,:); li_out(6:8,:); li_out(10:11,:)]*100;
pLabels = {'kgrow_G_V', 'kkill_G_V', 'kint_G_V', 'kmet', 'EC50_G_V'};
cmap1  = colormap(ax2, parula(11));
h1 = bar(li_f(1:5,:), 'FaceColor', 'flat');
for c = 1:size(squeeze(output_data(2,:,:)),2)
    h1(c).CData = cmap1(c,:);
end
set(gca, 'XTick',1:length(pLabels),'XTickLabel', pLabels)
legend('0.001x', '0.01x', '0.02x', '0.1x', '0.2x', '1x', '5x', '10x', '50x', '100x','1000x')
title('\itG. vaginalis \rm\bfParameters')
ylabel('%Max Growth \itL. iners')
set(gca, 'Fontsize', 14)

subplot(2,2,4)
gv_f = [gv_out(1:4,:); gv_out(6:8,:); gv_out(10:11,:)]*100;
li_f = [li_out(1:4,:); li_out(6:8,:); li_out(10:11,:)]*100;
pLabels = {'kgrow_L_I', 'kint_L_I', 'kkill_L_I', 'EC50_L_I'};
cmap1  = colormap(ax2, parula(11));
h1 = bar(li_f(6:end,:), 'FaceColor', 'flat');
for c = 1:size(squeeze(output_data(2,:,:)),2)
    h1(c).CData = cmap1(c,:);
end
set(gca, 'XTick',1:length(pLabels),'XTickLabel', pLabels)
legend('0.001x', '0.01x', '0.02x', '0.1x', '0.2x', '1x', '5x', '10x', '50x', '100x','1000x')
title(['\it L. iners\rm \bfParameters'])
ylabel('%Max Growth \itL. iners')
set(gca, 'Fontsize', 14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3) Sensitivity ratio (Fig. 2c)
clear;
% Inputs for the function:
t = 48; % time (hr)
dGvLi = [1/1000 1/100 1/50 1/10 1/5 1 5 10 50 100 1000];
MNZ = 500;
load('params.mat')

[output_dat, gv_killed, li_killed] = Sense_1D_ratio(t,dGvLi,params,MNZ);

%%%%%%%%%%%%%%%%%%%%%% Plot the data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
ax2 = figure(3);
int_LB = [gv_killed; li_killed]*100;

cmap1  = colormap(ax2, parula(11));
h1 = bar(int_LB, 'FaceColor', 'flat');
for c = 1:size(squeeze(output_dat(2,:,:)),2)
    h1(c).CData = cmap1(c,:);
end
set(gca, 'XTick',1:3,'XTickLabel', {'G. vaginalis', 'L. iners'})
legend('0.001x', '0.01x', '0.02x', '0.1x', '0.2x', '1x', '5x', '10x', '50x', '100x','1000x')
ylabel('% Max Growth')
set(gca, 'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



