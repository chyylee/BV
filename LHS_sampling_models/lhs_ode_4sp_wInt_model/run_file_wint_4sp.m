%% Compare final Gv and Li ratio trends
% May 11th, 2020
% Run file for 4 species model with interactions

%% Load PRCC result (parameter sets)
load('Model_LHS_fixed_Atassi.mat')

%%
tot_cell = 1;
p_BV = [0.01 0.1 0.5 0.9 0.99]*tot_cell;
p_LB = tot_cell - p_BV;
rat_BVLB = p_BV./p_LB;

dParams = [0.6 100];
yCompare = zeros(length(dParams),9,size(icMatrix,1));
% Call ODE function:

for j = 1:length(dParams)
    tot_cell = 2;
    foldx = dParams(j);
    lb = tot_cell ./ (foldx + 1);
    gv = tot_cell - lb;

    tspan = [0:0.1:48];
    options = [];

    yOut = zeros(9,size(paramMatrix,1));
    for i = 1:size(paramMatrix,1)
        params = paramMatrix(i,:);

        y0 = [lb gv 0 0 0 0 lb gv 0 0 0];
        [t,y] = ode45(@wint_4sp_ode,tspan, y0, options, params);

        y0 = [lb gv 0 0 0 0 lb gv 0 0 0];
        [tn,yn] = ode45(@wint_4sp_ode,tspan, y0, options, params);
        
        NY = y(end,:)./yn(end,:);
        
        
        yOut(:,i) = [NY(1);NY(2);NY(7);NY(8); NY(1)/(NY(1) + NY(2));
            (NY(1) + NY(7))/(NY(1) + NY(2) + NY(7) + NY(8)); 
            (y(end,1)+y(end,7))/(y(end,1)+y(end,2)+y(end,7)+y(end,8));
            (y(end,1))/(y(end,1)+y(end,2));
            (y(end,7))/(y(end,7)+y(end,8))];
        
    end
    yCompare(j,:,:) = yOut;
end
%% More streamlined results
graph_dat = yCompare(:,:,:);
ynames = {'ABX/NO ABX GROWTH'};

spnames = {'LI','GV','LJ','Other','Li/(Gv+Li)','Tot LB/(Total)','%LB','%LI','%Lj'};
xnames = {NaN(1,size(graph_dat,1))};

cmap = parula(size(graph_dat,2));
clear fp
for j = 1:size(graph_dat,2)
    subplot(2,ceil(size(graph_dat,2)/2),j)
    temp = squeeze(graph_dat(:,j,:));
    boxplot(temp')
    hold on
    for i = 1:size(graph_dat,1)
        scatter(i*ones(size(temp(i,:))).*(1+(rand(size(temp(i,:)))-0.5)/(i*5)),temp(i,:),'filled')
        xnames{i} = [num2str(dParams(i)),'x'];
        %[0.64,0.08,0.18]
        
    end
    [~,p] = ttest(temp(1,:),temp(2,:))
    fp(j,:) = [temp(1,:),temp(2,:)];
    xticklabels(xnames)
    ylabel(ynames)
    title(spnames(j))
    set(gca,'fontsize',14)
end
fpf = fp([1:4,7],:)*100;

%% Check parameter distributions:

% Compare parameters for Gv vs Li
subplot(2,3,1)
boxplot([paramMatrix(:,1) paramMatrix(:,7)]) % growth rates
xticklabels({'Gv','Li'})
ylabel('Growth Rate')

subplot(2,3,2)
boxplot([paramMatrix(:,3) paramMatrix(:,8)]) % internalization rates
xticklabels({'Gv','Li'})
ylabel('Internalization Rate')

subplot(2,3,3)
boxplot([paramMatrix(:,5) paramMatrix(:,9)]) % Carrying Capacity
xticklabels({'Gv','Li'})
ylabel('Carrying Capacity')

subplot(2,3,4)
boxplot([paramMatrix(:,6) paramMatrix(:,11)]) % EC50
xticklabels({'Gv','Li'})
ylabel('EC50')

subplot(2,3,5)
boxplot([paramMatrix(:,4)]) % Metabolism Rate
xticklabels({'Gv'})
ylabel('Metabolism Rate')

%% Calculate PRCC
paramnames = ["k_{grow-Gv}", "k_{kill-Gv}", "k_{int-Gv}",...
    "k_{met-Gv}", "K_{Gv}", "EC50_{Gv}", "k_{grow-Li}", ...
    "k_{int-Li}", "K_{Li}", "k_{kill-Li}",  "EC50_{Li}", ...
    "k_{grow-Gv2}","k_{kill-Gv2}", "k_{int-Gv2}",...
    "k_{met-Gv2}", "K_{Gv2}", "EC50_{Gv2}", "k_{grow-Lj}", ...
    "k_{int-Lj}", "K_{Lj}", "k_{kill-Li}", "EC50_{Lj}", ...
    "s_{Li->Lj}","s_{Lj->Li}","s_{Li->Gv}","s_{Gv->Li}", ...
    "s_{Lj->Gv}","s_{Gv->Lj}","s_{Li->BV2}","s_{BV2->Li}",...
    "s_{Lj->BV2}","s_{BV2->Lj}","s_{Gv->BV2}","s_{BV2->Gv}","Inital Ratio"];

var_id = 7;
% Y = squeeze(graph_dat(2,7,:));
Y = yblock(:,var_id);
k = size(paramMatrix,2);
% LHSmatrix = paramMatrix;
LHSmatrix = xblock;
prcc = zeros(k, 1);
prcc_significance = zeros(k, 1);
for i = 1:k
    % The LHS matrix with the column for parameter i deleted.
    z = LHSmatrix;
    z(:,i) = [];

    % LHSmatrix(:,i) is Nx1, z is Nx(k-1) and Y is N x timePointCount
    % rho and p are 1 x timePointCount.
    [rho, p] = partialcorr(LHSmatrix(:,i), Y, z, 'type', 'Spearman');
    
    prcc(i,:) = rho;
    prcc_significance(i, :) = p;
end

figure()
superbar(prcc,'P',prcc_significance)
xticklabels(paramnames)
ylabel('PRCC')
xticks(1:length(prcc))
xtickangle(270)
title(paramnames(var_id))
%title(['Time = ',num2str(odeSettings.tspan(end)),' NS = ',num2str(length(LHSmatrix(:,1)))])
set(gca,'Fontsize',16)

%% Prep for PLSR

clear yblock xblock
% Reformat YBLOCK
s1 = size(graph_dat,3);
for i = 1:size(graph_dat,1)
    sp = s1*i - (s1 - 1);
    ep = s1*i;
    yblock(sp:ep,:) = squeeze(graph_dat(i,:,:))';
    xblock(sp:ep,:) = [paramMatrix repmat(dParams(i),s1,1)];
end

%% Create PLSR model


var_id = 7;
fn ="";
classes = {1,2};
yname = spnames(var_id);
ttl = spnames(var_id);
X = xblock;
y = yblock(:,var_id);
[rsq,qsq,mse] = PLSR_TB(X,y,paramnames,fn,classes,yname,ttl,1);

%%
% %%
% function plot_PLSR(XS,XL,PCTVAR,pcn1,pcn2,y,names)
% subplot(1,2,1)
% scatter(XS(:,pcn1), XS(:,pcn2),40,y,'filled');
% hl = refline([0,0]);
% hl.Color = 'k';
% vl = xline(0);
% xlabel(['PC',num2str(pcn1),' - ', num2str(round(PCTVAR(pcn1),2)), '%'])
% ylabel(['PC', num2str(pcn2),' - ', num2str(round(PCTVAR(pcn2),2)), '%'])
% set(gca,'fontsize',14)
% title('Scores Plot')
% colorbar()
% 
% subplot(1,2,2)
% scatter(XL(:,pcn1),XL(:,pcn2))
% %text(XL(:,pcn1),XL(:,pcn2), names)
% xlabel(['PC',num2str(pcn1),' - ', num2str(round(PCTVAR(pcn1),2)), '%'])
% ylabel(['PC', num2str(pcn2),' - ', num2str(round(PCTVAR(pcn2),2)), '%'])
% hl = refline([0,0]);
% hl.Color = 'k';
% vl = xline(0);
% set(gca,'fontsize',14)
% title('Coefficient Plot')
% end

%%
Y = abs(fp(7,1:2000) - fp(7,2001:4000));

paramnames = ["k_{grow-Gv}", "k_{kill-Gv}", "k_{int-Gv}",...
    "k_{met-Gv}", "K_{Gv}", "EC50_{Gv}", "k_{grow-Li}", ...
    "k_{int-Li}", "K_{Li}", "k_{kill-Li}",  "EC50_{Li}", ...
    "k_{grow-Gv2}","k_{kill-Gv2}", "k_{int-Gv2}",...
    "k_{met-Gv2}", "K_{Gv2}", "EC50_{Gv2}", "k_{grow-Lj}", ...
    "k_{int-Lj}", "K_{Lj}", "k_{kill-Li}", "EC50_{Lj}", ...
    "s_{Li->Lj}","s_{Lj->Li}","s_{Li->Gv}","s_{Gv->Li}", ...
    "s_{Lj->Gv}","s_{Gv->Lj}","s_{Li->BV2}","s_{BV2->Li}",...
    "s_{Lj->BV2}","s_{BV2->Lj}","s_{Gv->BV2}","s_{BV2->Gv}"];

calc_mat = [paramMatrix];
% calc_mat = [LHSmatrix];
k = size(calc_mat,2);
prcc = zeros(k, 1);
prcc_significance = zeros(k, 1);
for i = 1:k
    % The LHS matrix with the column for parameter i deleted.
    z = calc_mat;
    z(:,i) = [];

    % LHSmatrix(:,i) is Nx1, z is Nx(k-1) and Y is N x timePointCount
    % rho and p are 1 x timePointCount.
    [rho, p] = partialcorr(calc_mat(:,i), Y', z, 'type', 'Spearman');
    
    prcc(i,:) = rho;
    prcc_significance(i, :) = p;
end

figure()
superbar(prcc,'P',prcc_significance*length(prcc))
xticklabels(paramnames)
ylabel('PRCC')
xticks(1:length(prcc))
xtickangle(270)
%title(spnames(var_id))
%title(['Time = ',num2str(odeSettings.tspan(end)),' NS = ',num2str(length(LHSmatrix(:,1)))])
set(gca,'Fontsize',16)