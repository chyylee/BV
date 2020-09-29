%% Compare final Gv and Li ratio trends
% May 11th, 2020
% Script to evaluate model results based on LHS of parameter ranges for the
% 4 species model with no interspecies interactions

%% Load PRCC result (parameter sets)

load('Model_LHS.mat')

%%
tot_cell = 1;
p_BV = [0.4 0.7 0.95]*tot_cell;
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

        [t,y] = ode45(@Noint_4sp_ode,tspan, y0, options, params);
        
        y0 = [lb gv 0 0 0 0 lb gv 0 0 0];
        [tn,yn] = ode45(@Noint_4sp_ode,tspan, y0, options, params);
        
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
%% See results

gv1000 = squeeze(yCompare(1,2,:));
gv0001 = squeeze(yCompare(2,2,:));

Oth1000 = squeeze(yCompare(1,4,:));
Oth0001 = squeeze(yCompare(2,4,:));

li1000 = squeeze(yCompare(1,1,:));
li0001 = squeeze(yCompare(2,1,:));

lj1000 = squeeze(yCompare(1,3,:));
lj0001 = squeeze(yCompare(2,3,:));

subplot(1,2,1)
boxplot([gv1000 gv0001 Oth1000 Oth0001])
hold on;
scatter(ones(size(gv1000)).*(1+(rand(size(gv1000))-0.5)/10),gv1000,'r','filled')
scatter(2*ones(size(gv0001)).*(1+(rand(size(gv0001))-0.5)/10),gv0001,'r','filled')
scatter(3*ones(size(Oth1000)).*(1+(rand(size(Oth1000))-0.5)/10),Oth1000,'b','filled')
scatter(4*ones(size(Oth0001)).*(1+(rand(size(Oth0001))-0.5)/10),Oth0001,'b','filled')

ylabel('%Max growth Non-Optimal')
xticklabels({'Gv 1000x','Gv 0.001x','Other 1000x','Other 0.001x'})
set(gca,'fontsize',14)

subplot(1,2,2)
boxplot([li1000 li0001 lj1000 lj0001])
hold on
scatter(ones(size(li1000)).*(1+(rand(size(li1000))-0.5)/10),li1000,'r','filled')
scatter(2*ones(size(li0001)).*(1+(rand(size(li0001))-0.5)/10),li0001,'r','filled')
scatter(3*ones(size(lj1000)).*(1+(rand(size(lj1000))-0.5)/10),lj1000,'b','filled')
scatter(4*ones(size(lj0001)).*(1+(rand(size(lj0001))-0.5)/10),lj0001,'b','filled')
%ranksum(li1000,li0001)
ylabel('%Max growth Optimal')
set(gca,'fontsize',14)
xticklabels({'Li 1000x','Li 0.001x','Lj 1000x','Lj 0.001x'})

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

%% PRCC with LB difference
% new Y block
Y = abs(- fp(7,1:2000) + fp(7,2001:4000));

paramnames = ["k_{grow-Gv}", "k_{kill-Gv}", "k_{int-Gv}",...
    "k_{met-Gv}", "K_{Gv}", "EC50_{Gv}", "k_{grow-Li}", ...
    "k_{int-Li}", "K_{Li}", "k_{kill-Li}",  "EC50_{Li}", ...
    "k_{grow-Gv2}","k_{kill-Gv2}", "k_{int-Gv2}",...
    "k_{met-Gv2}", "K_{Gv2}", "EC50_{Gv2}", "k_{grow-Lj}", ...
    "k_{int-Lj}", "K_{Lj}", "k_{kill-Li}", "EC50_{Lj}"];

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