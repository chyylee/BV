%% Compare final Gv and Li ratio trends
% May 11th, 2020
% Run file for 4 species model with interactions

%% Load PRCC result (parameter sets)
load('Model_LHS.mat')

%% Compare ratios

dParams = [0.6 100]; % 0.6x and 100 x
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

        y0 = [lb gv 500 0 0 0 lb gv 0 0 0]; % with ABX
        [t,y] = ode45(@wint_4sp_ode,tspan, y0, options, params);

        y0 = [lb gv 0 0 0 0 lb gv 0 0 0]; % "control"
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
%% Plot results
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

