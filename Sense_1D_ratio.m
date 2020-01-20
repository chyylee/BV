function [output_dat, gv_killed, li_killed] = Sense_1D_ratio(t,dGvLi,params,MNZ)
    % This function generates a 1D perturbatiton/sensitivity analysis where
    % the GvLi ratio is varied.
    
    % Input description:
    % t: time in hrs
    % GvLi: M-length vector of scaling factors (ie 1000 = 1000*{original input
    % value})
    % MNZ: dose of MNZ
    % params: model parameter values in N-length cell array
    
    % Outputs:
    % output_dat: 6xM matrix of all ODE outputs
    % gv_killed: scaled gv cell count to unperturbed simulation
    % li_killed: scaled li cell count to unperturbed simulation

    options = [];
    tspan = [0:0.1:t];
    dParams = dGvLi;
    for j = 1:length(dParams)
        tot_cell = 2;
        foldx = dParams(j);
        lb = tot_cell ./ (foldx + 1);
        gv = tot_cell - lb;
        y0 = [lb gv MNZ 0 0 0]; 
        [t,y] = ode45(@GvLi_ODE_function,tspan, y0, options, params);

        y0n = [lb gv 0 0 0 0];
        [tn,ynull] = ode45(@GvLi_ODE_function,tspan, y0n, options, params);

        output_dat(:,j) = y(end,:)./y(end,:);

        gv_killed(j) = y(end,2)./ynull(end,2); % Divide by carrying capacity
        li_killed(j) = y(end,1)./ynull(end,2); % Divide by carrying capacity
    end

    save('Sens_1D_ratio', 't', 'dGvLi', 'output_dat', 'gv_killed', 'li_killed', 'MNZ')
end
