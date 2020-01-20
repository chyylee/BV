function [surface_out, vpGvLi, vpMNZ] = surface_data(t,params,GvLi_range,MNZ_range)
    % This function creates a 2D perturbation analysis to analyze how
    % changing the MNZ dose and GvLi ratio influences Li or Gv growth.
    
    % Inputs:
    % t: time (hrs)
    % params: parameter values in cell array
    % GvLi_range: M-length vector with GvLi values 
    % MNZ_range: N-length vector with MNZ doses (ug/ml)
    
    % Outputs:
    % surface_out: 6xNxM matrix (6 ODEs, NxM different simulations)
    % vpGvLi: M-length vector GvLi ratios
    % vpMNZ: N-length vector of MNZ doses

    vparams1 = GvLi_range;
    vparams2 = MNZ_range;
    tspan = [0:0.1:t];

    for i = 1:length(vparams1)
        for j = 1:length(vparams2)
            tot_cell = 2;
            foldx = vparams1(i);
            lb = tot_cell ./ (foldx + 1);
            gv = tot_cell - lb;

            y0 = [lb gv vparams2(j) 0 0 0]; 
            options =[];
            [t,y] = ode45(@GvLi_ODE_function,tspan, y0, options, params);

            output_dat(:,j,i) = y(end,:);
        end
    end

    vparams2_null = zeros(1,length(vparams2));
    for i = 1:length(vparams1)
        for j = 1:length(vparams2_null)
            tot_cell = 2;
            foldx = vparams1(i);
            lb = tot_cell ./ (foldx + 1);
            gv = tot_cell - lb;

            y0 = [lb gv vparams2_null(j) 0 0 0]; 
            options =[];
            [t,y] = ode45(@GvLi_ODE_function,tspan, y0, options, params);

            output_dat_null(:,j,i) = y(end,:);
        end
    end
    

    surface_out = output_dat ./ output_dat_null;
    vpGvLi = vparams1;
    vpMNZ = vparams2;
    save('Sim_Exp_Surf_data', 'surface_out', 'vpGvLi', 'vpMNZ')
end
