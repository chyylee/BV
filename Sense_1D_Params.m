function [gv_out, li_out, output_data] = Sense_1D_Params(t,dParams,MNZ,params,GvLi_ratio)
    % This function generates a 1D perturbatiton/sensitivity analysis where
    % the parameters of the model are varied based on the dParams input.
    
    % Input description:
    % t: time in hrs
    % dParams: N-length vector of scaling factors (ie 1000 = 1000*{original parameter
    % value})
    % MNZ: dose of MNZ
    % params: model parameter values in M-length cell array
    % GvLi_ratio: ratio of Gv to Li (1000 = 1000x Gv to Li)
    
    % Output description:
    % gv_out: MxN matrix of scaled gv growth to unperturbed simulation
    % li_out: MxN matrix of scaled li growth to unperturbed simulation
    % output_data: 6xMxN matrix of all ODE outputs 
    
    tot_cell = 2;
    foldx = GvLi_ratio;
    lb = tot_cell ./ (foldx + 1);
    gv = tot_cell - lb;
    cot = 0;
    
    tspan = [0,t];
    options = [];
    
    % Run unperturbed in silico co-culture
    y0 = [lb gv 0 0 0 0]; 
    [t0, ynull] = ode45(@GvLi_ODE_function, tspan, y0, options, params);
    null_Li = ynull(end,1);
    null_Gv = ynull(end,2);

    % Run 1D perturbation analysis
    y0 = [lb gv MNZ 0 0 0]; 
    for i = 1:length(dParams)
        for j = 1:length(params)
            newP = params;
            newP{j} = dParams(i)*params{j};
                [t, y] = ode45(@GvLi_ODE_function, tspan, y0, options, newP);
                for ode_num = 1:length(y(1,:))
                    output_data(ode_num,j,i) = y(end,ode_num);
                end
            cot = cot + 1;
            check_p(cot,:) = newP;
            gv_out(j,i) = y(end,2)./null_Gv;
            li_out(j,i) = y(end,1)./null_Li;
        end
    end

    save('sense_1D_params', 'gv_out', 'li_out', 'dParams', 'MNZ', 'GvLi_ratio', 'output_data','null_Li','null_Gv')
end