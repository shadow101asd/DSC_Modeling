function f = DSC_TVG_Plotting_fromDSC4(figure_num, filename, idx, planet)
%DSC_TVG_Plotting Summary of this function goes here
%   Detailed explanation goes here
    load(filename)
    
    es_opt = linspace(emin_optN(idx), emax_optN(idx), numshells_optN(idx));
    shells_opt = ones(1,numshells_optN(idx))*numbranches_optN(idx);

    if planet == "Mars"
        f = DSC_TVG_Plotting(figure_num, XEa, XMa, "Earth", "Mars", "blue", "red", ...
            muSu, etR, Ki, es_opt, shells_opt, numcarts_optN(idx));
    elseif planet == "Venus"
        f = DSC_TVG_Plotting(figure_num, XEa, XVe, "Earth", "Venus", "blue", "green", ...
            muSu, etR, Ki, es_opt, shells_opt, numcarts_optN(idx));
    else
        error("Wrong planet name!!")
    end
end

