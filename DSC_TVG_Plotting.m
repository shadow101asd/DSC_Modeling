function f = DSC_TVG_Plotting(figure_num, X1, X2, name1, name2, color1, color2, muSu, etR, Ki, es_opt, shells_opt, cartwheels_opt)
%DSC_TVG_Plotting Summary of this function goes here
%   Detailed explanation goes here
    AU = 1.496e8;

    % Recompute XSats
    numsats_opt = sum(shells_opt)*cartwheels_opt;
    emin_opt = min(es_opt);
    emax_opt = max(es_opt);
    numshells_opt = length(shells_opt);
    numbranches_opt = max(shells_opt);
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells_opt,es_opt,cartwheels_opt);
    
    D12 = distanceBetweenXs(X1,X2);
    D1Sats = distanceBetweenXs(X1,XSats);
    D2Sats = distanceBetweenXs(X2,XSats);
    Xs = collateXs(X1,X2,XSats);
    [graphslist,~,paths,~,Dequivs] = networkAnalysis(X1,X2,name1,name2,XSats);
    [~,nT] = size(X1);
    
    % Plotting

    f = figure(figure_num);
    for t=1:nT
        clf(figure_num)
        viscircles([0,0], 2.2, Color="black");
        set(gca, "XLim", [-1.75, 1.75], "YLim", [-1.75, 1.75], "DataAspectRatio", [1, 1, 1]);
        hold on
        % Plot graphs
        p1 = plot(graphslist{t}, 'XData', reshape(Xs(1,t,:)/AU,numsats_opt+2,1,1), ...
            'YData', reshape(Xs(2,t,:)/AU,numsats_opt+2,1,1), "EdgeColor", [0.3010 0.7450 0.9330]);
        % Highlight Sun + Planets
        scatter(0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black');
        scatter(X1(1,t)/AU,X1(2,t)/AU,2e2, 'o', color1,'filled','MarkerEdgeColor', 'black');
        scatter(X2(1,t)/AU,X2(2,t)/AU,2e2, 'o', color2,'filled','MarkerEdgeColor', 'black');
        % Highlight paths
        highlight(p1, paths{t}, 'LineWidth', 6, 'EdgeColor', color2);
        hold off
        % Add axis labels and title
        xlabel('x [AU]');
        ylabel('y [AU]');
        title(sprintf(['Graph Plots Over Time, %s-%s MBSP. \n ' ...
            '%d Satellites: %d cartwheel formations, each with %d branches and %d shells. \n ' ...
            'Eccentricity range: [%f - %f]. Resulting ED: %f AU'], ...
            name2, name1, numsats_opt, cartwheels_opt, numbranches_opt, numshells_opt, emin_opt, emax_opt,mean(Dequivs)/AU));
        pause(0.1);
    end

end

