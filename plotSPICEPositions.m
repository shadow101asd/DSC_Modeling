function figure = plotSPICEPositions(X, fontstring, colorstring)
% Helper function for scatterplotting ephemeride data
    figure = scatter3(X(1,:),X(2,:),X(3,:), fontstring, colorstring);
end

