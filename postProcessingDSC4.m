%% Loading Data

clear all
load DSC4_min1max200spacing1.mat

%% Post-Processing of Optimal Configurations

Best_DequivMean = 10^40;
for i = 1:N
    if Best_DequivMean > DequivMean(i)
        Best_DequivMean = DequivMean(i);
    else
        DequivMean(i) = DequivMean(i-1);
        numsats_optN(i) = numsats_optN(i-1);
        emin_optN(i) = emin_optN(i-1);
        emax_optN(i) = emax_optN(i-1);
        numshells_optN(i) = numshells_optN(i-1);
        numbranches_optN(i) = numbranches_optN(i-1);
        numcarts_optN(i) = numcarts_optN(i-1);
    end
end

% Now that DequivMean is monotone decreasing, we backpropagate

current_max_sats = MAX_MAXSATS+10;
for i = N:-1:1
    if i >= current_max_sats
        DequivMean(i) = DequivMean(i+1);
        numsats_optN(i) = numsats_optN(i+1);
        emin_optN(i) = emin_optN(i+1);
        emax_optN(i) = emax_optN(i+1);
        numshells_optN(i) = numshells_optN(i+1);
        numbranches_optN(i) = numbranches_optN(i+1);
        numcarts_optN(i) = numcarts_optN(i+1);
    else
        current_max_sats = numsats_optN(i);
    end
end

%% Recompute Perfomances

DequivMean(N) = 0.0;
DequivMin(N) = 0.0;
DequivMax(N) = 0.0;
DrealMean(N) = 0.0;
DrealMin(N) = 0.0;
DrealMax(N) = 0.0;

for i = 1:N
    % Recalculate obtained orbits
    es_opt = linspace(emin_optN(i),emax_optN(i),numshells_optN(i));
    shells_opt = numbranches_optN(i)*ones(numshells_optN(i),1);
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,shells_opt,es_opt,numcarts_optN(i));
    
    % Extract Performance Parameters
    [graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats);
    DequivMean(i) = mean(Dequivs)/AU;
    DequivMin(i) = min(Dequivs)/AU;
    DequivMax(i) = max(Dequivs)/AU;
    DrealMean(i) = mean(Dpaths)/AU;
    DrealMin(i) = min(Dpaths)/AU;
    DrealMax(i) = max(Dpaths)/AU;
end

DEaMa = distanceBetweenXs(XEa,XMa);

%% Saving data

filename = "DSC4_PP/PP_DSC4_min"+ num2str(MIN_MAXSATS) + "max" + num2str(MAX_MAXSATS) + "spacing" + num2str(SPACING);
save(filename)

%% Plotting (Optional)

c = 3*10^5; % Speed of light, km/s
fit_type = "power2";
[fitequiv, gof] = fit(SATNUMS',DequivMean', fit_type);
Rsquared = gof.rsquare;

figure(1)
plot(SATNUMS,DequivMean, 'Color', 'blue');
hold on
plot(fitequiv);
plot(SATNUMS,DequivMax, 'Color', 'blue', 'LineStyle', '--');
plot(SATNUMS,DequivMin, 'Color', 'blue', 'LineStyle', '--');
hold off
legend('Mean', sprintf('%s Function Fit: R^2 = %f', fit_type, Rsquared), 'Bounds at each # of satellites')
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Distance [AU]')
title('Equivalent Transmission Distance')

figure(2)
plot(SATNUMS,DrealMean*AU/c/60, 'Color', 'red');
hold on
plot(SATNUMS,DrealMax*AU/c/60, 'Color', 'red', 'LineStyle', '--');
yline(mean(DEaMa)/c/60, 'green');
yline(max(DEaMa)/c/60, 'Color', 'green', 'LineStyle', '--');
yline(min(DEaMa)/c/60, 'Color', 'green', 'LineStyle', '--');
plot(SATNUMS,DrealMin*AU/c/60, 'Color', 'red', 'LineStyle', '--');
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Transmission Time [min]')
title('Real Transmission Time')
legend('Constellation Mean', 'Bounds at each # of satellites', ...
    'Earth-Mars direct transmission time (mean)', 'Earth-Mars direct transmission time (max/min throughout orbits)')

figure(3)
plot(SATNUMS, numshells_optN);
hold on
plot(SATNUMS, numbranches_optN);
plot(SATNUMS, numcarts_optN);
plot(SATNUMS, numsats_optN);
fplot(@(x) x, [0 max(SATNUMS)]);
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('#')
title('Cartwheel Constellation Configuration Numbers as a Function of Number of Satellites')
legend('Number of Shells per CFFC', 'Number of Branches per CFFC', 'Number of CFFCs', 'Number of satellites in optimal constellation', 'y=x')

figure(4)
plot(SATNUMS, emin_optN);
hold on
plot(SATNUMS, emax_optN);
hold off
xlabel('Maximum permissible # of satellites in constellation')
ylabel('Eccentricity')
title('Cartwheel Constellation Eccentricities as a Function of Number of Satellites')
legend('emin', 'emax')