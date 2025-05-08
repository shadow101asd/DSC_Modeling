a = AU;
e = 0.0;
i = 1.0;
Om = 5.0;
w = 4.0;
f0i = 0.0;

f0s = [];
Xs = [];
numdays = 365*2;
dts = 0:3600:(3600*24*numdays);

for dt = 0:3600:(3600*24*numdays)
    f0new = updateTrueAnomaly(a,e,i,Om,w,f0i,muSu,dt);
    f0s = [f0s, f0new];
    Xs = [Xs,Keplerian2Cartesian(a,e,i,Om,w,f0new,muSu)];
end


figure(3)
plot(dts,f0s)

figure(4)
plotSPICEPositions(Xs, '.', 'red');