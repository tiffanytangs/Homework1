#Generates the plot of Number of Identifiable paramters vs Time Skips

using PyPlot

#Data obtained from running 'UnitTestSensitivity.jl'
#x corresponds to the time skips and y and z are the number of identifiable
#paramters for measurements of mRNA1 and P3, and all species, respectively
x=[1,5,10,20];
y=[6,5,5,4];
z=[13,10,9,7];

#Plots the above data
scatter(x,y,label="mRNA1 and P3")
scatter(x,z,label="All Species")
xlabel("Time Skip")
ylabel("Number of Identifiable Parameters")
title("During Inducer Number of Identifiable Parameters")
legend();
axis([0, 22, 0, 20]);
