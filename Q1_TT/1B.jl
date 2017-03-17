#This code generates a plot of the mRNA levels as a function of time for two cases:
#when RNAP is pre-incubated with the DNA and when RNAP is not pre-incubated with DNA.
using ODE
using PyPlot

#Constants from JV's datadictionary
cell_diameter = 1.1                 # mum
number_of_rnapII = 4600                # copies/cells
doubling_time_cell = 100           # hrs
max_transcription_rate = 60.0       # nt/sec
average_transcript_length = 1200       # nt
fraction_nucleus = 0.0                 # dimensionless
av_number = 6.02e23                 # number/mol
gene_coding_length=1080;          #characteristic length
w=0.1;

#Calculated constants
#Cell volume
V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

#Calculate the rnapII_concentration and ribosome_concentration
rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM

# Maximum specific growth rate -
maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1

# How fast do my cells die?
death_rate_constant = 0.2*maximum_specific_growth_rate                            # hr^-1

#Calculated terms for the rate expression
kt=max_transcription_rate*average_transcript_length;
rt=rnapII_concentration;
Lt=gene_coding_length;
Ltj=average_transcript_length;
Gj=5; #amount of genes
Ktj=10; #gene dosage
kdj=death_rate_constant;
rtj_bar=kt*rt*(Lt/Ltj)*(Gj/(Ktj+Gj)); #kinetic rate

#ODE for pre-incubated case
function ode(t, m)
  uj=w/(1+w)
  rtj=rtj_bar*uj
  dmjdt=rtj-(kdj)*m;
end;

#ODE for non-incubated case
function ode2(t2,m2)
  uj=(w/(1+w))*(1-(1/(1+exp(t2-3))))
  rtj=rtj_bar*uj
  dmjdt2=rtj-(kdj)*m2;
end;

#Initial conditions (amount of mRNA at t=0)
IC = 0.0;

#Time Vector
time = 0:0.1:10;
time2 = 0:0.1:10;

#ODE solver
t, m = ode45(ode, IC, time);
t2, m2 = ode45(ode2,IC,time2);
x = map(m -> m[1], m)
y = map(m2 -> m2[1], m2)


using PyPlot
plot(t, x, label="Pre-incubated")
plot(t2, y, label="Non-incubation")
xlabel("Time")
ylabel("mRNA Levels")
title("mRNA Levels vs Time")
legend();
