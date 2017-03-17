# i dont know what i am doing
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
gene_coding_length=100;          #characteristic length
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
Ltj_1=10;
Ltj_2=100;
Ltj_3=1000;
Gj=5; #amount of genes
Ktj=10; #gene dosage
kdj=death_rate_constant;

#ODE for L_tj=10
function ode(t, m)
  rtj_bar=kt*rt*(Lt/Ltj_1)*(Gj/(Ktj+Gj))
  uj=w/(1+w)
  rtj=rtj_bar*uj
  dmjdt=rtj-(kdj)*m;
end;

#ODE for L_tj=100
function ode2(t2, m2)
  rtj_bar=kt*rt*(Lt/Ltj_2)*(Gj/(Ktj+Gj))
  uj=w/(1+w)
  rtj=rtj_bar*uj
  dmjdt2=rtj-(kdj)*m2;
end;

#ODE for L_tj=10
function ode3(t3, m3)
  rtj_bar=kt*rt*(Lt/Ltj_3)*(Gj/(Ktj+Gj))
  uj=w/(1+w)
  rtj=rtj_bar*uj
  dmjdt3=rtj-(kdj)*m3;
end;

#Initial conditions (amount of mRNA at t=0)
IC = 0.0;

#Time Vector
time = 0:0.1:10;
time2 = 0:0.1:10;
time3 = 0:0.1:10;

#ODE solver
t, m = ode45(ode, IC, time);
t2, m2 = ode45(ode2,IC,time2);
t3, m3 = ode45(ode3,IC,time3);
x = map(m -> m[1], m)
y = map(m2 -> m2[1], m2)
z = map(m3 -> m3[1], m3)

using PyPlot
plot(t, x, label="L_tj=10 bp")
plot(t2, y, label="L_tj=100 bp")
plot(t3, z, label="L_tj=1000 bp")
xlabel("Time")
ylabel("mRNA Levels")
title("Q1c. mRNA Levels vs Time")
legend();
