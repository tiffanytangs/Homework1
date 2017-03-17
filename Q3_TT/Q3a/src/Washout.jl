include("Include.jl")

# define some colors -
const shaded_color_value = "lightgray"
const mean_color_value = "dimgray"
const experimental_color_value = "black"

const P1_color = "blue"
const P1_shaded_color="skyblue"

const P2_color = "green"
const P2_shaded_color="lightgreen"

const P3_color = "orange"
const P3_shaded_color="navajowhite"

# Script to solve the balance equations -
time_start = 0.0
time_stop = 25.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)

# How many samples do we want to explore?
number_of_samples = 10
sigma = 0.20

# initialize storage -
time_array = []
P1_array = zeros(number_of_timesteps,number_of_samples)
P2_array = zeros(number_of_timesteps,number_of_samples)
P3_array = zeros(number_of_timesteps,number_of_samples)
X1=zeros(number_of_timesteps,number_of_samples)
X2=zeros(number_of_timesteps,number_of_samples)
X3=zeros(number_of_timesteps,number_of_samples)
X4=zeros(number_of_timesteps,number_of_samples)
X5=zeros(number_of_timesteps,number_of_samples)
X6=zeros(number_of_timesteps,number_of_samples)
X7=zeros(number_of_timesteps,number_of_samples)
X8=zeros(number_of_timesteps,number_of_samples)
X9=zeros(number_of_timesteps,number_of_samples)
X10=zeros(number_of_timesteps,number_of_samples)
X11=zeros(number_of_timesteps,number_of_samples)
X12=zeros(number_of_timesteps,number_of_samples)
# main loop -
for sample_index = 1:number_of_samples

  # Load the data dictionary (default parameter values) -
  data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

  # Pertub the RNAPs and Ribosome abundance -
  rnapII_concentration  = data_dictionary["rnapII_concentration"] # muM
  ribosome_concentration  = data_dictionary["ribosome_concentration"] # muM

  rnapII_concentration = abs(rnapII_concentration*(1+sigma*randn()))
  ribosome_concentration  = abs(ribosome_concentration*(1+sigma*randn()))

  data_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
  data_dictionary["ribosome_concentration"] = ribosome_concentration # muM

  # Run the washout simulation -
  (T,X) = washout_simulation(time_start,time_stop,time_step_size,data_dictionary)

  time_array = T
  # len_time_array=length(time_array)
  for (time_index,time_value) in enumerate(T)
    P1_array[time_index,sample_index] = X[time_index,9];
    P2_array[time_index,sample_index] = X[time_index,10];
    P3_array[time_index,sample_index] = X[time_index,11];
    X1[time_index,sample_index] = X[time_index,1];
    X2[time_index,sample_index] = X[time_index,2];
    X3[time_index,sample_index] = X[time_index,3];
    X4[time_index,sample_index] = X[time_index,4];
    X5[time_index,sample_index] = X[time_index,5];
    X6[time_index,sample_index] = X[time_index,6];
    X7[time_index,sample_index] = X[time_index,7];
    X8[time_index,sample_index] = X[time_index,8];
    X9[time_index,sample_index] = X[time_index,9];
    X10[time_index,sample_index] = X[time_index,10];
    X11[time_index,sample_index] = X[time_index,11];
    X12[time_index,sample_index] = X[time_index,12];
  end
end
X1_mean=mean(X1,2)
X2_mean=mean(X2,2)
X3_mean=mean(X3,2)
X4_mean=mean(X4,2)
X5_mean=mean(X5,2)
X6_mean=mean(X6,2)
X7_mean=mean(X7,2)
X8_mean=mean(X8,2)
X9_mean=mean(X9,2)
X10_mean=mean(X10,2)
X11_mean=mean(X11,2)
X12_mean=mean(X12,2)
len_time_array=length(X12_mean)
X_all=[X1_mean X2_mean X3_mean X4_mean X5_mean X6_mean X7_mean X8_mean X9_mean X10_mean X11_mean X12_mean]
X_s=zeros(len_time_array,1)

#Saves each state in a seperate .dat file
for state_index=1:12
  X_s=X_all[:,state_index]
  file_path = "./States/Washout-S"*string(state_index)*".dat"
  writedlm(file_path, X_s)
end

# Confidence interval?
SF = (1.96/sqrt(number_of_samples))

# Make some plots -
# P1 -
P1_mean = mean(P1_array,2)
P1_std = std(P1_array,2)
P1_lower_bound = P1_mean - SF*P1_std
P1_upper_bound = P1_mean + SF*P1_std
plot(time_array,P1_mean,lw=2,color=P1_color, label="P_1")
fill_between(time_array,vec(P1_lower_bound),vec(P1_upper_bound),color=P1_shaded_color,lw=3)

# P2 -
P2_mean = mean(P2_array,2)
P2_std = std(P2_array,2)
P2_lower_bound = P2_mean - SF*P2_std
P2_upper_bound = P2_mean + SF*P2_std
plot(time_array,P2_mean,lw=2,color=P2_color, label="P_2")
fill_between(time_array,vec(P2_lower_bound),vec(P2_upper_bound),color=P2_shaded_color,lw=3)

# P3 -
P3_mean = mean(P3_array,2)
P3_std = std(P3_array,2)
P3_lower_bound = P3_mean - SF*P3_std
P3_upper_bound = P3_mean + SF*P3_std
plot(time_array,P3_mean,lw=2,color=P3_color, label="P_3")
fill_between(time_array,vec(P3_lower_bound),vec(P3_upper_bound),color=P3_shaded_color,lw=3)


xlabel("Time")
ylabel("Protein Concentration")
title("Model A Concentration Profile")
legend();
# Dump to disk -
savefig("../figs/Washout-3G.pdf")
