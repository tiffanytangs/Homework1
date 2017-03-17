#Finds what identifiable parameters(IP) and how many (lIP).
#Change the CHANGE here to change the file paths (3 total) depending on which
#sampling window you want to use: pre, during, post inducer
#and for which time skip (1 total)
include("Include.jl")

data_dictionary=DataDictionary(0.0,0.0,0.0)
#CHANGE PATH HERE dependng on what you want to run
# path_to_sensitivity_files="./sensitivity_ss"
path_to_sensitivity_files="./sensitivity_post"
# path_to_sensitivity_files="./sensitivity_ss"
file_pattern="AdjSimulation-P"

#CHANGE HERE FOR TIME SKIPS
# time_skip=1;
# time_skip=5;
#time_skip=10;
time_skip=20;

#Identifiable parameters from all species
(T,SA)=calculate_sensitivity_array(path_to_sensitivity_files,file_pattern,time_skip,data_dictionary)
IP=estimate_identifiable_parameters(SA,.1)
lIP=length(IP) #tells you how many identifiable
#CHANGE HERE
file_path = "./sensitivity_post/IP-"*string(time_skip)*".dat"
writedlm(file_path, IP)

#Identifiable parameters from mRNA1 and P3
#Generating sensitivity array for mRNA1
ts=length(T)
SA_mRNA1=zeros(ts,24)
for ii=1:ts
  SA_mRNA1[ii,:]=SA[[4+9*(ii-1)],:]
end

#Generating sensitivity array for P3
SA_P3=zeros(ts,24)
for jj=1:ts
  SA_P3[jj,:]=SA[[9+9*(jj-1)],:]
end

#Merging mRNA1 and P3 sensitivity arrays
SA_mRNA1_P3=[SA_mRNA1;SA_P3];
#Identifiable parameters for mRNA1 and P3 measurements
IP_mRNA1_P3=estimate_identifiable_parameters(SA_mRNA1_P3,.1)
lIP_mRNA1_P3=length(IP_mRNA1_P3) #tell you how many identifiable
#CHANGE HERE FILE PATH
file_path = "./sensitivity_post/IP_mRNA1_P3-"*string(time_skip)*".dat"
writedlm(file_path, IP_mRNA1_P3)
