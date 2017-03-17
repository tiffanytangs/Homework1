#Nonlinear model analysis as discussed in Kremling et al.
#Change these if time or state arrays change size
#RUN WASHOUT FILE for each model and fetch the amount of times and number of states
time_array=collect(1:1:251);
state_array=collect(1:1:12);

#Setting up variables
len_time_array=length(time_array)
len_state_array=length(state_array)
Trapezoid_all=zeros(len_state_array,1)

#Calculate objective function for each state
for state_index=state_array
#Read the state files for model A and model B - CHANGE TO FILE PATHS OF THE STATE DAT FILES
file_path_A = "C:/Users/Tiffany/Desktop/Julia/Processed_Files/Q3a/src/States/Washout-S"*string(state_index)*".dat"
xa=readdlm(file_path_A)
file_path_B = "C:/Users/Tiffany/Desktop/Julia/Processed_Files/Q3b/src/States/Washout-S"*string(state_index)*".dat"
xb=readdlm(file_path_B)

#calculate objective function at one instance of time
y=zeros(len_time_array,1)
for i=1:len_time_array
y[i]=((xb[i]-xa[i])^2)/((1/(xb[i]-xa[i]))^2)
end

#Perform trapezoid rule
z=zeros(len_time_array-1,1)
for j=1:len_time_array-1
  z[j]=((y[j]+y[j+1])/2)*(time_array[j+1]-time_array[j])
end
Trapezoid_Rule=sum(z)

#Stores result of trapezoid rule
Trapezoid_all[state_index]=Trapezoid_Rule
end
