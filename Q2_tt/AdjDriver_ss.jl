# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
include("Include.jl")

# Setup the timescale of the simulation -
time_start = 0.0
time_stop = 240.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)

# Load the data dictionary (default parameter values) -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# What is the size of the system?
number_of_states = data_dictionary["number_of_states"]

# main loop -
initial_condition_array = data_dictionary["initial_condition_array"]
parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
average_scaled_sensitivity_array = zeros(number_of_states,1)
for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)

  # grab the dictionary -
  local_data_dictionary = deepcopy(data_dictionary)
  @show parameter_index

  # solve the Adj balances
  (T,X) = adj_washout_ss(time_start,time_stop,time_step_size,parameter_index,local_data_dictionary)
  data_array = [T X]
  file_path = "./sensitivity_ss/AdjSimulation-P"*string(parameter_index)*".dat"
  writedlm(file_path,data_array)

  #Calculating average scaled sensitivity array
  state_array = X[:,1:number_of_states]
  coolbeans = number_of_states+1
  sensitivity_array = X[:,coolbeans:end]
  scaled_sensitivity_array = scale_sensitivity_array(T,state_array,sensitivity_array,parameter_index,local_data_dictionary)

  # time average -
  average_sensitivity_col = time_average_array(T,scaled_sensitivity_array)

  # grab -
  average_scaled_sensitivity_array = [average_scaled_sensitivity_array average_sensitivity_col]
end

average_scaled_sensitivity_array = average_scaled_sensitivity_array[:,2:end]

#Producing heat map of the time average sensitivity array
abs_avg=abs(average_scaled_sensitivity_array)
pcolormesh(abs_avg, cmap="Reds")
colorbar()
xlabel("Parameters")
ylabel("States")
title("Pre-Inducer Time Averaged Sensitivity Array")

#Saving the time averaged sensitivity array to disk
file_path = "./sensitivity_ss/Avg_Scaled_Sensitivity_Array.dat"
writedlm(file_path, average_scaled_sensitivity_array)
