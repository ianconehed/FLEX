# FLEX
Code for Cone, Clopath and Shouval 2024 "Learning to Express Reward Prediction Error-like Dopaminergic Activity Requires Plastic Representations of Time"

The code for the main FLEX model is located in \main_model.

"Main_vectorized_simple_dopa.m" includes the main code for the simulations - a running the program will run a demo simulation and plot relevant figures afterwards. Comments are included in the file for further instruction.

"spiking_parameters_simple_dopa.m" includes the parameters for use in the main file.

The code for data analysis is located in \data_analysis.

Note that to access the data, it will have to be sourced from the appropriate repositories (see references in paper). For the data from Coddington and Dudman 2018, it is available from the authors upon request.

The code for the versions of the model with multiple cues and value encoding is located in \multiple_cues_and_value_model.

It is structured like the main_model folder, so running "Main_vectorized_simple_dopa_m_c.m" will run a demo simulation with distractor cues and plot relevant figures afterwards. Comments are included in the file for further instruction.

"spiking_parameters_simple_dopa.m" includes the parameters for this simulation, as before.

To test the value network, follow the instructions in the comments "in value_wrapper.m"

This code was written and run in MATLAB 2023b. Please direct any questions to ianconehed(at)gmail(dot)com. It was tested both on Mac and Windows platforms. It should take no more than 30 minutes to run the code.