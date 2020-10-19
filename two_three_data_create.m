%% program to create the data

clear all

radius = [5:5:200]; %  floe radius
floe_thickness = [0.25:0.25:2]; % floe thickness
period = [5:20]; % period



R_and_T = R_and_T_make(radius, floe_thickness,period);
Fourier_cosine_D = Fourier_cosine_D_make(radius, floe_thickness,period);


%%


save two_three_data

