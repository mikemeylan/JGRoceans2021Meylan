function R = elastic_plate_R(L,floe_thickness,H,omega,E,rho_w,rho_i,nu,number_of_roots)
% program to calculate the far field scattering in Fourier cosine basis for an elastic ice floe
% L =   floe length
% floe_thickness =  floe thickness
% H = water depth
% omega =  freqency in radians
% theta =  angles to calculate D
% E =  Youngs modulus,
% rho_w  density of the water,
% rho_i  density of the ice
% nu = 0.3; % Poissons ratios
% number_of_roots = number of vertical eigenmodes
% N = number of fourier modes
%
% D(theta) = sum D_coefficents*cos(n*theta); 

g = 9.81;
% code uses these non-dimensional variables
alpha = omega^2/g;
D = E*floe_thickness^3/(12*(1-nu^2));
beta = D/(g*rho_w);
gamma = rho_i*floe_thickness/rho_w;

% [mat_a,~]  = circle_solve(alpha,beta,gamma,nu,H,a,N,number_of_roots);
% k=-dispersion_free_surface(alpha,0,H);


   [a_s,b_s,a_a,b_a] = finite_plate_symmetry(alpha,beta,gamma,H,number_of_roots,L,0,nu);
   R(1) = 1/2*(a_s(1) + a_a(1));
   R(2) = 1/2*(a_s(1) - a_a(1));
   
end

