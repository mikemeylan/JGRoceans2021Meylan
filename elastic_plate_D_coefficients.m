function D_coefficents = elastic_plate_D_coefficients(a,floe_thickness,H,omega,E,rho_w,rho_i,nu,number_of_roots,N)
% program to calculate the far field scattering in Fourier cosine basis for an elastic ice floe
% a =   floe radius
% floe_thickness =  floe thickness
% H = 1 water depth
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

[mat_a,~]  = circle_solve(alpha,beta,gamma,nu,H,a,N,number_of_roots);
k=-dispersion_free_surface(alpha,0,H);


    D_coefficents = 2*sqrt(alpha)*sqrt(pi/(2*1i*k))*mat_a(1,N+1:end);
    D_coefficents(1) = D_coefficents(1)/2;
end

