function [mat_a,mat_b] = circle_solve(alpha,beta,gamma,nu,H,a,N,number_of_roots)
% [mat_a,mat_b] = circle_solve(alpha,beta,gamma,nu,H,a,N,number_of_roots)
% this program solves for the coefficents in the expansion of the circular plate.
% alpha, beta and gamma are the standard values.
% nu is poisson's ratio. H is the water depth. a is the
% radius of the plate. N is the number of angular modes and
% number_of_roots +1 is the number of roots of the dispersion equation 



% now we must solve for the coefficients of displacement for each 
% incoming mode
k=dispersion_free_surface(alpha,number_of_roots,H); k(1) = -k(1);
kappa=dispersion_elastic_surface(alpha,beta,gamma,number_of_roots+2,H);
kappa(3) = -kappa(3);

vec_e = 1.^(-N:N)/(1i*sqrt(alpha));% this is the vector of coefficents of the incoming waves
% vec_e = ones(size(vec_e));
%now we calculate vec_a and vec_b for each component
mat_a = zeros(number_of_roots+1,2*N+1);
mat_b = zeros(number_of_roots+3,2*N+1);
for n = -N:N
   [mat_a_temp,mat_b_temp] = circle_plate_matching_one_n(k,kappa,nu,H,a,n);
%    [mat_a_temp,mat_b_temp] = circle_plate_matching_one_n2(alpha,beta,gamma,nu,H,a,n,number_of_roots);
   
   
   mat_a(:,n+N+1) = vec_e(n + N + 1)*mat_a_temp;
   mat_b(:,n+N+1) = vec_e(n + N + 1)*mat_b_temp;
end


   