function [vec_a,vec_b]= circle_plate_matching_one_n(mroots_water,mroots_ice,nu,H,a,n)
% this program solves the circlular plate problem for one
% value of n (the coefficient in the expansion of e^(i n\theta)
%
% mroots_water and mroots_ice are the roots of the dispersion 
% equations with positive imaginary part (for the travelling)
% and positive part for the damped roots
% 
% mroots_ice = dispersion_ice_standard(alpha*H,beta/H^4,gamma/H,number_of_roots+2)/H;
%
% and the dispersion equation for the water
% mroots_water = dispersion_standard(alpha*H,number_of_roots)/H;
%
% nu is Poisson's ratio and H is the water depth. 
% a is the radius of the plate. 
%
% vec_a is the vector of coefficients in the expansion of the water
% vec_b is the vector of coefficients in the expansion of the plate


number_of_roots = length(mroots_water) - 1;

%
% This is the old method which produces the same result as the more efficent method we use
%
% % the  matrics of innerproducts
% M1 = zeros(number_of_roots+1,number_of_roots+1);
% for j = 1:number_of_roots+1;
%       M1(j,j) = besselk(n,mroots_water(j)*a)*innerproduct_semiinfinite(mroots_water(j),mroots_water(j),H);
% end
% 
% M2 = zeros(number_of_roots+1,number_of_roots+3);
% for j = 1:number_of_roots+1;
%    for k = 1:number_of_roots+3
%       M2(j,k) = besseli(n,mroots_ice(k)*a)*innerproduct_semiinfinite(mroots_water(j),mroots_ice(k),H);
%   end
% end
% 
% M3 = zeros(number_of_roots+1,number_of_roots+1);
% for j = 1:number_of_roots+1;
%       M3(j,j) = mroots_water(j)*(dbesselk(n,mroots_water(j)*a))*innerproduct_semiinfinite(mroots_water(j),mroots_water(j),H);
% end
% 
% M4 = zeros(number_of_roots+1,number_of_roots+3);
% for j = 1:number_of_roots+1;
%    for k = 1:number_of_roots+3
%       M4(j,k) = mroots_ice(k)*dbesseli(n,mroots_ice(k)*a)*innerproduct_semiinfinite(mroots_water(j),mroots_ice(k),H);
%   end
% end
% 
% % now we have to define the matrix A
A = zeros(2,number_of_roots+3);
for j = 1:number_of_roots+3;
    A(1,j) = mroots_ice(j)*tan(mroots_ice(j)*H)*(mroots_ice(j)^2*besseli(n,mroots_ice(j)*a) - ...
        (1-nu)/a*(mroots_ice(j)*dbesseli(n,mroots_ice(j)*a) - n^2/a*besseli(n,mroots_ice(j)*a)));
    A(2,j) = (mroots_ice(j))*tan(mroots_ice(j)*H)*(mroots_ice(j)^3*dbesseli(n,mroots_ice(j)*a) + ...
        (1-nu)/a^2*(-n^2*mroots_ice(j)*dbesseli(n,mroots_ice(j)*a) + n^2/a*besseli(n,mroots_ice(j)*a)));
end

% finally we need to define the vector of unknowns which is 
% f = zeros(number_of_roots+1 + number_of_roots+3 , 1);
% f(1) = -besseli(n,mroots_water(1)*a)*innerproduct_semiinfinite(mroots_water(1),mroots_water(1),H);
% f(number_of_roots + 1 + 1) = -mroots_water(1)*dbesseli(n,mroots_water(1)*a)*innerproduct_semiinfinite(mroots_water(1),mroots_water(1),H);
% 
% Now we solve the system of equations
% bigM = [M1,-M2;M3,-M4;zeros(2,number_of_roots+1),A];
% solution = bigM\f;
% bigM
% f
% vec_a = solution(1:number_of_roots+1);
% vec_b = solution(number_of_roots+2:length(solution));

% we can solve the system of equations more effiently by eliminating
% vec_a. This system is also more stable. 


B = zeros(number_of_roots+1,number_of_roots+3);
for p = 1:number_of_roots+1;
   for m = 1:number_of_roots+3
      B(p,m) = (mroots_ice(m)*dbesseli(n,mroots_ice(m)*a) - ...
         mroots_water(p)*dbesselk(n,mroots_water(p)*a)*besseli(n,mroots_ice(m)*a)/besselk(n,mroots_water(p)*a)) ... 
         *innerproduct_semiinfinite(mroots_water(p),mroots_ice(m),H);
  end
end
% B;
f = zeros(number_of_roots+3,1);

f(1) = (mroots_water(1)*dbesseli(n,mroots_water(1)*a) - ...
         mroots_water(1)*dbesselk(n,mroots_water(1)*a)*besseli(n,mroots_water(1)*a)/besselk(n,mroots_water(1)*a)) ... 
         *innerproduct_semiinfinite(mroots_water(1),mroots_water(1),H);
      
      vec_b = [B;A]\f;
      
      % now we will calculate vec_a
      
      C = zeros(number_of_roots+1,number_of_roots+3);
for p = 1:number_of_roots+1;
   for m = 1:number_of_roots+3
      C(p,m) = besseli(n,mroots_ice(m)*a) ... 
         *innerproduct_semiinfinite(mroots_water(p),mroots_ice(m),H) / ...
         (besselk(n,mroots_water(p)*a) ... 
         *innerproduct_semiinfinite(mroots_water(p),mroots_water(p),H));
   end
   
   vec_a = C*vec_b;vec_a(1) = vec_a(1) - besseli(n,mroots_water(1)*a)/besselk(n,mroots_water(1)*a) ;
end

  
function out = dbesseli(n,z)
% this function calcualtes the derivative of besselk(n,z)

out = (besseli(n+1,z) + besseli(n-1,z))/2;

function out = dbesselk(n,z)
% this function calcualtes the derivative of besselk(n,z)

out = -(besselk(n+1,z) + besselk(n-1,z))/2;

function out = innerproduct_semiinfinite(k,kappa,H)
% out = innerproduct_semiinfinite(k,kappa,H)
% calculates the innerproduct \int_{-H}^0 cos(k(z+H))/coskH cos(\kappa (z+H)/ cos \kappa H dz.
% there are two case depending on whether k= \kappa or not

if k == kappa
    if abs(k) < 1e-10
        out = H;
        
    else
		if real(k) < 1e-10
            out = (1/2)*((-1i*tanh(1i*k*H) + k*H*sech(-1i*k*H)^2)/k);

        else
            out = (1/2)*((tan(k*H) + k*H*sec(k*H)^2)/k);
		end
    end
    
else
    if real(k) < 1e-10
        tankH = -1i*tanh(1i*k*H);
    else
        tankH = tan(k*H);
    end
    
    if real(kappa) < 1e-10
        tankappaH = -1i*tanh(1i*kappa*H);
    else
        tankappaH = tan(kappa*H);
    end

    
     out =((tankH)*k -(tankappaH)*kappa)/(k^2-kappa^2);
%      out =((tan(k*H))*k -(tan(kappa*H))*kappa)/(k^2-kappa^2)
end



