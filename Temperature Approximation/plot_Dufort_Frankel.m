%!!! 100/Dx & 1000/Dt must be integers
%Dufort-Frankel scheme is unconditionally stable BUT Dt/Dx must approach 0 for convergence and precision
%Same as Crank-Nicolson i get the same result with larger steps
%Returns maximum temperature at t_max = 1000 : 31.63 Celsius same as Euler and Crank-Nicolson
%e.g. Dx = 5 & Dt = 1  or  Dx = 2 & Dt = 1

max_temp = Dufort_Frankel(Dx , Dt);
fprintf('Maximum Temperature : %f Celsius\n',max_temp);