%!!! 100/Dx & 1000/Dt must be integers
%If lambda <= 0.5 result is unconditionally stable
%If lambda >= 0.5 -> Dt/(Dx)^2 = 1/(2*(2*lambda-1)) e.g. lambda = 1.5, Dx = 2, Dt = 1
%Larger step Dt -> fewer computations -> same result
%Returns maximum temperature at t_max = 1000 : 31.63 Celsius same as Euler
%e.g. lambda = 0.2 , Dx = 1 , Dt = 1

max_temp = Crank_Nicolson(lambda , Dx , Dt);
fprintf('Maximum Temperature : %f Celsius\n',max_temp);
