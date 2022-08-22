%!!! 100/Dx & 1000/Dt must be integers
%Dx = 0.5 & Dt = 0.1 for stability and Dx = 0.5 & Dt = 0.2 for instability
%Returns maximum temperature at t_max = 1000.
%31.63 Celsius at stability and infinite temperature (divergence) at instability

max_temp = Euler(Dx , Dt);
fprintf('Maximum Temperature : %f Celsius\n',max_temp);
