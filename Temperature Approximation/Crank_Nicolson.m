function max_temp = Crank_Nicolson(lambda , Dx , Dt)

% Parameters
L = 100;
t_max = 1000;
lambda = 0.2;
Dx = 1;
Dt = 1;
r = Dt/Dx^2;
N = L/Dx + 1;
M = t_max/Dt + 1;
To = 100;
% Boundary Conditions
u = zeros(M,N);
u(1,:) = To*(sin(pi/L*(0:N-1)*Dx)).^2; %temp @t = 0

A = (1+2*(1-lambda)*r)*diag(ones(N-2,1)) - (1-lambda)*r*diag(ones(N-3,1),-1) - (1-lambda)*r*diag(ones(N-3,1),1);
B = zeros(N-2,1);
for j = 1:M-1
    for i = 2:N-1
        B(i-1) = lambda*r*u(j,i+1) + (1-2*lambda*r)*u(j,i) + lambda*r*u(j,i-1);
    end
    X = A\B;
    u(j+1,2:N-1) = X';
end

max_temp = max(u(M,:));


% Contour Plots

x = 0:Dx:L;
y  = 0:Dt:t_max;
[X , Y] = meshgrid(x,y);

subplot(2,1,1);
contour(X , Y , u , 'ShowText' , 'on');
title('Temperature (Crank-Nicolson)');
xlabel('Length x');
ylabel('Time');
colorbar;

subplot(2,1,2);
contour(X , Y , u , 20);
title('Temperature (Crank-Nicolson)');
xlabel('Length x');
ylabel('Time');
colorbar;

end