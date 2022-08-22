function max_temp = Euler(Dx , Dt)
% Parameters
%Dx = 0.5;
%Dt = 0.1;
L = 100;
t_max = 1000;
To = 100;
r = Dt/Dx^2;
N = L/Dx + 1;
M = t_max/Dt + 1;

% Known Boundary Conditions
u = zeros(N,M); %Later transpose
for i = 1:N
    u(i,1) = To*(sin(pi/L*(i-1)*Dx))^2; %temp @t=0
end

for j = 1:M-1
    for i = 2:N-1
        u(i,j+1) = r*u(i+1,j) + (1-2*r)*u(i,j) + r*u(i-1,j);
    end
end
u = u';

max_temp = max(u(M,:));

%Contour Plots
x = 0:Dx:L;
y = 0:Dt:t_max;
[X , Y] = meshgrid(x,y);

subplot(2,1,1);
contour(X , Y , u , 'ShowText' , 'on')
title('Temperature (Euler)');
xlabel('Length x');
ylabel('Time');
colorbar;

subplot(2,1,2);
contour(X , Y , u , 20)
title('Temperature(Euler)');
xlabel('Length x');
ylabel('Time');
colorbar;

end