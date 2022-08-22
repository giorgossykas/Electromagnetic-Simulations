function max_temp = Dufort_Frankel(Dx,Dt)
% Parameters
L = 100;
t_max = 1000;
%Dx = 2;
%Dt = 0.1;
r = Dt/Dx^2;
N = L/Dx + 1;
M = t_max/Dt + 1;
To = 100;

%Boundary Conditions
u = zeros(M,N);
u(1,:) = To*(sin(pi/L*(0:N-1)*Dx)).^2; %temp @t = 0

%Central Differences for t = Dt
for i = 2:N-1
    u(2,i) = r*u(1,i+1) + (1-2*r)*u(1,i) + r*u(1,i-1);
end

%Dufort-Frankel
for j = 2:M-1
    for i = 2:N-1
        u(j+1,i) = (1/(1+2*r))*(2*r*(u(j,i+1) + u(j,i-1)) + (1-2*r)*u(j-1,i));
    end
end

max_temp = max(u(M,:));

x = 0:Dx:L;
y = 0:Dt:t_max;
[X , Y] = meshgrid(x,y);

subplot(2,1,1);
contour(X , Y , u , 'ShowText' , 'on');
title('Tempetature (Dufort-Frankel)');
xlabel('Length x');
ylabel('Time');
colorbar;
subplot(2,1,2);
contour(X , Y , u , 20);
title('Tempetature (Dufort-Frankel)');
xlabel('Length x');
ylabel('Time');
colorbar;

end