function B = matrix_B(N,a,b,V)

for k = 1:101
    points(k) = a + ((b-a)/100)*(k-1);
end
r = points;

for i = 1:N
    for k = 1:101
        Y(k) = (r(k)-a)*(r(k)-b)*r(k)^i;
    end
    B(i) = ((-2*pi*V)/(a-b))*trapz(points,Y);
end