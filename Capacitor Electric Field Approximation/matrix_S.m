function S = matrix_S(N,a,b)

for k = 1:101
    points(k) = a + ((b-a)/100)*(k-1);
end
r = points;

for i = 1:N
    for j = 1:N
       for k = 1:101
           Y(k) = (r(k)-a)*(r(k)-b)*r(k)^(i+1)*((j+2)^2*r(k)^j - (j+1)^2*(a+b)*r(k)^(j-1) + j^2*a*b*r(k)^(j-2));
       end
       S(i,j) = 2*pi*trapz(points,Y);
    end
end
