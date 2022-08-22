for N = 1:10
    %%% Get values of E for capacitance
    a = 0.005;
    b = 0.05;
    V = 1;
    c = matrix_S(N,a,b)\transpose(matrix_B(N,a,b,V));
    %c = transpose(c);

    for k = 1:101
        points(k) = a + ((b-a)/100)*(k-1);
    end
    r = points;

    f0 = (r-b)/(a-b)*V;
    for i = 1:N
        for k = 1:101
            fn(i,k) = c(i)*(r(k)-a)*(r(k)-b)*r(k)^i;
        end
    end

    F = [f0 ; fn];
    F_tot = sum(F,1);

    %%% E = - Grad F
    for i = 1:101
        E0(i) = -V/(a-b);
    end

    for i = 1:N
        for k = 1:101
            En(i,k) = -c(i)*((i+2)*r(k)^(i+1) - (i+1)*(a+b)*r(k)^i + a*b*i*r(k)^(i-1));
        end
    end

    E = [E0 ; En];
    E_tot = sum(E,1);

    for i = 1:101
        E_real(i) = 1/(log(10)*r(i));
    end
    
    %Capacitance
    C_ex = (2*pi*8.854*10^(-12))/log(10);
    for i = 1:101
        E_dens(i) = E_tot(i)^2*r(i);
    end
    W = (pi*8.854*10^(-12))*trapz(r,E_dens);
    C_est(N) = 2*W/V^2;
    error(N) = abs(C_est(N)-C_ex)/C_ex*100;
    

end

x = [1:10];
plot(x,error)
xlabel('N')
ylabel('error %')
