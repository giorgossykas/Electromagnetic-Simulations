for N = [2 5 10]
    a = 0.002;
    b = 0.02;
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

    %%%%%%Plots

    plot(r,E_tot)
    xlabel('Radial Distance')
    ylabel('Electric Field')
    hold on
end

plot(r,E_real,'r')
legend({'N=2','N=5','N=10','Real'},'Location','northeast')
hold off

