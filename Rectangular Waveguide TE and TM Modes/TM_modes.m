%% Wavenumber k cutoff theoretical values
k_th = [351.23 444.27 566.34 647.63]';
%%% Dimensions M, N either 50x25 or 100x50
M = 100;
N = 50;
h = 1/N/100;
A = 4*diag(ones((N-1)*(M-1),1)) - diag(ones((N-1)*(M-1)-1,1),1) - diag(ones((N-1)*(M-1)-1,1),-1) - diag(ones((N-1)*(M-1)-(M-1),1),M-1) - diag(ones((N-1)*(M-1)-(M-1),1),-(M-1));
for i = M:M-1:(M-1)*(N-1)
    A(i,i-1) = 0;
    A(i-1,i) = 0;
end
lambda = eigs(A,4,'smallestabs');
%% cut-off wavenumbers
k_cutoff = sqrt(lambda)/h;
%% errors
error = abs(k_th - k_cutoff)./k_th*100;
fprintf(' TM_11 error : %f %% \n TM_21 error : %f %% \n TM_31 error : %f %% \n TM_12 error : %f %% \n',error(1),error(2),error(3),error(4));
%% Eigenvalues, Eigenvectors & Electric Field

[V,D] = eig(A);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

E1 = reshape(Vs(:,1),M-1,N-1)';
E2 = reshape(Vs(:,2),M-1,N-1)';
E3 = reshape(Vs(:,3),M-1,N-1)';
E4 = reshape(Vs(:,4),M-1,N-1)';

E_1 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        E_1(j+1,i+1) = E1(j,i);
    end
end

E_2 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        E_2(j+1,i+1) = E2(j,i);
    end
end

E_3 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        E_3(j+1,i+1) = E3(j,i);
    end
end

E_4 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        E_4(j+1,i+1) = E4(j,i);
    end
end

%% Contour Plots

x = 0:h*100:2;
y = 0:h*100:1;
[X , Y] = meshgrid(x,y);

subplot(2,2,1);
contour(X , Y , E_1 , 20);
xlabel('x');
ylabel('y');
title('Electric field, TM 11 mode');
colorbar;

subplot(2,2,2);
contour(X , Y , E_2 , 20);
xlabel('x');
ylabel('y');
title('Electric field, TM 21 mode');
colorbar;

subplot(2,2,3);
contour(X , Y , E_3 , 20);
xlabel('x');
ylabel('y');
title('Electric field, TM 31 mode');
colorbar;

subplot(2,2,4);
contour(X , Y , E_4 , 20);
xlabel('x');
ylabel('y');
title('Electric field, TM 12 mode');
colorbar;
