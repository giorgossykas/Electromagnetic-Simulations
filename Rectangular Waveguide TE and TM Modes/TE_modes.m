%% Theoretical values
k_th = [0 157.08 314.15 314.15 351.23]';
%%% Dimensions M, N either 50x25 or 100x50
M = 100;
N = 50;
h = 1/N/100;
%% matrix A
A = 4*diag(ones((M-1)*(N-1),1)) - diag(ones((M-1)*(N-1)-1,1),1) - diag(ones((M-1)*(N-1)-1,1),-1) - diag(ones((M-1)*(N-1)-(M-1),1),M-1) - diag(ones((M-1)*(N-1)-(M-1),1),-(M-1));
%% zeros
for i = M:M-1:(M-1)*(N-1)
    A(i-1,i) = 0;
    A(i,i-1) = 0;    
end
%% 2s and first last 3s
for i = 1:(M-1)*(N-1)
    if i == 1 || i == M-1 || i == (M-1)*(N-2)+1 || i == (M-1)*(N-1)
        A(i,i) = 2;
    end
    
    if (i > 1 && i < M-1) || (i > (M-1)*(N-2)+1 && i < (M-1)*(N-1)) 
        A(i,i) = 3;
    end
end
% in between 3s
for i = M:M-1:(M-1)*(N-3)+1
    A(i,i) = 3;
    A(i+(M-2),i+(M-2)) = 3;
end
lambda = eigs(A,5,'smallestabs'); %% first 5 because 1st is zero(false).
%% cut-off wavenumbers
k_cutoff = sqrt(lambda)/h;
%% errors
error = abs(k_th - k_cutoff)./k_th*100;
fprintf(' TE_10 error : %f %% \n TE_01 error : %f %% \n TE_20 error : %f %% \n TE_11 error : %f %% \n',error(2),error(3),error(4),error(5));
%% Eigenvalues, Eigenvectors and Magnetic Field
[V,D] = eig(A);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

H1 = reshape(V(:,2),M-1,N-1)';
H2 = reshape(V(:,3),M-1,N-1)';
H3 = reshape(V(:,4),M-1,N-1)';
H4 = reshape(V(:,5),M-1,N-1)';
%%
H_1 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        H_1(j+1,i+1) = H1(j,i);
    end
end
for j = 2:M
    H_1(1,j) = H_1(2,j);
    H_1(N+1,j) = H_1(N,j);
end
for i = 1:N+1
    H_1(i,1) = H_1(i,2);
    H_1(i,M+1) = H_1(i,M);
end

%%
H_2 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        H_2(j+1,i+1) = H2(j,i);
    end
end
for j = 2:M
    H_2(1,j) = H_2(2,j);
    H_2(N+1,j) = H_2(N,j);
end
for i = 1:N+1
    H_2(i,1) = H_2(i,2);
    H_2(i,M+1) = H_2(i,M);
end
%%
H_3 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        H_3(j+1,i+1) = H3(j,i);
    end
end
for j = 2:M
    H_3(1,j) = H_3(2,j);
    H_3(N+1,j) = H_3(N,j);
end
for i = 1:N+1
    H_3(i,1) = H_3(i,2);
    H_3(i,M+1) = H_3(i,M);
end
%%
H_4 = zeros(N+1,M+1);
for i = 1:M-1
    for j = 1:N-1
        H_4(j+1,i+1) = H4(j,i);
    end
end
for j = 2:M
    H_4(1,j) = H_4(2,j);
    H_4(N+1,j) = H_4(N,j);
end
for i = 1:N+1
    H_4(i,1) = H_4(i,2);
    H_4(i,M+1) = H_4(i,M);
end

%% Contour Plots
x = 0:h*100:2;
y = 0:h*100:1;
[X,Y] = meshgrid(x,y);

subplot(2,2,1);
contour(X , Y , H_1 , 30);
xlabel('x');
ylabel('y');
title('Magnetic field, TE 10 mode');
colorbar;

subplot(2,2,3);
contour(X , Y , H_2 , 30);
xlabel('x');
ylabel('y');
title('Magnetic field, TE 20 mode');
colorbar;

subplot(2,2,2);
contour(X , Y , H_3 , 30);
xlabel('x');
ylabel('y');
title('Magnetic field, TE 01 mode');
colorbar;

subplot(2,2,4);
contour(X , Y , H_4 , 30);
xlabel('x');
ylabel('y');
title('Magnetic field, TE 11 mode');
colorbar;
