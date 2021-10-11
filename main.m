%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                        TP3 AMS304                        %%%%%%%%%
%%%%%%%%% Accélération du code BEM pour l’équation de Helmholtz 2D %%%%%%%%%
%%%%%%%%%               Farah CHAABA & Valentin MICHEL             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear();

%% 0. SVD
n=400;
err = [];
for N=[n/2:n]
    N
    A = rand(n);
    [U,D,V] = svd(A);
    D= diag(D);
    D = [D(1:N);zeros(n-N,1)];
    Abis = U*diag(D)*V.';
    err = [err;sum(sum(A.*A-Abis.*Abis))];
end
figure()
plot([n/2:n],err);