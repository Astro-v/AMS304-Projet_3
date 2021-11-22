function [secondMembre] = B(s)
    global N;
    global n;
    global w;
    global ksi;
    global k;
    secondMembre = zeros(N,1);
    for i=[1:1:N]
        a = s(:,i);
        b = s(:,mod(i,N)+1);
        for ind=[1:1:n]
            x = ksi(ind)*((b-a)/2)+(b+a)/2;
            secondMembre(i) = secondMembre(i) - (norm(b-a)*w(ind)/2.0)*exp(-1i*k*x(1));
        end
    end
end