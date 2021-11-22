function [u]=resolve(s,p,x)
    global k;
    global N;
    global w;
    global ksi;
    global n;
    hank = zeros(1,N);
    for i=[1:1:N]
        a = s(:,i);
        b = s(:,mod(i,N)+1);
        for ind=[1:1:n]
            y = ksi(ind)*((b-a)/2)+(b+a)/2;
            hank(i)=hank(i)+(norm(b-a)*w(ind)/2.0)*(1i/4)*besselh(0,k*norm(x-y));
        end
    end
    u = hank*p;
end