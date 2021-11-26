function [p]=traceAnalytique(x,y)
    global R;
    global k;
    r = sqrt(x^2+y^2); 
    theta = -i*log((x+i*y)/r);
    n = [-20:1:20];
    p = (k/2)*(-i).^(n).*(besselj(n,k*R)./besselh(n,k*R)).*(besselh(n-1,k*r)-besselh(n+1,k*r))*exp(i*theta*n).'+i*k*cos(theta)*exp(-i*k*r*cos(theta));
end