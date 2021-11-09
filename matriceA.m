function [Aij]=matriceA(i,j)
    global N;
    global n;
    global w;
    global ksi;
    global k;
    global s;
    Aij = 0;
    ai = s(:,i);
    bi = s(:,mod(i,N)+1);
    aj = s(:,j);
    bj = s(:,mod(j,N)+1);
    if (i~=j)
        for ind1=[1:1:n]
            x1 = ksi(ind1)*((bi-ai)/2)+(bi+ai)/2;
            for ind2=[1:1:n]
                x2 = ksi(ind2)*((bj-aj)/2)+(bj+aj)/2;
                Aij = Aij + (norm(bi-ai)*w(ind1)/2)*(norm(bj-aj)*w(ind2)/2)*(1i/4)*besselh(0,k*norm(x1-x2));
            end
        end
    else
        for ind=[1:1:n]
            x = ksi(ind)*((bi-ai)/2)+(bi+ai)/2;
            dj1=bi-x;
            dj=ai-x;
            taue=(bi-ai)/norm(bi-ai);
            I = (-1/(2*pi))*(dj1'*taue*log(norm(dj1))-dj'*taue*log(norm(dj))-norm(bi-ai))+(1i/4-(1/(2*pi))*(log(k/2)+0.5772156649))*norm(bi-ai);
            Aij = Aij + (norm(bi-ai)*w(ind)/2)*I;
        end
    end
    
end