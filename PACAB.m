function [U0,V0]=PACAB(ind2,ind1)
    global eps;
    k=1;
    n1 = length(ind1);
    n2 = length(ind2);
    i = 1;
    I = [];
    U0 = [];
    V0 = [];
    bool = true;
    firstStep = true;
    frobeniusAk2 = 0;
    nUk = 0;
    nVk = 0;
    vk = zeros([1,n1]);
    uk = zeros([n2,1]);
    UV = uk*vk;
    while (firstStep || (bool && frobeniusAk2*eps*eps>=nUk*nVk))
        for l=[1:n1]
            vk(l) = matriceA(ind2(i),ind1(l))-UV(i,l);
        end
        [~, j] = max(abs(vk));
        delta = vk(j);
        if (abs(delta)<=1e-4)
            if (length(U0)>=n2-1)
                bool = 0;
            end
        else
            for l=[1:n2]
                uk(l) = matriceA(ind2(l),ind1(j))-UV(l,j);
            end
            vk = vk./delta;
            k=k+1;
        end
        uktmp = uk;
        I = [I,i];
        uktmp(I) = 0;
        [~, i] = max(abs(uktmp));
        nUk = sum(uk.*uk);
        nVk = sum(vk.*vk);
        frobeniusAk2=frobeniusAk2+nUk*nVk+2*uk'*UV*vk';
        U0 = [U0 uk];
        V0 = [V0;vk];
        UV = UV + uk*vk;
        firstStep = false;
    end
end 