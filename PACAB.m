function [U0,V0]=PACAB(ind1,ind2)
    global eps;
    k=1;
    n1 = length(ind1);
    n2 = length(ind2);
    i = 1;
    I = [];
    U0 = [];
    V0 = [];
    bool = true;
    frobeniusAk2 = 0;
    nUk = 0;
    nVk = 0;
    vk = zeros([1,n2]);
    uk = zeros([n1,1]);
    UV = uk*vk;

    while (bool && frobeniusAk2*eps<=nUk*nVk)
        for l=[1:n2]
            vk(l) = matriceA(ind1(i),ind2(l))-UV(i,l);
        end
        [~, j] = max(abs(vk));
        delta = vk(j); 
        if (abs(delta)<=1e-12)
            if (length(I)>=n2-2)
                bool = 0;
            end
        else
            for l=[1:n1]
                uk(l) = matriceA(ind1(l),ind2(j))-UV(l,j);
            end
            uk = uk./delta;
            vk = vk;
            k=k+1;
            U0 = [U0 uk];
            V0 = [V0;vk];
            UV = UV + uk*vk;
        end
        uktmp = uk;
        I = [I,i];
        uktmp(I) = 0;
        [~, i] = max(abs(uktmp)); 
        nUk = uk'*uk;
        nVk = vk*vk';
        frobeniusAk2=frobeniusAk2+nUk*nVk+2*uk.'*UV*vk.';

    end
end 