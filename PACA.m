function [U0,V0]=PACA(A)
    global eps;
    k=1;
    R=A;
    i=1;
    I = [i];
    U0 = [];
    V0 = [];
    bool = 1;
    frobeniusR=frobenius(R);
    while (bool && frobenius(R)>=eps*frobenius(A))
        [~, j] = max(abs(R(i,:)));
        delta = R(i,j);
        if (abs(delta)<=1e-4)
            if (length(U0)>=length(A)-1)
                bool = 0;
            end
        else
            uk = R(:,j); 
            vk = R(i,:)./delta;
            R = R-uk*vk;
            k=k+1;
            U0 = [U0 uk];
            V0 = [V0;vk];
        end
        utmp = uk;
        utmp(I)=0;
        [~, i] = max(abs(utmp));
        I = [I,i];
        frobeniusR=sqrt(frobeniusR^2+(norm(uk)*norm(vk))^2+2*(uk'*U0)*V0*vk');
    end
end 