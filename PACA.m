function [Ar]=PACA(A)
    global eps;
    k=1;
    R=A;
    i=1;
    U0 = [];
    V0 = [];
    bool = 1;
    while (bool)
        [~, j] = max(abs(R(i,:)));
        delta = R(i,j);
        if (abs(delta)<=1e-5)
            if (length(U0)==length(A)-1)
                bool = 0;
            end
        else
            uk = R(:,j); 
            vk = R(i,:)./delta;
            R = R-uk*vk;
            k=k+1;
        end
    U0 = [U0 uk];
    V0 = [V0;vk];
    uk(i)=0;
    [~, i] = max(abs(uk));
    end
    Ar = U0*V0-A;
end 