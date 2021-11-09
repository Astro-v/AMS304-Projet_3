function [U0,V0]=ACA(A)
    global eps;
    k=0;
    R=A;
    U0 = [];
    V0 = [];
    while (frobenius(R)>=eps*frobenius(A))
        k=k+1;
        [val, i] = max(abs(R));
        [~, j] = max(val);
        i = i(j);
        delta = R(i,j);
        uk = R(:,j); 
        vk = R(i,:)./delta;
        R = R-uk*vk;
        U0 = [U0 uk];
        V0 = [V0;vk];
    end
end 