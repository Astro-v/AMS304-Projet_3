function [block,n]=admi(L1,L2)
    global eta;
    global validation;
    block.I1 = L1((length(L1)+1)/2);
    block.I2 = L2((length(L2)+1)/2);
    if (min(block.I1.diam,block.I2.diam)<eta*dist(block.I1,block.I2))
        block.adm = true;
        block.leaf = true;
        [block.U0,block.V0] = PACAB(block.I1.index,block.I2.index);
        n = size(block.U0,1)*size(block.U0,2)+size(block.V0,1)*size(block.V0,2);
        if (validation)
        	patch([block.I1.index(1) block.I1.index(end) block.I1.index(end) block.I1.index(1)],-[block.I2.index(1) block.I2.index(1) block.I2.index(end) block.I2.index(end)],"green");
        	text((block.I1.index(1)+block.I1.index(end))/2,-(block.I2.index(1)+block.I2.index(end))/2,int2str(rank(block.U0*block.V0)));
        end
    elseif (length(L1)>1)
        block.adm = false;
        block.leaf = false;
        L11 = L1(1:(length(L1)-1)/2);
        L12 = L1((length(L1)+3)/2:end);
        L21 = L2(1:(length(L2)-1)/2);
        L22 = L2((length(L2)+3)/2:end);
        [block.child1,n1] = admi(L11,L21);
        [block.child2,n2] = admi(L11,L22);
        [block.child3,n3] = admi(L12,L22);
        [block.child4,n4] = admi(L12,L21);
        n = n1+n2+n3+n4;
    else
        block.adm = false;
        block.leaf = true;
        block.A = zeros(length(block.I1.index),length(block.I2.index));
        for i = [1:length(block.I1.index)]
            for j = [1:length(block.I2.index)]
                block.A(i,j) = matriceA(block.I1.index(i),block.I2.index(j));
            end
        end
        n = size(block.A,1)*size(block.A,2);
        if (validation)
            patch([block.I1.index(1) block.I1.index(end) block.I1.index(end) block.I1.index(1)],-[block.I2.index(1) block.I2.index(1) block.I2.index(end) block.I2.index(end)],"red");
        end
    end
end