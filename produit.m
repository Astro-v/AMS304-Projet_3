function [V] = produit(A,U)
    V = zeros(size(U));
    N = length(U);
    if (A.leaf == false)
        V = V + produit(A.child1,U) + produit(A.child2,U) + produit(A.child3,U) + produit(A.child4,U);
    else if (A.leaf)
        if (A.adm)
            V(A.I1.index) = A.U0*A.V0*U(A.I2.index);
        else
            V(A.I1.index) = A.A*U(A.I2.index);
        end
    end
end