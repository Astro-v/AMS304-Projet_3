function [V] = produit(A,U)
    V = zeros(size(U));
    N = length(U);
    if (A.leaf == false)
        V = V + produit(A.child1) + produit(A.child2) + produit(A.child3) + produit(A.child4);
    else if (A.leaf)
    end
end