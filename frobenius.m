function [norme]=frobenius(M)
    norme = sqrt(sum(sum(M.*M)));
end