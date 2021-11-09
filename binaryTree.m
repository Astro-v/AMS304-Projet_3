function [L]=binaryTree(s,index,Nleaf,level)
    I.xmax = max(s(1,:));
    I.ymax = max(s(2,:));
    I.xmin = min(s(1,:));
    I.ymin = min(s(2,:));
    I.diam = sqrt((I.xmax-I.xmin)^2+(I.ymax-I.ymin)^2);
    I.level = level;
    I.index = index;
    I.s = s;
    L = [I];
    if (length(s)>Nleaf)
        L = [binaryTree(s(:,1:floor(length(s)/2)),index(1:floor(length(s)/2)),Nleaf,level+1) L binaryTree(s(:,(floor(length(s)/2)+1):end),index((floor(length(s)/2)+1):end),Nleaf,level+1)];
    end
end