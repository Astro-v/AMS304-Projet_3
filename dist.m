function [d]=dist(I1, I2)
    d=0;
    d=d+(max(0,min(I1.s(1,:))-max(I2.s(1,:)))^2+max(0,min(I2.s(1,:))-max(I1.s(1,:)))^2);
    d=d+(max(0,min(I1.s(2,:))-max(I2.s(2,:)))^2+max(0,min(I2.s(2,:))-max(I1.s(2,:)))^2);
    d=sqrt(d);% changer par I.xmax ..
end