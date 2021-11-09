function [a,c]=mesh()
    global N;
    global R;
    global cas;
    if (strcmp(cas,"circle"))
        pas=2*pi/N;
        theta = [-pi:pas:pi];
        x = R*cos(theta);
        y = R*sin(theta);
        x = x(1:length(x)-1);
        y = y(1:length(y)-1);

        a = [x; y];
        b = [x(2:N); y(2:N)];
        b(:,N)=[x(1),y(1)];

        %milieux des éléments
        c(1,1:N)= a(1,:)+(b(1,:)-a(1,:))/2;
        c(2,1:N)= a(2,:)+(b(2,:)-a(2,:))/2;
    elseif (strcmp(cas,"square"))
        pas=8/N;
        n=length([-1:pas:1])-1;
        
        a = [ones(1,n) [(1-pas):-pas:-1] -ones(1,n) [(-1+pas):pas:1];[(-1+pas):pas:1] ones(1,n) [(1-pas):-pas:-1] -ones(1,n)];
        b = [a(:,2:end) a(:,1)];
        %milieux des éléments
        c(1,1:length(a))= a(1,:)+(b(1,:)-a(1,:))/2;
        c(2,1:length(a))= a(2,:)+(b(2,:)-a(2,:))/2;
    end
end