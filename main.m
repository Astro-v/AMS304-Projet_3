%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                        TP3 AMS304                        %%%%%%%%%
%%%%%%%%% Accélération du code BEM pour l’équation de Helmholtz 2D %%%%%%%%%
%%%%%%%%%               Farah CHAABAN & Valentin MICHEL             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear();

profile on;

%% 0. Initialisation
global k;k = 2*pi;
global eps; eps=1e-5;
global eta; eta=3;
global R;R = 1;
global dens;dens = 10;
global N; N = floor(dens*2*k*R+1)

global cas; cas = "circle";
global display; display = false;

global s; [s,c] = mesh();

global N; N = length(s);

% Quadrature
global w; w = [128/225;(322+13*sqrt(70))/900;(322+13*sqrt(70))/900;(322-13*sqrt(70))/900;(322-13*sqrt(70))/900];
global ksi; ksi = [0;+1/3*sqrt(5-2*sqrt(10/7));-1/3*sqrt(5-2*sqrt(10/7));+1/3*sqrt(5+2*sqrt(10/7));-1/3*sqrt(5+2*sqrt(10/7))];
global n; n=length(ksi);

%% 1. ACA

%% 2. PACA

%% 3. Binary tree

L=binaryTree(s,[1:N],10,0);

%% 4. admissibilité

A = admi(L,L);

%% Seconds Membre

b = B(s);

%% trace de p

p = gmres(@(U)produit(A,U),b,10,1e-06,N);

figure()
plot(real(-1i.*log((c(1,:)+1i.*c(2,:))./R)),real(p),real(-1i.*log((c(1,:)+1i.*c(2,:))./R)),imag(p))
legend({'Partie réel de la trace de p','Partie réel de la trace de p'})
xlabel('\theta')
ylabel('Trace de p')


%% DISPLAY SOLUTION SQUARE

if (display && strcmp(cas,"circle"))
    x = [-5:0.1:5];
    y = [-5:0.1:5];
    sol=[];
    for m = [1:1:size(x,2)]
        solb = []; 100*m/size(x,2)
        for l = [1:1:size(y,2)]
            if (x(m)^2+y(l)^2>=1)
                g = resolve(s,p,[x(m),y(l)]);
                solb = [solb (real(g)+1)/2];
            else
                solb = [solb 0];
            end
        end
        sol = [sol;solb];
    end

    figure()
    imshow(sol)
end

if (display && strcmp(cas,"square"))
    x = [-5:0.05:5];
    y = [-5:0.05:5];
    sol=[];
    for m = [1:1:size(x,2)]
        solb = []; 100*m/size(x,2)
        for l = [1:1:size(y,2)]
            if (max(abs(x(m)),abs(y(l)))>=R)
                g = resolve(s,p,[x(m),y(l)]);
                solb = [solb (real(g)+1)/2];
            else
                solb = [solb 0];
            end
        end
        sol = [sol;solb];
    end

    figure()
    imshow(sol)
end

profile viewer

%% TIME COMPLEXITY

NBR = [26,51,101,202,403];
T_MATA = [0.107,0.379,1.489,5.104,21.080];
