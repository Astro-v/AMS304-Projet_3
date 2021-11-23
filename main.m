%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                        TP3 AMS304                        %%%%%%%%%
%%%%%%%%% Accélération du code BEM pour l’équation de Helmholtz 2D %%%%%%%%%
%%%%%%%%%               Farah CHAABAN & Valentin MICHEL             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear();

profile on;

%% 0. Initialisation
global k;k = 2*pi;
global eps; eps=1e-6;
global eta; eta=3;
global R;R = 1;
global dens;dens = 10;
global N; N = floor(dens*2*k*R+1)

global cas; cas = "circle";
global display; display = false;
global validation; validation = false;
global complexite; complexite = true;

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

if (validation)
    figure()
    plot(real(-1i.*log((c(1,:)+1i.*c(2,:))./R)),real(p),real(-1i.*log((c(1,:)+1i.*c(2,:))./R)),imag(p))
    legend({'Partie réel de la trace de p','Partie réel de la trace de p'})
    xlabel('\theta')
    ylabel('Trace de p')
end


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

% time complexity for admi
NBR = [26,51,101,202,403];
T_MATA = [0.199,0.487,1.312,3.174,7.487];
T_MATA = T_MATA./log(NBR);
P = polyfit(log(NBR),log(T_MATA),1);P(1)
if (complexite)
    figure()
    loglog(NBR,T_MATA)
    xlabel('Nombre de points (log N)')
    ylabel('Temps d execution (log (T/log N))')
end

% time complexity for ACA

NBR = [100,200,400,800]
T_ACA = [0.029,0.153,0.967,8.828]./log(NBR);
P = polyfit(log(NBR),log(T_ACA),1);P(1)
if (complexite)
    figure()
    loglog(NBR,T_ACA)
    xlabel('Nombre de points (log N)')
    ylabel('Temps d execution (log T)')
end

