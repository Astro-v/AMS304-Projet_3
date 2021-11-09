%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                        TP3 AMS304                        %%%%%%%%%
%%%%%%%%% Accélération du code BEM pour l’équation de Helmholtz 2D %%%%%%%%%
%%%%%%%%%               Farah CHAABAN & Valentin MICHEL             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear();

%% 0. Initialisation
global k;k = 2*pi;
global eps; eps=1e-5;
global eta; eta=3;
global R;R = 1;
global dens;dens = 20;
global N; N = floor(dens*2*k*R+1);

global cas; cas = "circle";


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