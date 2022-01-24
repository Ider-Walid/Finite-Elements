%% TRAVAIL RÉALISÉ PAR IDER WALID - ZZ2 F4

%% TP FINAL EF


% Programme de calcul de la solution du probleme variationnel
% Méthodes des EF P3


% Contexte du problème:


% Résolution sur l'intervalle [0,1] de :
% -beta*u'' + alpha*u = f
% ici : Alpha = 0 et Beta= 1 ( donc -u''=f)
% CL de Dirichlet : u(a) = -1 & u(b) = 1



%  Initialisations de nos variables:

beta = 1;
alpha = 0; 
a = 0;
b = 1;
if a>=b, error('On a besoin de a<b'),end
n = 2;
h = (b-a)/n; % Initialisation du Pas
xm = a:h:b; %  Maillage dans a,b




% Assemblage du membre de droite :




F = sparse(2*n+2,1); % Initialisation
for k = 1:n % Considération de l'intervalle[xm(k),xm(k+1)]
    for il = 1:4 
        ig = 2*(k-1) + il; 
        F(ig) = F(ig) + integsimpson('phibasef',xm(k) ,xm(k+1),il,k,xm);
    end
end


% Assemblage des matrices :


K = sparse(2*(n+1),2*(n+1)); % Initialisation
M = K; % Initialisation
for k = 1:n % Considération de l'intervalle [xm(k),xm(k+1)] 
    for il =1:4 
        ig = 2*(k-1) + il; 
        for jl = 1:4 
            jg = 2*(k-1) + jl; 
            M(ig,jg) = M(ig,jg) + integsimpson('phiphibase',xm(k),xm(k+1),il,jl,k,xm);
            K(ig,jg) = K(ig,jg) + integsimpson('phiderphiderbase',xm(k),xm(k+1),il,jl,k,xm);
        end
    end
end

K = beta*K + alpha*M; 



% Matrice du membre de gauche




%% Init
ca =1;
cb =-1;
  
for i=2:4
  F(i) = F(i) - K(1,i)*ca;
end
for i = 2*n-1:2*n+2
  if (i~=2*n+1)
    F(i)=F(i) - K(2*n+1,i)*cb;
  end
end

Ff = F([2:2*n 2*n+2]);
Kf=K([2:2*n 2*n+2],[2:2*n 2*n+2]);

uh = Kf\Ff;

u = zeros(length(uh)+2,1);
u(2:2*n+1) = uh;
u(1)=ca;
u(2*n+1)=cb;    




X = [];% INIT
Y = []; % INIT
u = [1;u(2:2*n+1);-1];
for k = 1:n 
    x = xm(k):0.01:xm(k + 1);
    y = evaluation(x,xm,u,k);
    X =horzcat(X,x);%% Horizontal concatenation of arrays (MATLAB)
    Y =horzcat(Y,y);
end
plot(X,Y,'m*');
hold on;





%% Fonction théorique




c1 = -3/2;
c2 = 1;
XXX = 0:0.01:1;
YYY = -0.5*XXX.*XXX + c1*XXX + ones(size(XXX)); % Fonction réelle
plot(XXX,YYY,'b');

%% +++ 
legend('Fonction approchée', 'Fonction théorique')
title('Méthode des EF P3 pour les CL de Dirichlet -> Alpha =0')
