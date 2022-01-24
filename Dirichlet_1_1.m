%% TRAVAIL RÉALISÉ PAR IDER WALID - ZZ2 F4

%% TP EF FINAL (Semaine 4)


% Programme de calcul de la solution du probleme variationnel
% Méthodes des EF P3


% Contexte du problème:


% Résolution sur l'intervalle [0,1] de :
% -beta*u'' + alpha*u = f
% ici : Alpha = Beta= 1 (donc -u''+u=f)
% CL de Dirichlet : u(a) = -1 & u(b) = 1



%  Initialisation de variables:

beta = 1;
alpha = 1; 
a = 0;
b = 1;
if a>=b, error('On a besoin de a<b'),end
n = 2;
h = (b-a)/n; % Initialisation du Pas
xm = a:h:b; %  Maillage dans a,b




% Assemblage du membre de droite :




F = sparse(2*n+2,1); % Initialisation
for k = 1:n % Considération de l'intervalle [xm(k),xm(k+1)] 
    for il = 1:4 
        ig = 2*(k-1) + il; 
        F(ig) = F(ig) + integsimpson('phibasef',xm(k) ,xm(k+1),il,k,xm);
    end
end





% Assemblage de nos matrices :



K = sparse(2*(n+1),2*(n+1)); % Initialisation
M = K; % Initialisation
for k = 1:n % Considération de l'intervalle[xm(k),xm(k+1)] 
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

X = [];
Y = [];
u = [1;u(2:2*n+1);-1];
for k = 1:n 
    x = xm(k):0.01:xm(k + 1);
    y = evaluation(x,xm,u,k);
    X =horzcat(X,x);
    Y =horzcat(Y,y);
end
plot(X,Y,'m*');
hold on;




%% Fonction théorique:




c1 = 0;
c2 = -2/sinh(1);
XXX = 0:0.01:1;
YYY = c1*cosh(XXX) + c2*sinh(XXX) + ones(size(XXX)); % FONCTON RÉELLE
plot(XXX,YYY,'b');
%axis equal


% +++ 
title('Méthode des EF P3 pour les CL de Dirichlet -> Alpha =1')
legend('Fonction approchée', 'Fonction théorique')
