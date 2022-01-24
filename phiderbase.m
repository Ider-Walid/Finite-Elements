function y = phiderbase(x,il,k,xm)
% les derivees phi_i' des deux fonctions de base, il=1,2,3,4 sur l'intervalle [a,b]

a=xm(k);
b=xm(k+1);

if a>=b, error('a < b !!'),end
if (il<1 || il>4),error('il n''y a que deux fonctions de base'),end
if (il==1)
    y= (6*(x-b).*(x-a))./(b-a).^3;
end
if (il==2)
    y= ((x-b).*(3*x-2*a-b))./(a-b).^2 ;
end
if (il==3)
    y=  (6*(x-a)./(a-b).^3).*(x-b);
end
if (il==4)
    y=  ((x-a).*(3*x-2*b-a))./(b-a).^2;
end
end
