function [y]=phibase(x,il,k,xm)
% Fonctions de base EF P3 (x4)

% Constantes:
a=xm(k);
b=xm(k+1);


if (il<1 || il>4)
  error('Erreur!')
end
if (il==1)
    y= (((x-b)./(a-b)).^2).*(2*(x-a)./(b-a)+1);
end
if (il==2)         
    y= (((x-b)./(a-b)).^2).*(x-a);
end
if (il==3)
    y=(((x-a)./(b-a)).^2).*(2*(x-b)./(a-b)+1);
end
if (il==4)
    y= (((x-a)./(b-a)).^2).*(x-b);
end
end



