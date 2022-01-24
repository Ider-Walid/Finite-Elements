function y = evaluation(x,xm,u,k)
    y = 0;
    for il = 1:4
        y = y + u(2*(k-1)+il)*phibase(x, il, k, xm);
    end
end