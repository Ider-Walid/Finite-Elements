function y = phiphibase(x,il,jl,k,xm)
y = phibase(x,il,k,xm) .* phibase(x,jl,k,xm);
end