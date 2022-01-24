function y = phiderphiderbase(x,il,jl,k,xm)
y = phiderbase(x,il,k,xm).*phiderbase(x,jl,k,xm);
end