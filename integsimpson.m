function val = integsimpson(fun,a,b,varargin)
c = (a+b) / 2;
y = feval(fun,[a c b],varargin{:});
fa = y(1);
fc = y(2);
fb = y(3);
val = (b-a) * (fa+4*fc+fb) / 6;
end