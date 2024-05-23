function y =LTP_D_func(temp,a,b,c,d,e,f,g,lam,p_d)
if p_d == 1
    y = (a*tanh((temp/b)-c)).*(temp>=b*c);
else
    y = d*exp(-(temp+e)/lam) + f*temp - g;
end