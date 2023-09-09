function y = phi_4(temp,a,b,c)
y = (a*tanh((temp/b)-c)).*(temp>=b*c);
        
end