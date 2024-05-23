function y =dopa_func(temp,thresh_size,d_0)
y = (temp-(d_0+thresh_size/2)).*(temp>=(d_0+thresh_size/2))...
    + (temp -(d_0-thresh_size/2)).*(temp<=(d_0-thresh_size/2)) ;
end