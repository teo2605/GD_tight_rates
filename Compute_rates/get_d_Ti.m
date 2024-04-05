function dTi_dh = get_d_Ti(k,h,i)    
    dTi_dh = -get_Ti(k,h,i)/h + 2*i/h * ( (1-k*h)^(-2*i-1) - (1-h)^(-2*i-1)  );
end