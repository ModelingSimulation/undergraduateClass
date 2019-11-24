function vv = vel_car(d)
global d_min d_max v_max;
if (d < d_min)
  vv=0;
elseif (d < d_max) 
  vv=v_max*log(d/d_min)/log(d_max/d_min);
else
  vv=v_max;
end

