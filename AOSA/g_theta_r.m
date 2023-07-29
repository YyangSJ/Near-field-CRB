function f= g_theta_r(x,theta,DD,Dd)

f=sin(theta)/DD/Dd*((x/cos(theta)-tan(theta))*atan(x/cos(theta)-tan(theta))-0.5*log(abs((x/cos(theta)-tan(theta))^2+1)))...
    +x/2/DD/Dd*log(abs(1-2*x*sin(theta)+x^2))-x/DD/Dd-sin(theta)/2/DD/Dd*log(abs(1-2*x*sin(theta)+x^2))+...
    cos(theta)/DD/Dd*atan(x/cos(theta)-tan(theta));
end

