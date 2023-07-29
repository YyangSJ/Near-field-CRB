function f = g_theta(x,theta,DD,Dd)

f=(x-sin(theta))/2/DD/Dd*sqrt(1-2*sin(theta)*x+x^2)+cos(theta)^2/2/DD/Dd*log(abs(sqrt(1-2*sin(theta)*x+x^2)-sin(theta)+x))...
    +sin(theta)*(x-sin(theta))/DD/Dd*atanh((x-sin(theta))/sqrt(1-2*sin(theta)*x+x^2))-sin(theta)/DD/Dd*sqrt(1-2*sin(theta)*x+x^2);


end

