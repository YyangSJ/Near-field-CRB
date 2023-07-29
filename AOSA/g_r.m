function f = g_r(x,theta,DD,Dd)

f=1/DD/Dd*(x-sin(theta))*log(abs(sqrt(x^2-2*x*sin(theta)+1)+x-sin(theta)))-1/DD/Dd*sqrt(x^2-2*x*sin(theta)+1);
end

