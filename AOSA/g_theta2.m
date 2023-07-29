function f = g_theta2(x,theta,DD,Dd)

f=x^2/DD/Dd/2+x*sin(theta)/DD/Dd*log(abs(1-2*sin(theta)*x+x^2))-2*x*sin(theta)/DD/Dd-sin(theta)^2/DD/Dd*log(abs(1-2*x*sin(theta)+x^2))...
    +2*cos(theta)*sin(theta)/DD/Dd*atan(x/cos(theta)-tan(theta))-cos(2*theta)/DD/Dd*(x/cos(theta)-tan(theta))*atan(x/cos(theta)-tan(theta))...
    +cos(2*theta)/2/DD/Dd*log(abs((x/cos(theta)-tan(theta))^2+1));

end
