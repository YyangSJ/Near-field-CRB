function [CRB_r,CRB_theta] = WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)
% 
DD=D/r;
Dd=d/r; 
phi_theta=(r*cos(theta)*(R^2+r^2-R*r*cos(theta))-2*R*r^2*sin(theta)^2)/(R^2+r^2-2*R*r*cos(theta))^(3/2);
phi_r=R*sin(theta)*(R-r*cos(theta))/(R^2+r^2-2*R*r*cos(theta))^(3/2);

chit=4*pi^2*cos(theta)^2*r^2/lambda^2;
chir2=pi^2*d^2*(Nr^2-1)/3/lambda^2

x1=-0.5*K*DD-0.5*M*Dd;x2=-0.5*K*DD+0.5*M*Dd;x3=0.5*K*DD-0.5*M*Dd;x4=0.5*K*DD+0.5*M*Dd;
S_theta2=g_theta2(x4,theta,DD,Dd)-g_theta2(x3,theta,DD,Dd)-g_theta2(x2,theta,DD,Dd)+g_theta2(x1,theta,DD,Dd)
S_theta=g_theta(x4,theta,DD,Dd)-g_theta(x3,theta,DD,Dd)-g_theta(x2,theta,DD,Dd)+g_theta(x1,theta,DD,Dd)
S_r=sin(theta)*S_theta-(g_r(x4,theta,DD,Dd)-g_r(x3,theta,DD,Dd)-g_r(x2,theta,DD,Dd)+g_r(x1,theta,DD,Dd))
S_r2=K*M-cos(theta)^2*S_theta2
S_theta_r=sin(theta)*S_theta2-(g_theta_r(x4,theta,DD,Dd)-g_theta_r(x3,theta,DD,Dd)-g_theta_r(x2,theta,DD,Dd)+g_theta_r(x1,theta,DD,Dd))

 

Q11=chit*(S_theta2/K/M-S_theta^2/K^2/M^2)+chir2*phi_theta^2
Q22=chit/r^2/cos(theta)^2*(S_r2/K/M-S_r^2/K^2/M^2)+chir2*phi_r^2
Q12=chit/r/cos(theta)*(S_theta_r/K/M-S_theta*S_r/K^2/M^2)+chir2*phi_theta*phi_r
 

DQ=Q11*Q22-Q12^2 

CRB_r=sqrt(Q11/DQ/(K*M*Nr)/2)
CRB_theta=sqrt(Q22/DQ/(K*M*Nr)/2)
 
end

