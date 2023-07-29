% Validation 
clear all;
K=3;
M=32;
Nr=1;
theta=60/180*pi
for rr=1:1
r=10
R=40
lambda=3e8/100e9;
d=lambda/2;
D0=1/2*2^5*lambda;
D=d*(M-1)+D0;
 
% d=(D*(K-1)+(M-1)*d)/(K*M-1); 
% D=d*M;
DD=D/r;
Dd=d/r; 
D_array=D*(K-1)+d*(M-1);
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

S_theta_sum=0;
S_r_sum=0;
S_theta2_sum=0;
S_theta_r_sum=0;
S_r2_sum=0;
for k=-(K-1)/2:(K-1)/2
    for m=-(M-1)/2:(M-1)/2
        S_theta2_sum=S_theta2_sum+(k*D+m*d)^2/(r^2-2*(k*D+m*d)*r*sin(theta)+(k*D+m*d)^2);
        S_theta_sum=S_theta_sum+((m)*d+(k)*D)/sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2);
        S_r_sum=S_r_sum+(((m)*d+(k)*D)*sin(theta)-r)/sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2);
        S_theta_r_sum=S_theta_r_sum+(m*d+k*D)*((m*d+k*D)*sin(theta)-r)/(r^2-2*(m*d+k*D)*r*sin(theta)+(m*d+k*D)^2)
        S_r2_sum=S_r2_sum+(r^2-2*(m*d+k*D)*sin(theta)*r+(m*d+k*D)^2*sin(theta)^2)/(r^2-2*(m*d+k*D)*r*sin(theta)+(m*d+k*D)^2)
    end
end
S_theta_sum
S_r_sum

Q11=chit*(S_theta2/K/M-S_theta^2/K^2/M^2)+chir2*phi_theta^2
Q22=chit/r^2/cos(theta)^2*(S_r2/K/M-S_r^2/K^2/M^2)+chir2*phi_r^2
Q12=chit/r/cos(theta)*(S_theta_r/K/M-S_theta*S_r/K^2/M^2)+chir2*phi_theta*phi_r

Q11_sum=chit*(S_theta2_sum/K/M-S_theta_sum^2/K^2/M^2)+chir2*phi_theta^2
Q22_sum=chit/r^2/cos(theta)^2*(S_r2_sum/K/M-S_r_sum^2/K^2/M^2)+chir2*phi_r^2
Q12_sum=chit/r/cos(theta)*(S_theta_r_sum/K/M-S_theta_sum*S_r_sum/K^2/M^2)+chir2*phi_theta*phi_r


DQ=Q11*Q22-Q12^2
DQ_sum=Q11_sum*Q22_sum-Q12_sum^2

CRB_r=sqrt(Q11/DQ/(K*M*Nr)/2)
CRB_theta=sqrt(Q22/DQ/(K*M*Nr)/2)


CRB_r_sum=sqrt(Q11_sum/DQ_sum/(K*M*Nr)/2)
CRB_theta_sum=sqrt(Q22_sum/DQ_sum/(K*M*Nr)/2)
end