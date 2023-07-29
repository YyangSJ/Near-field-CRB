function [CRB_r,CRB_theta] = WSMS_HSPW2(theta,lambda,r,R,K,M,Nr,D,d)
DD=D/r; 
phi_theta=(r*cos(theta)*(R^2+r^2-R*r*cos(theta))-2*R*r^2*sin(theta)^2)/(R^2+r^2-2*R*r*cos(theta))^(3/2);
phi_r=R*sin(theta)*(R-r*cos(theta))/(R^2+r^2-2*R*r*cos(theta))^(3/2);

chiK=4*pi^2*cos(theta)^2*r^2/K/lambda^2;
chiM=pi^2*d^2/3/lambda^2*(M^2-1);
chir2=pi^2*d^2*(Nr^2-1)/3/lambda^2;

T_S_theta2=K+sin(theta)/DD*log(abs((1-K*DD*sin(theta)+K^2*DD^2/4)/(1+K*DD*sin(theta)+K^2*DD^2/4)))...
    -cos(2*theta)/DD/cos(theta)*atan(K*DD/2/cos(theta)-tan(theta))+cos(2*theta)/DD/cos(theta)*atan(-K*DD/2/cos(theta)-tan(theta));

T_S_theta=sin(theta)/DD*atanh((K*DD/2-sin(theta))/sqrt(K^2*DD^2/4-K*DD*sin(theta)+1))-sqrt(K^2*DD^2/4+K*DD*sin(theta)+1)/DD...
    +sqrt(K^2*DD^2/4-K*DD*sin(theta)+1)/DD-sin(theta)/DD*atanh((-K*DD/2-sin(theta))/sqrt(K^2*DD^2/4+K*DD*sin(theta)+1));

T_S_r2=K-cos(theta)^2*T_S_theta2;

T_S_r=sin(theta)*T_S_theta-1/DD*log(abs((sqrt(K^2*DD^2/4-K*DD*sin(theta)+1)+K*DD/2-sin(theta))/(sqrt(K^2*DD^2/4+K*DD*sin(theta)+1)-K*DD/2-sin(theta))));

T_S_theta_r=sin(theta)*T_S_theta2-tan(theta)/DD*atan(K*DD/2/cos(theta)-tan(theta))+tan(theta)/DD*atan(-K*DD/2/cos(theta)-tan(theta))...
    -1/2/DD*log(abs((K^2*DD^2/4-K*DD*sin(theta)+1)/(K^2*DD^2/4+K*DD*sin(theta)+1)));
 
Q11=chiK*(T_S_theta2-T_S_theta^2/K)+chiM*cos(theta)^2+chir2*phi_theta^2;
Q22=chiK/r^2/cos(theta)^2*(T_S_r2-T_S_r^2/K)+chir2*phi_r^2;
Q12=chiK/r/cos(theta)*(T_S_theta_r-T_S_theta*T_S_r/K)+chir2*phi_theta*phi_r;



DQ=Q11*Q22-Q12^2
CRB_r=sqrt(Q11/DQ/(K*M*Nr)/2)
CRB_theta=sqrt(Q22/DQ/(K*M*Nr)/2)


end

