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

DD=D/r;
Dd=d/r;
D_array=D*(K-1)+d*(M-1);
phi_theta=(r*cos(theta)*(R^2+r^2-R*r*cos(theta))-2*R*r^2*sin(theta)^2)/(R^2+r^2-2*R*r*cos(theta))^(3/2);
phi_r=R*sin(theta)*(R-r*cos(theta))/(R^2+r^2-2*R*r*cos(theta))^(3/2);
for K=K
    b=[];b_theta=[];b_r=[];
    kk=0;
    b_theta_sum=0;
for k=-(K-1)/2:(K-1)/2
    kk=kk+1;
    mm=0;
    for m=-(M-1)/2:(M-1)/2
        mm=mm+1;
    b(kk,mm)=sqrt(1/K/M)*exp(-1i*2*pi/lambda*sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2));
    b_theta(kk,mm)=sqrt(1/K/M)*1i*2*pi/lambda*exp(-1i*2*pi/lambda*sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2))...
        *((m)*d+(k)*D)*r*cos(theta)/sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2);
    b_r(kk,mm)=sqrt(1/K/M)*1i*2*pi/lambda*exp(-1i*2*pi/lambda*sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2))...
        *(((m)*d+(k)*D)*sin(theta)-r)/sqrt(r^2-2*((m)*d+(k)*D)*r*sin(theta)+((m)*d+(k)*D)^2);  
    b_theta_sum=b_theta_sum+b_theta(kk,mm)'*b_theta(kk,mm);
    end
end

b=b(:);
b_theta=b_theta(:);
b_r=b_r(:);
phi_theta=(r*cos(theta)*(R^2+r^2-R*r*cos(theta))-2*R*r^2*sin(theta)^2)/(R^2+r^2-2*R*r*cos(theta))^(3/2);
phi_r=R*sin(theta)*(R-r*cos(theta))/(R^2+r^2-2*R*r*cos(theta))^(3/2);
a=[];a_theta=[];a_r=[];
nrr=0;
for nr=-(Nr-1)/2:(Nr-1)/2
    nrr=nrr+1;
    a(nrr)=sqrt(1/Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)));
    a_theta(nrr)=1i*2*pi*(nr)*d/lambda/sqrt(Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)))...
        *phi_theta;
    a_r(nrr)=1i*2*pi*(nr)*d/lambda/sqrt(Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)))...
    *phi_r;
end
a=a.';
a_theta=a_theta.';
a_r=a_r.';
h=kron(conj(b),a);
h_theta=kron(conj(b_theta),a)+kron(conj(b),a_theta);
h_r=kron(conj(b_r),a)+kron(conj(b),a_r);

h_theta2=h_theta'*h_theta;
h_r2=h_r'*h_r;

Q11=h_theta2-abs(h_theta'*h)^2
Q12=real(h_theta'*h_r)-real(h'*h_theta*h_r'*h)
Q22=h_r2-abs(h_r'*h)^2

DQ=Q11*Q22-Q12^2
CRB_r(rr)=sqrt(Q11/DQ/(K*M*Nr)/2)
CRB_theta(rr)=sqrt(Q22/DQ/(K*M*Nr)/2)

end
end



