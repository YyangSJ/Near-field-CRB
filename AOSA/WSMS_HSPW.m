function [CRB_r,CRB_theta] = WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d)

phi_theta=(r*cos(theta)*(R^2+r^2-R*r*cos(theta))-2*R*r^2*sin(theta)^2)/(R^2+r^2-2*R*r*cos(theta))^(3/2);
phi_r=R*sin(theta)*(R-r*cos(theta))/(R^2+r^2-2*R*r*cos(theta))^(3/2);

at=[]; at_theta=[];
mm=0;
for m=-(M-1)/2:(M-1)/2
    mm=mm+1;
    at(mm)=sqrt(1/M)*exp(1i*2*pi/lambda*(m)*d*sin(theta));
    at_theta(mm)=1i*2*pi*(m)*d*cos(theta)/lambda/sqrt(M)*exp(1i*2*pi/lambda*(m)*d*sin(theta));
end
w=[];w_theta=[];w_r=[];
kk=0;
for k=-(K-1)/2:(K-1)/2
    kk=kk+1;
    w(kk)=sqrt(1/K)*exp(-1i*2*pi/lambda*sqrt(r^2-2*(k)*D*r*sin(theta)+(k)^2*D^2));
    w_theta(kk)=1i*2*pi/lambda*(k)*D*r*cos(theta)/sqrt(r^2-2*(k)*D*r*sin(theta)+(k)^2*D^2)*w(kk);
    w_r(kk)=1i*2*pi/lambda*((k)*D*sin(theta)-r)/sqrt(r^2-2*(k)*D*r*sin(theta)+(k)^2*D^2)*w(kk);
end
at=at.';
at_theta=at_theta.';
w=w.';
w_theta=w_theta.';
w_r=w_r.';

ar=[];ar_theta=[];ar_r=[];
nrr=0;
for nr=-(Nr-1)/2:(Nr-1)/2
    nrr=nrr+1;
    ar(nrr)=sqrt(1/Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)));
    ar_theta(nrr)=1i*2*pi*(nr)*d/lambda/sqrt(Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)))...
        *phi_theta;
    ar_r(nrr)=1i*2*pi*(nr)*d/lambda/sqrt(Nr)*exp(1i*2*pi/lambda*((nr)*d*r*sin(theta))/sqrt(R^2+r^2-2*R*r*cos(theta)))...
        *phi_r;
end
ar=ar.';
ar_theta=ar_theta.';
ar_r=ar_r.';

h=kron(kron(conj(w),conj(at)),ar);
h_theta=kron(kron(conj(w_theta),conj(at)),ar)+kron(kron(conj(w),conj(at_theta)),ar)+kron(kron(conj(w),conj(at)),ar_theta);
h_r=kron(kron(conj(w_r),conj(at)),ar)+kron(kron(conj(w),conj(at)),ar_r);

h_theta2=h_theta'*h_theta;
h_r2=h_r'*h_r;

Q11=h_theta2-abs(h_theta'*h)^2;
Q12=real(h_theta'*h_r)-real(h'*h_theta*h_r'*h);
Q22=h_r2-abs(h_r'*h)^2;
DQ=Q11*Q22-Q12^2;

CRB_r=sqrt(Q11/DQ/(K*M*Nr)/2);
CRB_theta=sqrt(Q22/DQ/(K*M*Nr)/2);


end

