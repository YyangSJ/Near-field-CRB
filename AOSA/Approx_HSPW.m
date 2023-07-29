clear all;
close all;

M=128;
Nr=1;
theta=30/180*pi
    R=80
lambda=3e8/100e9;
d=lambda/2;
len=30
for k=1:3
    for rr=1:len
        K=3+6*(k-1);
  
 r=rr*1
    I=6;
    D0=1/2*2^I*lambda;
    D=d*(M-1)+D0;
    D_array=D*(K-1)+d*(M-1); 
%     [CRB_r_SW(k,rr),CRB_theta_SW(k,rr)]=WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)    
%     [CRB_r_SW2(k,rr),CRB_theta_SW2(k,rr)]=WSMS_SW2(theta,lambda,r,R,K,M,Nr,D,d) 
    [CRB_r_HSPW(k,rr),CRB_theta_HSPW(k,rr)]=WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d)    
    [CRB_r_HSPW2(k,rr),CRB_theta_HSPW2(k,rr)]=WSMS_HSPW2(theta,lambda,r,R,K,M,Nr,D,d) 
    end
end
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255

figure
semilogy(1:30,CRB_r_HSPW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_HSPW2(1,:),'ok--', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_r_HSPW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_HSPW2(2,:),'ok--', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_r_HSPW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_HSPW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on 
% semilogy(1:30,CRB_r_HSPW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_HSPW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
% hold on 
grid on

semilogy(1:30,CRB_theta_HSPW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW2(1,:),'ok--', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_theta_HSPW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW2(2,:),'ok--', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_theta_HSPW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on 
grid on