clear all;
close all;

M=128; 
theta=0/180*pi

lambda=3e8/100e9;
d=lambda/2;
len=30
for nr=1:3
    for rr=1:len
        K=12;
    r=rr*1
    R=31
    Nr=1+(nr-1)*17
    I=10;
    D0=1/2*2^I*lambda;
    D=d*(M-1)+D0;
    D_array=D*(K-1)+d*(M-1); 
    [CRB_r_SW(nr,rr),CRB_theta_SW(nr,rr)]=WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)    
    [CRB_r_HSPW(nr,rr),CRB_theta_HSPW(nr,rr)]=WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d)
        [CRB_theta_PW(nr,rr)]=WSMS_PW(theta,lambda,r,R,K,M,Nr,D,d)
    end
end

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
% figure
% semilogy(1:30,CRB_r_SW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_HSPW(1,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_SW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_HSPW(2,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_HSPW(3,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
% hold on
% % semilogy(-1.5:0.125:1.5,CRB_r_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% % hold on
% % semilogy(-1.5:0.125:1.5,CRB_r_HSPW(3,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.7)
% % hold on
% grid on
% figure
 
semilogy(1:30,CRB_theta_SW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW(1,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_PW(1,:),'sk:', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on
 semilogy(1:30,CRB_theta_SW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW(2,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_PW(2,:),'sk:', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on
 semilogy(1:30,CRB_theta_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_HSPW(3,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_PW(3,:),'sk:', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on
% semilogy(-1.5:0.125:1.5,CRB_theta_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(-1.5:0.125:1.5,CRB_theta_HSPW(3,:),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
% hold on
% semilogy(-1.5:0.125:1.5,CRB_theta_PW(3,:),'sk:', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
% hold on
grid on
axis([1,30,2e-9,6e-5])

lgh=legend('SW-WSMS, $$N_r=1$$', 'HSPW-WSMS, $$N_r=1$$','PW-WSMS, $$N_r=1$$', 'SW-WSMS, $$N_r=18$$', 'HSPW-WSMS, $$N_r=18$$',...
    'PW-WSMS, $$N_r=18$$','SW-WSMS, $$N_r=35$$', 'HSPW-WSMS, $$N_r=35$$','PW-WSMS, $$N_r=35$$');
set(lgh,'interpreter','latex');
 
xlabel('Range $$r$$ (meters)','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{\theta}$$','interpreter','latex','fontsize',12)