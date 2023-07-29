clear all;
close all;
 
M=128;
Nr=1;
theta=45/180*pi

lambda=3e8/100e9;
d=lambda/2;
len=25
for i=1:2
    for k=1:2
    for rr=1:len
    r=rr*2
    R=80
    I=(i-1)*9+3;
    K=3+(k-1)*9
    D0=1/2*2^I*lambda;
    D=d*(M-1)+D0;
    D_array=D*(K-1)+d*(M-1); 
    [CRB_r_SW(i,k,rr),CRB_theta_SW(i,k,rr)]=WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)    
    [CRB_r_HSPW(i,k,rr),CRB_theta_HSPW(i,k,rr)]=WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d)
        [CRB_theta_PW(i,k,rr)]=WSMS_PW(theta,lambda,r,R,K,M,Nr,D,d)
    end
    end
end
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
semilogy(2:2:50,squeeze(CRB_r_SW(1,1,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_HSPW(1,1,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_SW(2,1,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_HSPW(2,1,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co2,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_SW(1,2,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_HSPW(1,2,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.7)
hold on
semilogy(2:2:50,squeeze(CRB_r_SW(2,2,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_r_HSPW(2,2,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.7)
hold on
% semilogy(2:2:50,squeeze(CRB_r_SW(3,1,:)),'dk-.', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(2:2:50,squeeze(CRB_r_HSPW(3,1,:)),'sk-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.7)
% hold on
% semilogy(2:2:50,squeeze(CRB_r_SW(3,2,:)),'sk-.', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(2:2:50,squeeze(CRB_r_HSPW(3,2,:)),'sk-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.7)
% hold on
grid on
axis([2,50,1e-5,50])
 
lgh=legend('SW-WSMS, $$I=3, K=3$$', 'HSPW-WSMS, $$I=3, K=3$$', 'SW-WSMS, $$I=12, K=3$$', 'HSPW-WSMS, $$I=12, K=3$$',...
    'SW-WSMS, $$I=3, K=12$$', 'HSPW-WSMS, $$I=3, K=12$$', 'SW-WSMS, $$I=12, K=12$$', 'HSPW-WSMS, $$I=12, K=12$$');
set(lgh,'interpreter','latex');
xlabel('Range $$r$$ (meters)','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{r}$$','interpreter','latex','fontsize',12)

figure
 
semilogy(2:2:50,squeeze(CRB_theta_SW(1,1,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_HSPW(1,1,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_PW(1,1,:)),'sk:', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on
 semilogy(2:2:50,squeeze(CRB_theta_SW(2,1,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_HSPW(2,1,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co2,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_PW(2,1,:)),'sk:', 'linewidth', 1, 'markerfacecolor',co2,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_SW(1,2,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_HSPW(1,2,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_PW(1,2,:)),'sk:', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_SW(2,2,:)),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_HSPW(2,2,:)),'ok-.', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on
semilogy(2:2:50,squeeze(CRB_theta_PW(2,2,:)),'sk:', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on
grid on
axis([2,50,2e-7,5e-4])

lgh=legend('SW-WSMS, $$I=3, K=3$$', 'HSPW-WSMS, $$I=3, K=3$$','PW-WSMS, $$I=3, K=3$$', 'SW-WSMS, $$I=12, K=3$$', 'HSPW-WSMS, $$I=12, K=3$$',...
    'PW-WSMS, $$I=12, K=3$$','SW-WSMS, $$I=3, K=12$$', 'HSPW-WSMS, $$I=3, K=12$$','PW-WSMS, $$I=3, K=12$$', 'SW-WSMS, $$I=12, K=12$$', 'HSPW-WSMS, $$I=12, K=12$$',...
    'PW-WSMS, $$I=12, K=12$$');
set(lgh,'interpreter','latex');
 
xlabel('Range $$r$$ (meters)','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{\theta}$$','interpreter','latex','fontsize',12)