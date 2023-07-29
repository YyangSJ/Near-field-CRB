clear all;
close all;

M=128;
Nr=1;
theta=45/180*pi
    R=40
lambda=3e8/100e9;
d=lambda/2;
len=30
for k=1:4
    for rr=1:len
        K=3+(k-1)*3;
  
 r=rr*1
    I=6;
    D0=1/2*2^I*lambda;
    D=d*(M-1)+D0;
    D_array=D*(K-1)+d*(M-1); 
    [CRB_r_SW(k,rr),CRB_theta_SW(k,rr)]=WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)    
    [CRB_r_SW2(k,rr),CRB_theta_SW2(k,rr)]=WSMS_SW2(theta,lambda,r,R,K,M,Nr,D,d) 
%     [CRB_r_HSPW(k,rr),CRB_theta_HSPW(k,rr)]=WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d)    
%     [CRB_r_HSPW2(k,rr),CRB_theta_HSPW2(k,rr)]=WSMS_HSPW2(theta,lambda,r,R,K,M,Nr,D,d) 
    end
end
co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255

figure
semilogy(1:30,CRB_r_SW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_SW2(1,:),'ok--', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_r_SW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_SW2(2,:),'ok--', 'linewidth', 1, 'markerfacecolor',co2,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_r_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_SW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_r_SW(4,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_r_SW2(4,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on 
% semilogy(1:30,CRB_r_HSPW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
% hold on
% semilogy(1:30,CRB_r_HSPW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
% hold on 
grid on
axis([1,30,2.5e-5,3])
lgh=legend('SW-WSMS, $$K=3$$', 'SW-WSMS Approx., $$K=3$$','SW-WSMS, $$K=6$$', 'SW-WSMS Approx., $$K=6$$'...
    , 'SW-WSMS, $$K=9$$','SW-WSMS Approx., $$K=9$$','SW-WSMS, $$K=12$$','SW-WSMS Approx., $$K=12$$');
set(lgh,'interpreter','latex');
 
xlabel('Range $$r$$ (meters)','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{r}$$','interpreter','latex','fontsize',12)

figure
semilogy(1:30,CRB_theta_SW(1,:),'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_SW2(1,:),'ok--', 'linewidth', 1, 'markerfacecolor',co1,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_theta_SW(2,:),'^k-', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_SW2(2,:),'ok--', 'linewidth', 1, 'markerfacecolor',co2,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_theta_SW(3,:),'^k-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_SW2(3,:),'ok--', 'linewidth', 1, 'markerfacecolor',co4,'markersize', 6.5)
hold on 
semilogy(1:30,CRB_theta_SW(4,:),'^k-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
hold on
semilogy(1:30,CRB_theta_SW2(4,:),'ok--', 'linewidth', 1, 'markerfacecolor',co5,'markersize', 6.5)
hold on 
grid on
axis([1,30,9.5e-6,1.5e-4])
lgh=legend('SW-WSMS, $$K=3$$', 'SW-WSMS Approx., $$K=3$$','SW-WSMS, $$K=6$$', 'SW-WSMS Approx., $$K=6$$'...
    , 'SW-WSMS, $$K=9$$','SW-WSMS Approx., $$K=9$$','SW-WSMS, $$K=12$$','SW-WSMS Approx., $$K=12$$');
set(lgh,'interpreter','latex');
 
xlabel('Range $$r$$ (meters)','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{\theta}$$','interpreter','latex','fontsize',12)