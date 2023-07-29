clear all;
close all;

M=128;
Nr=1;


lambda=3e8/100e9;
d=lambda/2;
len=13
for rr=1:1
    for tt=1:len
        
        K=3;
        r=rr*10
        theta=-1.5+3/(len-1)*(tt-1);
        theta=0;
        R=40
        I=tt;
        D0=1/2*2^I*lambda; 
        D=d*(M-1)+D0;
        D_array=D*(K-1)+d*(M-1);
        [CRB_r_SW(rr,tt),CRB_theta_SW(rr,tt)]=WSMS_SW(theta,lambda,r,R,K,M,Nr,D,d)
        [CRB_r_UA(rr,tt),CRB_theta_UA(rr,tt)]=UA(theta,lambda,r,R,K,M,Nr,D,d)
        [CRB_r_DUA(rr,tt),CRB_theta_DUA(rr,tt)]=DUA(theta,lambda,r,R,K,M,Nr,D,d)
    end
end

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure

semilogy(1:len,CRB_r_SW,'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:len,CRB_r_UA,'sk-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:len,CRB_r_DUA,'ok-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)
grid on
axis([1,len,7e-5,1])
lgh=legend('SW-WSMS','SW-UA, $$d^\prime=(D(K-1)+(M-1)d)/(KM-1)$$','SW-DUA, $$d^\prime=\lambda/2$$');
set(lgh,'interpreter','latex','fontsize',11)
xlabel('Inter-Subarray Spacing setting $$I$$','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{r}$$','interpreter','latex','fontsize',12)
figure
semilogy(1:len,CRB_theta_SW,'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(1:len,CRB_theta_UA,'sk-', 'linewidth', 1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(1:len,CRB_theta_DUA,'ok-', 'linewidth', 1, 'markerfacecolor', co5,'markersize', 6.5)

lgh=legend('SW-WSMS','SW-UA, $$d^\prime=(D(K-1)+(M-1)d)/(KM-1)$$','SW-DUA, $$d^\prime=\lambda/2$$');
set(lgh,'interpreter','latex','fontsize',11)

grid on
axis([1,len,2e-6,3e-4])
xlabel('Inter-Subarray Spacing setting $$I$$','interpreter','latex','fontsize',12)
ylabel('Root $$\textbf{CRB}_{\theta}$$','interpreter','latex','fontsize',12)
