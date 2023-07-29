clear all;
close all;

M=128;
Nr=1;
theta=0/180*pi

lambda=3e8/100e9;
d=lambda/2;
len=21

    for i=1:len
        K=2;
    r=10
    R=50
    Nr=12
   % psi0=pi/(len-1)*(p-1);
    I=i-1
    D0=2^I*lambda;
    D=d*(M-1)+D0;
%     D_array=D*(K-1)+d*(M-1);  
    [~,CRB_theta_SW(i)]=WSMS_HSPW(theta,lambda,r,R,K,M,Nr,D,d) 
    Upper(i)= sqrt(lambda^2/8/pi^2/(d^2*(M^2-1)/12)/K/M/Nr)
    Lower(i)= sqrt(lambda^2/8/pi^2/(r^2+d^2*(M^2-1)/12)/K/M/Nr)
    end

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
semilogy(0:1:20,CRB_theta_SW,'^k-', 'linewidth', 1, 'markerfacecolor', co1,'markersize', 6.5)
hold on
semilogy(0:1:20,Upper,'>k:', 'linewidth', 1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(0:1:20,Lower,'<k:', 'linewidth', 1, 'markerfacecolor', co3,'markersize', 6.5)
grid on
axis([1,20,2e-7,3e-4])

lgh=legend('HSPW-WSMS, $$\theta=0$$','Asymptotic case of $$\psi_0\rightarrow 0$$','Asymptotic case of $$\psi_0\rightarrow \pi$$');
set(lgh,'interpreter','latex','fontsize',12);
xlabel('Inter-Subarray Spacing Setting $$I$$','interpreter','latex','fontsize',14)
ylabel('Root $$\textbf{CRB}_{\theta}$$','interpreter','latex','fontsize',14)