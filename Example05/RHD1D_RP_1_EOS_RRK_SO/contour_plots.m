% Use MATLAB to plot numerical solutions based on data obtained by C/C++ code 
clear,clc
% data = load('RM1DP2T6.dat');
data = load('RP_1_priU_RK3_EOS_1_ES5_400.dat');
data_num = load('RP_1_priU_RK3_EOS_1_LF_100000.dat');

%% parameters
xmin = 0; xmax = 1;
T = 0.4;
Nx = 400;
Nx_num = 100000;
gpn = 3;
dx = (xmax - xmin)/Nx;
dx_num = (xmax - xmin)/Nx_num;
xe = (xmin-(gpn-0.5)*dx):dx:(xmax+(gpn-0.4)*dx); % 共N+6个点，且除去首尾6个虚拟单元的点，其他均为单元中点，即有限体积法的单元平均值
xe_num = (xmin-(gpn-0.5)*dx_num):dx_num:(xmax+(gpn-0.4)*dx_num);

%% numerical solutions through C/C++ data
r = data(:,1); 
u = data(:,2);
p = data(:,3);
r_num = data_num(:,1);
u_num = data_num(:,2);
p_num = data_num(:,3);

gamma = 4/3;

% r = r(1+gpn+1:end-gpn-1);
% u = u(1+gpn+1:end-gpn-1);
% p = p(1+gpn+1:end-gpn-1);

%% Plots
figure(1)
plot(xe(1+gpn:end-gpn),r,'bo','markersize',4)
hold on
plot(xe_num(1+gpn:end-gpn),r_num,'k-','linewidth',2)
axis square
% legend('ES5 solution','reference solution','location','northeast','fontsize',6);
% xlim([0.5 0.53])
% print -depsc2 RHD1D_density_CMP_EOS1_density.eps
figure(2)
plot(xe(1+gpn:end-gpn),u,'bo','markersize',4)
hold on
plot(xe_num(1+gpn:end-gpn),u_num,'k-','linewidth',2)
axis square
% legend('ES5 solution','reference solution','location','northeast','fontsize',6)
% xlim([0.5 0.53])
% print -depsc2 RHD1D_density_CMP_EOS1_velocity.eps
figure(3)
plot(xe(1+gpn:end-gpn),p,'bo','markersize',4)
hold on
plot(xe_num(1+gpn:end-gpn),p_num,'k-','linewidth',2)
axis square
% legend('ES5 solution','reference solution','location','northeast','fontsize',6)
% xlim([0.5 0.53])
% print -depsc2 RHD1D_density_CMP_EOS1_pressure.eps