% Use MATLAB to plot numerical solutions based on data obtained by C/C++ code 
clear,clc
% data = load('RM1DP2T6.dat');
% data = load('Isentropic_RK3_LF_100000.dat');
data_ref = load('1DIsentropic_priU_RK3_IDEOS_EC6_200_T08.dat');
data_num = load('1DIsentropic_priU_RK3_IDEOS_ES5_200_T08.dat');
data_EC = load('1DIsentropic_priU_RK3_IDEOS_nonES5_200_T08.dat');

%% parameters
xmin = -0.4; xmax = 2;
T = 0.4;
Nx = 200;
Nx_ref = 200;
Nx_num = 200;
gpn = 3;
dx = (xmax - xmin)/Nx;
dx_ref = (xmax - xmin)/Nx_ref;
dx_num = (xmax - xmin)/Nx_num;
xe = (xmin-(gpn-0.5)*dx):dx:(xmax+(gpn-0.4)*dx); % 共N+6个点，且除去首尾6个虚拟单元的点，其他均为单元中点，即有限体积法的单元平均值
xe_ref = (xmin-(gpn-0.5)*dx_ref):dx_ref:(xmax+(gpn-0.4)*dx_ref);
xe_num = (xmin-(gpn-0.5)*dx_num):dx_num:(xmax+(gpn-0.4)*dx_num);

%% numerical solutions through C/C++ data
% r = data(:,1); 
% u = data(:,2);
% p = data(:,3);
r_num = data_num(:,1);
u_num = data_num(:,2);
p_num = data_num(:,3);
r_EC = data_EC(:,1);
u_EC = data_EC(:,2);
p_EC = data_EC(:,3);
r_ref = data_ref(:,1);
u_ref = data_ref(:,2);
p_ref = data_ref(:,3);

gamma = 4/3;

% r = r(1+gpn+1:end-gpn-1);
% u = u(1+gpn+1:end-gpn-1);
% p = p(1+gpn+1:end-gpn-1);

%% Plots
figure(1)
plot(xe_num(1+gpn:end-gpn),r_ref,'k-','linewidth',2)
hold on
plot(xe_num(1+gpn:end-gpn),r_num,'b:','linewidth',2)
% plot(xe(1+gpn:end-gpn),r,'b.','markersize',1)
plot(xe_ref(1+gpn:end-gpn),r_EC,'ro','linewidth',0.5,'markersize',3)
% axis square
legend('EC6 solution','ES5 solution','non-ES5 solution','location','northwest');
xlim(gca,[-0.4 1.1]);
% ylim(gca,[0.95 1.4]);
set(gca,'TickLength',[0.01 0.5],'XMinorTick','on','XTick',...
    [-0.4 0 0.4 0.8],'XTickLabel',{'-0.4','0','0.4','0.8'},'YMinorTick','on',...
    'YTick',[0.8 1 1.2 1.4]);
% print -depsc2 RHD1D_Isentropic_CMP_EOS1_density.eps
return
figure(2)
plot(xe(1+gpn:end-gpn),u,'bo','markersize',4)
% hold on
% plot(xe_num(1+gpn:end-gpn),u_num,'k-','linewidth',2)
% axis square
% legend('ES5 solution','reference solution','location','northeast','fontsize',6)
% xlim([0.5 0.53])
% print -depsc2 RHD1D_density_CMP_EOS1_velocity.eps
figure(3)
plot(xe(1+gpn:end-gpn),p,'bo','markersize',4)
% hold on
% plot(xe_num(1+gpn:end-gpn),p_num,'k-','linewidth',2)
% axis square
% legend('ES5 solution','reference solution','location','northeast','fontsize',6)
% xlim([0.5 0.53])
% print -depsc2 RHD1D_density_CMP_EOS1_pressure.eps