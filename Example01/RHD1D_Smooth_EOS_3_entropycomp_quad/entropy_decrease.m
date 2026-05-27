% 画熵-DS随时间变化的曲线
clear,clc
data1 = load('tolEntr_1DSmooth80_EOS_3_EC4_quad.txt');
t1 = data1(:,1); S1 = data1(:,3);
data2 = load('tolEntr_1DSmooth80_EOS_3_ES4_quad.txt');
t2 = data2(:,1); S2 = data2(:,3);
data3 = load('tolEntr_1DSmooth160_EOS_3_ES4_quad.txt');
t3 = data3(:,1); S3 = data3(:,3);
data4 = load('tolEntr_1DSmooth320_EOS_3_ES4_quad.txt');
t4 = data4(:,1); S4 = data4(:,3);
plot(t1,S1,'k-','linewidth',2)
hold on
plot(t2,S2,'g-.','linewidth',2)
plot(t3,S3,'b:','linewidth',2)
plot(t4,S4,'r--','linewidth',3)
% ax = gca();
% ax.YRuler.Exponent = -7;
% set(gca,'YTick',...
%     [-1.6e-07 -1.4e-07 -1.2e-07 -1e-07 -8e-08 -6e-08 -4e-08 -2e-08 0 2e-08]);
% ylim([-1.6e-7,0.2e-7]);
% ylim([-11.27542650,-11.27542630]); % EC6ES5
ylim([-11.2754272,-11.2754262]); % EC4ES4
% legend("EC6: $320$ grids","ES5: $320$ grids",'location','southwest','interpreter','latex')
legend("SSP-RK3-EC4: $80$ grids","SSP-RK3-ES4: $80$ grids","SSP-RK3-ES4: $160$ grids","SSP-RK3-ES4: $320$ grids",'location','southwest','interpreter','latex')
% legend('$400$ grids','$200$ grids','$100$ grids','interpreter','latex')
print -depsc2 RHD_tolEntr_1DSmooth_EC4ES4_EOS_3.eps