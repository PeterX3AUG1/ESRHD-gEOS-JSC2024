% 画熵-DS随时间变化的曲线
clear,clc
data1 = load('tolEntr_1DRP1_400_RRK_ES5_SO.txt');
t1 = data1(:,1); S1 = data1(:,2);
data2 = load('tolEntr_1DRP1_200_RRK_ES5_SO.txt');
t2 = data2(:,1); S2 = data2(:,2);
data3 = load('tolEntr_1DRP1_100_RRK_ES5_SO.txt');
t3 = data3(:,1); S3 = data3(:,2);
% data4 = load('tolEntr_1DSmooth320_EOS_3_ES5.txt');
% t4 = data4(:,1); S4 = data4(:,2);
plot(t1,S1,'k-','linewidth',2)
hold on
plot(t2,S2,'r-.','linewidth',2)
plot(t3,S3,'b:','linewidth',2)
% plot(t4,S4,'r--','linewidth',3)
xlim([0,0.4])
% legend("EC6: $320$ grids","ES5: $320$ grids",'location','southwest','interpreter','latex')
% legend("EC6: $80$ grids","ES5: $80$ grids","ES5: $160$ grids","ES5: $320$ grids",'location','southwest','interpreter','latex')
legend('$400$ grids','$200$ grids','$100$ grids','location','northeast','interpreter','latex')
print -depsc2 RHD_tolEntr_1DRP1_ES5_RRK_EOS_1.eps