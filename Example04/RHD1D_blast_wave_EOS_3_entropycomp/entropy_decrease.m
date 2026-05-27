% 画熵-DS随时间变化的曲线
clear,clc
data1 = load('tolEntr_1Dblast4000_ES5_T43.txt');
t1 = data1(:,1); S1 = data1(:,2);
data2 = load('tolEntr_1Dblast2000_ES5_T43.txt');
t2 = data2(:,1); S2 = data2(:,2);
data3 = load('tolEntr_1Dblast1000_ES5_T43.txt');
t3 = data3(:,1); S3 = data3(:,2);
% data4 = load('tolEntr_1DSmooth320_EOS_3_ES5.txt');
% t4 = data4(:,1); S4 = data4(:,2);
plot(t1,S1,'k-','linewidth',2)
hold on
plot(t2,S2,'r-.','linewidth',2)
plot(t3,S3,'b:','linewidth',2)
% plot(t4,S4,'r--','linewidth',3)
% ylim([-3.6e-5,0.3e-5])
% yt=-3.5e-5:0.5e-5:0;
% set(gca,'ytick',yt);
% legend("EC6: $320$ grids","ES5: $320$ grids",'location','southwest','interpreter','latex')
% legend("EC6: $80$ grids","ES5: $80$ grids","ES5: $160$ grids","ES5: $320$ grids",'location','southwest','interpreter','latex')
legend('$4000$ grids','$2000$ grids','$1000$ grids','interpreter','latex')
print -depsc2 RHD_tolEntr_1Dblast_ES5_4000_EOS_3.eps