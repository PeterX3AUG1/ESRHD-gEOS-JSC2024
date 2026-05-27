% 画熵-DS随时间变化的曲线
clear,clc
data1 = load('tolEntr_2DSmooth80_EOS_2_EC6.txt');
t1 = data1(:,1); S1 = data1(:,2);
data2 = load('tolEntr_2DSmooth80_EOS_2_ES5.txt');
t2 = data2(:,1); S2 = data2(:,2);
data3 = load('tolEntr_2DSmooth160_EOS_2_ES5.txt');
t3 = data3(:,1); S3 = data3(:,2);
data4 = load('tolEntr_2DSmooth320_EOS_2_ES5.txt');
t4 = data4(:,1); S4 = data4(:,2);

plot(t1,S1,'k-','linewidth',2)
hold on
plot(t2,S2,'g-.','linewidth',2)
plot(t3,S3,'b:','linewidth',2)
plot(t4,S4,'r--','linewidth',2)

legend("EC6: $80\times80$ grids","ES5: $80\times80$ grids","ES5: $160\times160$ grids","ES5: $320\times320$ grids",'location','southwest','interpreter','latex')
% legend('$400$ grids','$200$ grids','$100$ grids','interpreter','latex')
print -depsc2 RHD_tolEntr_2DSmooth_EOS_2.eps