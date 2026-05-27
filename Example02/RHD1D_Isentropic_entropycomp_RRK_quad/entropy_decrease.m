% 画熵-DS随时间变化的曲线
clear,clc
data1 = load('tolEntr_1DIsentropic_200_RRK_EC6_T08_SO_dt.txt');
t1 = data1(:,1); S1 = data1(:,3);
data2 = load('tolEntr_1DIsentropic_200_RRK_ES5_T08_SO.txt');
t2 = data2(:,1); S2 = data2(:,3);
data3 = load('tolEntr_1DIsentropic_400_RRK_ES5_T08_SO.txt');
t3 = data3(:,1); S3 = data3(:,3);
data4 = load('tolEntr_1DIsentropic_800_RRK_ES5_T08_SO.txt');
t4 = data4(:,1); S4 = data4(:,3);
% data5 = load('tolEntr_1DSmooth800_ES5_T16.txt');
% t5 = data5(:,1); S5 = data5(:,2);
plot(t1,S1,'k-','linewidth',2)
hold on
plot(t2,S2,'g-.','linewidth',2)
plot(t3,S3,'b:','linewidth',2)
plot(t4,S4,'r--','linewidth',2)
% plot(t5,S5,'y.-','linewidth',2)
ylim([-11.72684,-11.72673])
% yt=-3.5e-5:0.5e-5:0;
% set(gca,'ytick',yt);
% legend("ES5: $400$ grids","ES5: $200$ grids",'location','southwest','interpreter','latex')
legend("EC6: $200$ grids","ES5: $200$ grids","ES5: $400$ grids","ES5: $800$ grids",'location','southwest','interpreter','latex')
% legend('$400$ grids','$200$ grids','$100$ grids','interpreter','latex')
print -depsc2 RHD_tolEntr_1DIsentr_EC6ES5_IDEOS_RRK_T08.eps