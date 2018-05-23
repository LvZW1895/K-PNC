clear;
% data1=textread('pnc_result.txt');
% data2=textread('comp_result.txt');
% X1=data1(:,1);
% Y1=data1(:,2);
% X2=data2(:,1);
% Y2=data2(:,2);
% semilogy(X1,Y1,'-ob');
% hold on;
% semilogy(X2,Y2,'-^r');
% axis([0,10,0,1]);
% grid on;
data1=textread('comp_result_ber.txt');
data2=textread('nidealcomp_result4_ber.txt');
data3=textread('nidealcomp_result8_ber.txt');
data4=textread('pnc_result_ber.txt');
X1=data1(:,1);
Y1=data1(:,2);
X2=data2(:,1);
Y2=data2(:,2);
X3=data3(:,1);
Y3=data3(:,2);
X4=data4(:,1);
Y4=data4(:,2);
semilogy(X1,Y1,'-ob');
hold on;
semilogy(X2,Y2,'-^r');
semilogy(X3,Y3,'-xg');
%semilogy(X4,Y4,'-s');
legend('idealCoMP','4bit nonidealCoMP','8bit nonidealCoMP','PNC')
xlabel('E_b/N_0(dB)')
grid on;
