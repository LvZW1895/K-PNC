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
data1=textread('data1.txt');
data2=textread('data2.txt');
data3=textread('clures1.txt');
data4=textread('clures2.txt');
[m,n]=size(data3)
%C = linspecer(m)
X1=data1(:,1);
Y1=data1(:,2);
X2=data2(:,1);
Y2=data2(:,2);
Z1=data3(:,1);
Z2=data4(:,1);
figure(1);
hold on;
c1=3*Z1
c2=3*Z2
scatter(X1,Y1,'filled','cdata',c1)  
figure(2);
hold on;
scatter(X2,Y2,'filled','cdata',c2)  




