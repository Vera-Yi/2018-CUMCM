clear;clc;
format short e
a=1/(300*1377);
h=0.0006;
k=0.082;
M=89;
ot=0.1;
n=5400;
ox=h/M;
r=a*ot/(ox)^2;
x0=zeros(M+1,1);
 
for i=1:M
    x0(i+1)=i*ox;
end
p=37*ones(M+1,1);
q=zeros(M+1,5400);
T=[p,q];

for i=1:n   %数据的输入
    B=zeros(M-1,1);   %存放系数矩阵主对角线元素
    A=zeros(M-2,1);   %存放系数矩阵主对角线元素下方次对角线的元素
    C=zeros(M-2,1);   %存放系数矩阵主对角线元素下方次对角线的元素
    S=zeros(M-1,1);   %存放右端的常数项
    
    for j=1:M-2
        B(j)=1+2*r*k;A(j)=-r*k;C(j)=-r*k;
        S=T(2:M,i);
    end
   
    B(M-1)=1+2*r*k;
    S(M-1)=T(M,1);
    T(1,i)=75;
    T(M+1,i)=48;
    S(1,1)=S(1,1)+r*k*T(1,i);
    S(M-1,1)=S(M-1,1)+r*k*T(M+1,i);  %追赶法
    S(1)=S(1)/B(1);
    tt=B(1);
    k1=2;
 
    while k1~=M
        B(k1-1)=C(k1-1)/tt;
        tt=B(k1)-A(k1-1)*B(k1-1);
        S(k1)=(S(k1)-A(k1-1)*S(k1-1))/tt;
        k1=k1+1;
    end
    k1=1;
    while k1~=M-1
        S(M-1-k1)=S(M-1-k1)-B(M-1-k1)*S(M-k1);
        k1=k1+1;
    end
 
    T(2:M,i+1)=S;   %把结果放入矩阵T中
   
end

%  fid=fopen('test.txt','w');
%  fprintf(fid,'format long\n',T);
% fclose(fid);
