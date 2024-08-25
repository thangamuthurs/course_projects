%For N=4
%%N- No. of elements, l- length of system, 
clc
clear all
close all
[Q,N,l]=deal(1,4,1);      
n=N+1;  
he=l/N;      %an element's size
syms x s q;
f= @(x) Q*x;         %No. of nodes
%given
[g,h]=deal(0,0);
%Initializing K & F to zero, N no. of equations
K=zeros(N,N);
F=zeros(N,1);
a=0:he:1;

%constructing K:
for i=1:N-1
    K(i,i) = K(i,i)+1;
    K(i,i+1)=K(i,i+1)-1;
    K(i+1,i)=K(i+1,i)-1;
    K(i+1,i+1)=K(i+1,i+1)+1;
    F(i,1)=F(i, 1) + (2*f(a(i)) + f(a(i+1)));
    F(i+1,1)=F(i+1,1) + f(a(i)) + 2*f(a(i+1));
end
K(N,N)=K(N,N)+1;
F(N,1)=F(N,1)+(2*f(a(N)) + f(a(N+1)));
F=sym((he/6)*F);
K=(1/he)*K;
d=(linsolve(K,F));
D=sym(d);

NA=sym(zeros(n,N));
X=(0:N)/N;

%creating shape functions:
for i=2:N
    for j=i
        NA(i,j-1)=(x-X(j-1))/he;
        NA(i,j)=(X(j+1)-x)/he;
    end
end
NA(n,N)=(x-X(N))/he;
NA(1,1)=(X(2)-x)/he;

s=sym(zeros(1,N));
for i=1:N   
    S(i)=s(i);
    for j=1:N
        S(i)=D(j)*NA(j,i)+S(i);
    end
end

R=sym(zeros(1,N));
figure(1);
%Exact Solution:
u= @(Q) 1/6-Q^3/6;      %Strong form
fplot(u,[0,1],'k','LineWidth',1,'LineStyle','-',Marker='*');
hold on ;
for i=1:N
    x1=a(i);
    x2=a(i+1);
    if i==1
        fplot(S(i),[x1,x2],'r','LineWidth',1.5,'DisplayName','Appprox');
        hold on
    else
        fplot(S(i),[x1,x2],'r','LineWidth',2);
        hold on
    end
end
hold off
xlabel('Nodes');
ylabel('Solution');
title('Approximate vs Exact Solution')
legend('u','approx')
hold off
figure(2);
ux=diff(sym(u),Q);      %Slope of strong form

q=linspace(0,1,N);
v=sym(zeros(1,N));
fplot(ux,[0,1],'k','LineWidth',2);
hold on;
R=sym(zeros(1,N));
for i=1:N;
    x1=a(i);
    x2=a(i+1);
    x3=0.5*(x1+x2);
    v(i)=x3;
    eq=S(i);
    eq_approx_slope=diff(eq,x);
    R1(i)=subs(eq,'x',x3);
    R2(i)=subs(u,'Q',x3);
    approx_slope(i)=subs(eq_approx_slope,'x',x3);
    exact_slope(i)=subs(ux,'Q',x3);
    fplot(approx_slope(i),[x1,x2],'r','LineWidth',2);
    hold on;
end
hold off;
Relative_error=abs(round(vpa(approx_slope-exact_slope)*2,8))
legend('ux','approx_slope')
xlabel('Nodes');
ylabel('Slope');
title('Approximate slope vs Exact Slope')
hold off