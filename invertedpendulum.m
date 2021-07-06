%% Mehmet KILIÇ - 498 TP
%% System Specs.
T=1e-3;
M=0.5; %kg
m=0.2; %kg
b=0.1; %N/m/sec
I=0.006; %kg.^2
g=9.8;
l=0.3; %m
%% From CTSS to DTSS Model
denum=(I*(m+m)+M*m*(l^2));
A22=-((I+m*(l^2))*b)/denum;
A23=(m^2*g*l^2)/denum;
A42=-m*l*b/denum;
A43=m*g*l*(M+m)/denum;
B21=(I+m*(l^2))/denum;
B41=m*l/denum;

Ac=[0 1 0 0;
    0 A22 A23 0;
    0 0 0 1;
    0 A42 A42 0];
Bc=[0;B21;0;B41];
Cc=[1 0 0 0;
    0 0 1 0];
Dc=[0;0];
state={'x' 'x_dot' 'theta' 'theta_dot'};
input={'u'};
output={'x';'theta'};
ctss=ss(Ac,Bc,Cc,Dc,'inputname',input,'statename',state,'outputname',output);
dtss=c2d(ctss,T,'zoh');

A=dtss.A;
B=dtss.B;
C=dtss.C;

%% Plant Properties
stability=eig(A);
controllability=rank(ctrb(A,B));
observability=rank(obsv(A,C));
% Plant is stable, controllable and observable!

[num,den]=ss2tf(Ac,Bc,Cc,Dc);
tf1=tf(num(1,:),den);
tf2=tf(num(2,:),den);

%% LQR Design for N=10

Q=C'*C;
Qf=Q;
R=0.001;
x0=[1;5;1;5];
x=[x0;zeros(36,1)];
u=zeros(9,1);
P=zeros(44,4);
P(41:44,1:4)=Qf;
K=zeros(10,4);
for i=11:-1:2
    Pp=P((4*(i+1)-7):(4*(i+1)-4),1:4);
    P((4*i-7):(4*i-4),1:4)=A'*Pp*A+Q-A'*Pp*B*((B'*Pp*B)^-1)*B'*Pp*A;
end
i=1;
for i=0:1:9
    Pp=P(4*i+5:4*i+8,:);
    K(i+1,:)=((R+B'*Pp*B)^-1)*B'*Pp*A;
    u(i+1,1)=-K(i+1,:)*x(4*i+1:4*i+4,:)
    x(4*i+5:4*i+8,:)=A*x(4*i+1:4*i+4,:)+B*u(i+1,1);
end