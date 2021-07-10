function [x1next,x2next,x3next,x4next] = Hau(x1,x2,x3,x4,h,f,g,m,n)

b1=1/8; b2=3/8;b3=3/8;b4=1/8;
a21=1/3;a31=-1/3; a32 = 1;a41 = 1; a42 = -1; a43 = 1;

k11 = f(x2);
k21 = g(x1,x3,x4);
k31 = m(x4);
k41 = n(x1,x2,x3);

k12 = f(x2+h*a21*k21);
k22 = g(x1+h*a21*k11,x3+h*a21*k31,x4+h*a21*k41);
k32 = m(x4+h*a21*k41);
k42 = n(x1+h*a21*k11,x2+h*a21*k21,x3+h*a21*k31);

k13 = f(x2+h*a31*k21+h*a32*k22);
k23 = g(x1+h*a31*k11+h*a32*k12,x3+h*a31*k31+h*a32*k32,x4+h*a31*k41+h*a32*k42);
k33 = m(x4+h*a31*k41+h*a32*k42);
k43 = n(x1+h*a31*k11+h*a32*k12,x2+h*a31*k21+h*a32*k22,x3+h*a31*k31+h*a32*k32);

k14 = f(x2+h*a41*k21+h*a42*k22+h*a43*k23);
k24 = g(x1+h*a41*k11+h*a42*k12+h*a43*k13,x3+h*a41*k31+h*a42*k32+h*a43*k33,x4+h*a41*k41+h*a42*k42+h*a43*k43);
k34 = m(x4+h*a41*k41+h*a42*k42+h*a43*k43);
k44 = n(x1+h*a41*k11+h*a42*k12+h*a43*k13,x2+h*a41*k21+h*a42*k22+h*a43*k23,x3+h*a41*k31+h*a42*k32+h*a43*k33);

x1next=x1 + h*(b1*k11 + b2*k12 + b3*k13 + b4*k14);
x2next=x2 + h*(b1*k21 + b2*k22 + b3*k23 + b4*k24);
x3next=x3 + h*(b1*k31 + b2*k32 + b3*k33 + b4*k34);
x4next=x4 + h*(b1*k41 + b2*k42 + b3*k43 + b4*k44);