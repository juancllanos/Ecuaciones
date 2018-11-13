% Metodo de Euler

%(1) y' = x + y
disp("y' = x + y");

% Solución análitica:
g = @(x) -x-1+(2/exp(-x));

% h = 0.1
xi =0;
xf =1;
h = 0.1;
f = @(x,y) x+y;
y0=1;

[W,X]=euler(xi,xf,h,f,y0);
[WR,XR] = Runge_Kutta(xi,xf,h,f,y0);
[Y,T] = sol(xi,xf,h,g);

disp("h = 0.1");
q = maxE(Y,W);
disp(q);
p = maxE(Y,WR);
disp(p);

%plot(X,W,XR,WR,T,Y);

% h = 0.05

h1 = 0.05;
[W1,X1]=euler(xi,xf,h1,f,y0);
[WR1,XR1] = Runge_Kutta(xi,xf,h1,f,y0);
[Y1,T1] = sol(xi,xf,h1,g);

disp("h = 0.05");
q1 = maxE(Y1,W1);
disp(q1);
p1 = maxE(Y1,WR1);
disp(p1);

%plot(X1,W1,XR1,WR1,T1,Y1);

% h = 0.02

h2 = 0.02;
[W2,X2]=euler(xi,xf,h2,f,y0);
[WR2,XR2] = Runge_Kutta(xi,xf,h2,f,y0);
[Y2,T2] = sol(xi,xf,h2,g);

disp("h = 0.02");
q2 = maxE(Y2,W2);
disp(q2);
p2 = maxE(Y2,WR2);
disp(p2);

%plot(X2,W2,XR2,WR2,T2,Y2);

% h = 0.01

h3 = 0.01;
[W3,X3]=euler(xi,xf,h3,f,y0);
[WR3,XR3] = Runge_Kutta(xi,xf,h3,f,y0);
[Y3,T3] = sol(xi,xf,h3,g);

disp("h = 0.01");
q3 = maxE(Y3,W3);
disp(q3);
p3 = maxE(Y3,WR3);
disp(p3);

%plot(X3,W3,XR3,WR3,T3,Y3);


% h = 0.005
h4 = 0.005;
[W4,X4]=euler(xi,xf,h4,f,y0);
[WR4,XR4] = Runge_Kutta(xi,xf,h4,f,y0);
[Y4,T4] = sol(xi,xf,h4,g);

disp("h = 0.005");
q4 = maxE(Y4,W4);
disp(q4);
p4 = maxE(Y4,WR4);
disp(p4);

%plot(X4,W4,XR4,WR4,T4,Y4);

% h = 0.001
h5 = 0.001;
[W5,X5]=euler(xi,xf,h5,f,y0);
[WR5,XR5] = Runge_Kutta(xi,xf,h5,f,y0);
[Y5,T5] = sol(xi,xf,h5,g);

disp("h = 0.001");
q5 = maxE(Y5,W5);
disp(q5);
p5 = maxE(Y5,WR5);
disp(p5);

%plot(X5,W5,XR5,WR5,T5,Y5);

%plot(X,W,X1,W1,X2,W2,X3,W3,X4,W4,X5,W5);

disp("---------------------------------");

%(2) y' + 3*x^2 * y = 6*x^2
disp("y' + 3*x^2 * y = 6*x^2");

% solución análitica:
gp = @(x) 2+exp(-x^3);

% h = 0.1
xip =0;
xfp =1;
hp = 0.1;
fp = @(x,y) 3*x^2*(2-y);
y0p =3;

[Wp,Xp]=euler(xip,xfp,hp,fp,y0p);
[WRp,XRp] = Runge_Kutta(xip,xfp,hp,fp,y0p);
[Yp,Tp] = sol(xip,xfp,hp,gp);

disp("h = 0.1");
qp = maxE(Yp,Wp);
disp(qp);
pp = maxE(Yp,WRp);
disp(pp);

%plot(Xp,Wp,XRp,WRp,Tp,Yp);

% h = 0.05

h1p = 0.05;
[W1p,X1p]=euler(xip,xfp,h1p,fp,y0p);
[WR1p,XR1p] = Runge_Kutta(xip,xfp,h1p,fp,y0p);
[Y1p,T1p] = sol(xip,xfp,h1p,gp);

disp("h = 0.05");
q1p = maxE(Y1p,W1p);
disp(q1p);
p1p = maxE(Y1p,WR1p);
disp(p1p);

%plot(X1p,W1p,XR1p,WR1p,T1p,Y1p);

% h = 0.02

h2p = 0.02;
[W2p,X2p]=euler(xip,xfp,h2p,fp,y0p);
[WR2p,XR2p] = Runge_Kutta(xip,xfp,h2p,fp,y0p);
[Y2p,T2p] = sol(xip,xfp,h2p,gp);

disp("h = 0.02");
q2p = maxE(Y2p,W2p);
disp(q2p);
p2p = maxE(Y2p,WR2p);
disp(p2p);

%plot(X2p,W2p,XR2p,WR2p,T2p,Y2p);

% h = 0.01

h3p = 0.01;
[W3p,X3p]=euler(xip,xfp,h3p,fp,y0p);
[WR3p,XR3p] = Runge_Kutta(xip,xfp,h3p,fp,y0p);
[Y3p,T3p] = sol(xip,xfp,h3p,gp);

disp("h = 0.01");
q3p = maxE(Y3p,W3p);
disp(q3p);
p3p = maxE(Y3p,WR3p);
disp(p3p);

%plot(X3p,W3p,XR3p,WR3p,T3p,Y3p);


% h = 0.005
h4p = 0.005;
[W4p,X4p]=euler(xip,xfp,h4p,fp,y0p);
[WR4p,XR4p] = Runge_Kutta(xip,xfp,h4p,fp,y0p);
[Y4p,T4p] = sol(xip,xfp,h4p,gp);

disp("h = 0.005");
q4p = maxE(Y4p,W4p);
disp(q4p);
p4p = maxE(Y4p,WR4p);
disp(p4p);

%plot(X4p,W4p,XR4p,WR4p,T4p,Y4p);

% h = 0.001
h5p = 0.001;
[W5p,X5p]=euler(xip,xfp,h5p,fp,y0p);
[WR5p,XR5p] = Runge_Kutta(xip,xfp,h5p,fp,y0p);
[Y5p,T5p] = sol(xip,xfp,h5p,gp);

disp("h = 0.001");
q5p = maxE(Y5p,W5p);
disp(q5p);
p5p = maxE(Y5p,WR5p);
disp(p5p);

%plot(X5p,W5p,XR5p,WR5p,T5p,Y5p);

%plot(Xp,Wp,X1p,W1p,X2p,W2p,X3p,W3p,X4p,W4p,X5p,W5p);

disp("---------------------------------");



%Funciones (Euler, Runge-Kutta y Error maximo):

function [W,X] = euler(xi,xf,h,f,y0)
X = xi:h:xf;
W(1) = y0;
for i=1:1:length(X)-1
    W(i+1)= W(i)+h*f(X(i),W(i));
end
end

function [W,t] = Runge_Kutta(xi,xf,h,f,y0)
t = xi:h:xf;
W(1)= y0;
for i = 1:1:length(t)-1
    k1 = f(t(i),W(i));
    k2 = f(t(i)+(1/2)*h,W(i)+(1/2)*h*k1);
    k3 = f(t(i)+(1/2)*h,W(i)+(1/2)*h*k2);
    k4 = f(t(i)+h,W(i)+h*k3);
    
    W(i+1) = W(i)+h*(1/6)*(k1+2*k2+2*k3+k4);
end
end

function [W,X] = sol(xi,xf,h,g)
X = xi:h:xf;
for i = 1:1:length(X)
    W(i) = g(X(i));
end
end

function x = maxE(X,Y)
for i = 1:1:length(X)
    E(i) = abs(X(i)-Y(i));
end
x = max(E);
end
    
    
    

