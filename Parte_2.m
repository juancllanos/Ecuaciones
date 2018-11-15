% Parte 2, ejercicio 9.5 de la página 541.

%(a) a = 1, b = 1
a = 1;
b = 1;
h = 0.005;
x0 = 5;
y0 = 2;

[X,Y,t] = RungeKutta9(a,b,h,x0,y0);

% Grafica x(t) & y(t) vs t
%plot(t,X,t,Y);
%xlabel('t')
%ylabel('X(t)    |    Y(t)');

%Plano de fase 
%plot(X,Y);
%xlabel('X');
%ylabel('Y');

%(b) a = 3 y a = 1/3 con b = 1

% a = 3
a1 = 3;

[X1,Y1,t1] = RungeKutta9(a1,b,h,x0,y0);

% Grafica x(t) & y(t) vs t
%plot(t1,X1,t1,Y1);
%xlabel('t')
%ylabel('X(t)    |    Y(t)');

% Plano de fase 
%plot(X1,Y1);
%xlabel('X');
%ylabel('Y');

% a = 1/3
a2 = 1/3;

[X2,Y2,t2] = RungeKutta9(a2,b,h,x0,y0);

% Grafica x(t) & y(t) vs t
%plot(t2,X2,t2,Y2);
%xlabel('t')
%ylabel('X(t)    |    Y(t)');

% Plano de fase 
%plot(X2,Y2);
%xlabel('X');
%ylabel('Y');

%(c) b = 3 y b = 1/3 con a = 1

% b = 3
[X3,Y3,t3] = RungeKutta9(a,a1,h,x0,y0);

% Grafica x(t) & y(t) vs t
%plot(t3,X3,t3,Y3);
%xlabel('t')
%ylabel('X(t)    |    Y(t)');

% Plano de fase 
%plot(X3,Y3);
%xlabel('X');
%ylabel('Y');

% b = 1/3

[X4,Y4,t4] = RungeKutta9(a,a2,h,x0,y0);

% Grafica x(t) & y(t) vs t
%plot(t4,X4,t4,Y4);
%xlabel('t')
%ylabel('X(t)    |    Y(t)');

% Plano de fase 
plot(X,Y,X1,Y1,X2,Y2,X3,Y3,X4,Y4);
xlabel('X');
ylabel('Y');
legend()

%-----------------------------------------------------
% Función para este ejercicio :
function [X,Y,t] = RungeKutta9(a,b,h,x0,y0)
t = 0:h:30;

xp = @(x,y) a*x*(1-(y/2));
yp = @(x,y) b*y*(-1+(x/3));

X(1)=x0;
Y(1)=y0;

for i = 1:1:length(t)-1
    
    k1x = xp(X(i),Y(i));
    k1y = yp(X(i),Y(i));
    
    k2x = xp(X(i)+(1/2)*h*k1x,Y(i)+(1/2)*h*k1y);
    k2y = yp(X(i)+(1/2)*h*k1x,Y(i)+(1/2)*h*k1y);
    
    k3x = xp(X(i)+(1/2)*h*k2x,Y(i)+(1/2)*h*k2y);
    k3y = yp(X(i)+(1/2)*h*k2x,Y(i)+(1/2)*h*k2y);
    
    k4x = xp(X(i)+h(1)*k3x,Y(i)+h*k3y);
    k4y = yp(X(i)+h(1)*k3x,Y(i)+h*k3y);
    
    X(i+1) = X(i)+h*(1/6)*(k1x+2*k2x+2*k3x+k4x);
    Y(i+1) = Y(i)+h*(1/6)*(k1y+2*k2y+2*k3y+k4y);
end
end