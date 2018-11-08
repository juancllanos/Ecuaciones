%PRIMERA PARTE 

%(1):

h = [0.1,0.05,0.02,0.01,0.005,0.001];
x = 0:h(1):1;
y0 = 1;
f = @(x,y) x+y;
w0 = y0;
W(1) = y0;
for i= 1:1:10
    W(i+1) = W(i)+h(1)*f(x(i),W(i));
end
%plot(x,W)

%(2):

E(1) = y0;
for i= 1:1:10
    E(i+1) = E(i)+h(2)*f(x(i),E(i));
end
R(1) = y0;
for i= 1:1:10
    R(i+1) = R(i)+h(3)*f(x(i),R(i));
end
T(1) = y0;
for i= 1:1:10
    T(i+1) = T(i)+h(4)*f(x(i),T(i));
end
U(1) = y0;
for i= 1:1:10
    U(i+1) = U(i)+h(5)*f(x(i),U(i));
end
O(1) = y0;
for i= 1:1:10
    O(i+1) = O(i)+h(6)*f(x(i),O(i));
end
g = @(x) -x-1+(2/exp(-x));
A=[];
for i = 1:1:11
    A(i) = g(x(i));
end

%plot(x,W,x,E,x,R,x,T,x,U,x,O,x,A)

%(3):
E=[];
for i = 1:1:11
    E(i) = abs(A(i)-O(i));
end
disp(E);

%(4):

%SEGUNDA PARTE :   

%(5):
Z(1)=y0;
for i = 1:h(1):10
    k1 = f(x(i),Z(i));
    k2 = f(x(i)+(1/2)*h(1),Z(i)+(1/2)*h(1)*k1);
    k3 = f(x(i)+(1/2)*h(1),Z(i)+(1/2)*h(1)*k2);
    k4 = f(x(i)+h(1),Z(i)+h(1)*k3);
    
    Z(i+1) = Z(i)+h(1)*(1/6)*(k1+2*k2+2*k3+k4);
end
