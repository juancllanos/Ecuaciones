 %PRIMERA PARTE 

%(1):

h = [0.1,0.05,0.02,0.01,0.005,0.001];
x1 = 0:h(1):1;
x2 = 0:h(2):1;
x3 = 0:h(3):1;
x4 = 0:h(4):1;
x5 = 0:h(5):1;
x6 = 0:h(6):1;
y0 = 1;
f = @(x,y) x+y;
w0 = y0;
W(1) = y0;
for i= 1:1:1/h(1)
    W(i+1) = W(i)+h(1)*f(x1(i),W(i));
end
%plot(x1,W)

%(2):

E(1) = y0;
for i= 1:1:1/h(2)
    E(i+1) = E(i)+h(2)*f(x2(i),E(i));
end
R(1) = y0;
for i= 1:1:1/h(3)
    R(i+1) = R(i)+h(3)*f(x3(i),R(i));
end
T(1) = y0;
for i= 1:1:1/h(4)
    T(i+1) = T(i)+h(4)*f(x4(i),T(i));
end
U(1) = y0;
for i= 1:1:1/h(5)
    U(i+1) = U(i)+h(5)*f(x5(i),U(i));
end
O(1) = y0;
for i= 1:1:1/h(6)
    O(i+1) = O(i)+h(6)*f(x6(i),O(i));
end
g = @(x) -x-1+(2/exp(-x));
A=[];
for i = 1:1:11
    A(i) = g(x(i));
end

%plot(x1,W,x2,E,x3,R,x4,T,x5,U,x6,O,x1,A);


%(3):
%E=[];
%for i = 1:1:11
%    E(i) = abs(A(i)-O(i));
%end
%disp(E);

%(4):

%SEGUNDA PARTE :   

%(5):


Q(1) = y0;
for i = 1:1:10
    disp(i);
    k1 = f(x1(i),Q(i));
    k2 = f(x1(i)+(1/2)*h(1),Q(i)+(1/2)*h(1)*k1);
    k3 = f(x1(i)+(1/2)*h(1),Q(i)+(1/2)*h(1)*k2);
    k4 = f(x1(i)+h(1),Q(i)+h(1)*k3);
    
    Q(i+1) = Q(i)+h(1)*(1/6)*(k1+2*k2+2*k3+k4);
end
disp(length(Q)+ " <---> " + length(x1));
plot(x1,W,x1,Q);

%(6):