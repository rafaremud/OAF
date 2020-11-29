%-----------------------------------------%
%   UNIVERSIDADE FEDERAL DE MATO GROSSO   %
%           INSTITUTO DE FÍSICA           %
%           FÍSICA COMPUTACIONAL          %
%        Aluno: Rafael Alves Dumer        %
%                                         &    
%      OSCILADOR AMORTECIDO FORÇADO       % 
%-----------------------------------------%
clc
clear
clf
close all

%primeira configuração: b=-0.1
%segunda configuração: b=+0.1
%terceira configuração: b=0

m=1;              %massa do objeto (Kg)
k=5;              %constante elastica da mola (N/m)
b=0;              %constante da força de atrito/arrasto (N/(m/s))
F=10;             %força impulsionadora (N)

w20=sqrt(k./m);   %frequencia angular caracteristica sem amortecimento

B = b./(2.*m);    %parametro de amortecimento

A = F./m;         %paramentro de impulso
w=5;              %frequencia angular da força de impulsão


A1 = 0.5;         %amplitude positiva da solução analitica     
A2 = 0.5;         %amplitude negativa da solução analitica

%função para a velocidade (dx/dt) = z
f = @(z) z;       
%função para a derivada da velocidade (dz/dt)
x = @(x,z,t) -2.*B.*z - w20.*x + A.*cos(w.*t);

%função oscilador amortecido
g = @(x,z) -2.*B.*z - w20.*x;
%função impulso
Forca = @(t) A.*cos(w.*t);

%Solução analitica = Solução complementar + Solução particular;
%oscilador amortecido
Sc = @(t) exp(-B.*t).*(A1.*exp(sqrt((B.^2)-(w20.^2)).*t) + A2.*exp(-sqrt((B.^2)-(w20.^2)).*t));
%força de impulsão externa
del = atan((2.*w.*B)./((w20.^2) - (w.^2)));
Sp = @(t) (A./sqrt(((w20.^2) - (w.^2)).^2 + (4.*(w.^2).*(B.^2)))).*cos(w.*t - del);
%posição
a = @(t) Sc(t) + Sp(t);
%velocidade analitiva derivada da posição
%primeira parte da derivada de Sc
DelSc1 = @(t) (-B.*exp(-B.*t).*(A1.*exp(sqrt((B.^2)-(w20^2)).*t) + A2.*exp(-sqrt((B.^2)-(w20^2)).*t)));
%segunda parte da derivada de Sc
DelSc2 = @(t) (exp(-B.*t).*(A1.*sqrt((B.^2)-(w20^2)).*exp(sqrt((B.^2)-(w20^2)).*t) - A2.*sqrt((B.^2)-(w20^2)).*exp(-sqrt((B.^2)-(w20^2)).*t)));
%parte da derivada de Sp
DelSp1 = @(t) ((-A.*w.*sin(w.*t - del))./(sqrt(((w20.^2) - (w.^2)).^2 + (4.*(w.^2).*(B.^2)))));
%velocidade
Va= @(t) DelSc1(t) + DelSc2(t) + DelSp1(t);

ind = 1; 
%condições iniciais da solução analitica
matriza(ind) = a(0);
matrizSc(ind) = Sc(0);
matrizSp(ind) = Sp(0);
matrizVa(ind) = Va(0);

%condições iniciais do método, igualadas com a analitica
z0 = Va(0);           %velocidade inicial
x0 = a(0);            %posição inicial
F0 = Sp(0);           %Força impulcionadora inicial
z10= Va(0);           %velocidade inicial para a solução amortecida
g0 = a(0);            %solução temporaria/amortecido

t0 = 0;               %tempo inicial
dt = 0.1;             %passo do tempo
tn = 50;              %tempo final


oscilacoesF = tn/(2.*pi./w);

%adicionando condições as respectivas matrizes linha de dados;   
matrizz(ind) = z0;
matrizx(ind) = x0;
matrizt(ind) = t0;
matrizg(ind) = g0;
matrizF(ind) = F0;

%loop para método Runge-Kutta de quarta ordem para a EDO de segunda ordem
for i = t0:dt:tn
  
ind = ind+1;

%runge-kutta para impulso+amortecido
f1 = dt.*f(z0);
x1 = dt.*x(x0,z0,t0);
f2 = dt.*f(z0+(x1.*0.5));
x2 = dt.*x((x0+(f1.*0.5)),(z0+(x1.*0.5)),(t0+(dt.*0.5)));
f3 = dt.*f(z0+(x2.*0.5));
x3 = dt.*g((x0+(f2.*0.5)),(z0+(x2.*0.5)),(t0+(dt.*0.5)));
f4 = dt.*f(z0+x3);
x4 = dt.*g((x0+f3),(z0+x3),(t0+dt));

%runge-kutta para oscilador amortecido
k1 = dt.*f(z10);
l1 = dt.*g(g0,z10);
k2 = dt.*f(z10+(l1.*0.5));
l2 = dt.*g((g0+(k1.*0.5)),(z10+(l1.*0.5)));
k3 = dt.*f(z10+(l2.*0.5));
l3 = dt.*g((g0+(k2.*0.5)),(z10+(l2.*0.5)));
k4 = dt.*f(z10+l3);
l4 = dt.*g((g0+k3),(z10+l3));

%runge-kutta para força de impulso
F1 = dt.*Forca(t0);
F2 = dt.*Forca((t0+(dt.*0.5)));
F3 = dt.*Forca((t0+(dt.*0.5)));
F4 = dt.*Forca((t0+dt));

zn = z0 + (1./6).*(x1 + 2.*x2 + 2.*x3 + x4); %velocidade amortecido+impulso
xn = x0 + (1./6).*(f1 + 2.*f2 + 2.*f3 + f4); %posição amortecido+impulso 
z1n= z10+ (1./6).*(l1 + 2.*l2 + 2.*l3 + l4); %velocidade amortecido 
gn = g0 + (1./6).*(k1 + 2.*k2 + 2.*k3 + k4); %posição amortecido 
Fn = F0 + (1./6).*(F1 + 2.*F2 + 2.*F3 + F4); %força de impulsão externa
tn = t0 + dt;

z0 = zn;
x0 = xn;
z10= z1n;
g0 = gn;
F0 = Fn;
t0 = tn;

matrizz(ind) = z0;
matrizx(ind) = x0;
matrizt(ind) = i;
matrizg(ind) = g0;
matrizF(ind) = F0;

matriza(ind) = a(i+dt);
matrizSc(ind) = Sc(i+dt);
matrizSp(ind) = Sp(i+dt);
matrizVa(ind) = Va(i+dt);

endfor

%Poincaré
matPCx(1) = 0;    %dados Poincaré eixo x analitico
matPCy(1) = 0;    %dados Poincaré eixo y analitico
matPCmx(1) = 0;   %dados Poincaré eixo x RK4
matPCmy(1) = 0;   %dados Poincaré eixo y RK4
%indices
ond = 0;
und = -1;
ant = 0;
%limitador do loop para não passar a quantidade de dados
tnn = oscilacoesF-2.*pi;
%loop para capturar ponto no periodo de oscilação da força externa
for j = (2.*pi./w):(2.*pi./w):tnn
und = und+1;
%loop de contagem para captar dados da posição ond
for k = ((2.*pi./w).*und):dt:j
ond = ond+1;
endfor
%armazenando ponto capturado na matriz de dados de Poincaré
ant = ant+1;
matPCx(ant) = matriza(ond);
matPCy(ant) = matrizVa(ond);
matPCmx(ant) = matrizx(ond);
matPCmy(ant) = matrizz(ond);
endfor

%plotar diagrama de fase tridimensional
subplot(2,2,1);
plot3(matrizt,matrizx,matrizz,'g-');
xlabel('t'); ylabel('x(t)'); zlabel('v(t)');
title('Diagrama de fase tridimensional RK4');
subplot(2,2,2);
plot3(matrizt,matriza,matrizVa,'k-');
xlabel('t'); ylabel('x(t)'); zlabel('v(t)');
title('Diagrama de fase tridimensional analitico');
%plotar diagrama de fase bidimensional
subplot(2,2,3);
plot(matrizx,matrizz,'g-');
xlabel('x(t)'); ylabel('v(t)');
title('Diagrama de fase bidimensional RK4');
subplot(2,2,4);
plot(matriza,matrizVa,'k-');
xlabel('x(t)');
ylabel('v(t)');
title('Diagrama de fase bidimensional analitico');

figure
%plotar posição em função do tempo
subplot(2,1,1);
plot(matrizt,matrizx);
hold on
plot(matrizt,matrizg);
hold on
plot(matrizt,matrizF,'y-');
xlabel('t');
ylabel('x(t)');
legend('amortecido forçado','amorteciemnto','força');
title('Posição do oscilador RK4');
subplot(2,1,2);
plot(matrizt,matriza);
hold on
plot(matrizt,matrizSc);
hold on
plot(matrizt,matrizSp,'y-');
xlabel('t');
ylabel('x(t)');
legend('amortecido forçado','amorteciemnto','força');
title('Posição do oscilador analitico');

figure
%plotar seção de Poincaré
subplot(2,1,1)
plot(matPCmx,matPCmy,'g*');
xlabel('x(t)');ylabel('v(t)');
title('Seção de Poincaré RK4');
subplot(2,1,2)
plot(matPCx,matPCy,'k*');
xlabel('x(t)');ylabel('v(t)');
title('Seção de Poincaré analitico');



