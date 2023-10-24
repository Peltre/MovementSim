clear all
close all


% Coordenadas Iniciales1
x1 = 10;
y1 = 170;

% Cordenadas Finales
x2 = 280;
y2 = 20;

% Suposición de coordenadas
x3 = 100;
y3 = 75;
x4 = 130;
y4 = 100;

Y = [y1; y2; y3; y4];
X = [x1^3 x1^2 x1 1; x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1];
A = inv(X)*Y;

x = linspace(x1,x2,201);

fcurva = A(1)*x.^3 + A(2)*x.^2 + A(3)*x + A(4);


% Definición de vectores iniciales
dx = [];
dy = [];
ds = [];
der2x = [];
der2y =[];
R = [];
c = [];

% longitud de curva y primera derivada

for i = 1:length(x)-1
    dx = [dx,x(i+1)-x(i)];
    dy = [dy,fcurva(i+1)-fcurva(i)];
    ds = [ds,sqrt(dx(i)^2 + dy(i)^2)];
    der1x = dx./ds;
    der1y = dy./ds;
end

% segunda derivada y radio de curvatura

for i = 1:length(der1x)-1
    der2x = [der2x, (der1x(i+1)-der1x(i))/ds(i)];
    der2y = [der2y, (der1y(i+1)-der1y(i))/ds(i)];
% Radio de curvatura
    R = [R,1/sqrt(der2x(i).^2+der2y(i).^2)];
end

longitud = sum(ds);
Radmin = min(R);

disp (['La longitud de la curva es: ', num2str(longitud)])
disp (['El radio de curvatura es: ', num2str(Radmin)])

% Aquí buscamos los puntos donde la primera derivada cambia de signo y pasa por su valor mas cercano a 0:
for j = 2:length(der1y)-1
    if (der1y(j)>0 && der1y(j-1)<0) | (der1y(j)>0 && der1y(j+1)<0)
        c=[c,j];% Guarda los índices donde esto ocurre
    end
end

% En el radio buscamos dónde alrededor de los índices en c, el radio es
% mínimo
indice1 = find(R==min(R(c(1)-5:c(1)+5)));
indice2 = find(R==min(R(c(2)-5:c(2)+5)));
R1 = R(indice1);
R2 = R(indice2);

% % Buscamos donde el valor de x a la izquierda desplazado aproximadamente un radio mínimo:
% primer tramo
t1=find(x>(x(indice1)-R1-0.6) & x<(x(indice1)-R1+0.6));
t2=find(x>(x(indice2)-R2-0.6) & x<(x(indice2)-R2+0.6));% Este debe ser solo un índice

% % HALLAMOS LA ECUACIÓN DE LA RECTA TANGENTE A ESTE PUNTO:
% % Recta 1
m1 = der1y(t1)/der1x(t1);
b1 = fcurva(t1) - m1*x(t1); % Punto de corte
coorx1 = [x(t1-20):x(t1+20)];
tangente1 = m1*coorx1 + b1;% REcta tangente al punto critico
%Recta perpendicular a la recta tangente
perpendicular1 = -(1/m1)*coorx1 + (fcurva(t1)+x(t1)/m1);


%Queremos una grada a 20 m de esta recta con un ancho de 10 m y un largo de
%80 m:
paralela1 = tangente1-20;% Recta a 20 unidades de la recta tangente al punto critico
crucex = ((fcurva(t1)+x(t1)/m1)-(b1-20))/(m1+(1/m1));
crucey = m1*crucex + b1-20;

paralela2 = tangente1-30;
crucex2 = ((fcurva(t1)+x(t1)/m1)-(b1-30))/(m1+(1/m1));
crucey2 = m1*crucex2 + b1-30;

distancia1 = sqrt((coorx1-crucex).^2+(paralela1-crucey).^2);
distancia2 = sqrt((coorx1-crucex).^2+(paralela2-crucey).^2);
largo1 = find(distancia1<21 & distancia1>19);% Las condiciones deberán ser tales que devuelva dos valores
largo2 = find(distancia2<21 & distancia2>19);% Las condiciones deberán ser tales que devuelva dos valores
gradax1 = coorx1(largo1(1));
graday1 = paralela1(largo1(1));
gradax2 = coorx1(largo1(2));
graday2 = paralela1(largo1(2));

gradax3 = coorx1(largo2(1));
graday3 = paralela2(largo2(1));
gradax4 = coorx1(largo2(2));
graday4 = paralela2(largo2(2));

% % HALLAMOS LA ECUACIÓN DE LA RECTA TANGENTE A ESTE PUNTO:
% % Recta 2
m2 = der1y(t2)/der1x(t2);
b2 = fcurva(t2) - m2*x(t2); % Punto de corte
coorx2 = [x(t2-35):x(t2+35)];
tangente2 = m2*coorx2 + b2;% REcta tangente al punto critico
%Recta perpendicular a la recta tangente
perpendicular12 = -(1/m2)*coorx2 + (fcurva(t2)+x(t2)/m2);

%Queremos una grada a 20 m de esta recta con un ancho de 10 m y un largo de
%80 m:
paralela12 = tangente2-20;% Recta a 20 unidades de la recta tangente al punto critico
crucex12 = ((fcurva(t2)+x(t2)/m2)-(b2-20))/(m2+(1/m2));
crucey12 = m2*crucex12 + b2-20;

paralela22 = tangente2-30;
crucex22 = ((fcurva(t2)+x(t2)/m2)-(b2-30))/(m2+(1/m2));
crucey22 = m2*crucex22 + b2-30;

distancia11 = sqrt((coorx2-crucex12).^2+(paralela12-crucey12).^2);
distancia22 = sqrt((coorx2-crucex12).^2+(paralela22-crucey12).^2);
largo11 = find(distancia11<23 & distancia11>21);% Las condiciones deberán ser tales que devuelva dos valores
largo22 = find(distancia22<23 & distancia22>21);% Las condiciones deberán ser tales que devuelva dos valores
gradax12 = coorx2(largo11(1));
graday12 = paralela12(largo11(1));
gradax22 = coorx2(largo11(2));
graday22 = paralela12(largo11(2));

gradax32 = coorx2(largo22(1));
graday32 = paralela22(largo22(1));
gradax42 = coorx2(largo22(2));
graday42 = paralela22(largo22(2));

%Calculemos la velocidad maxima (Sin peralte)
mk = 0.9;
g = 9.81;
vmax = sqrt(mk*g*R(t1));
vmax2 = sqrt(mk*g*R(t2));
xmax = [x(t1:indice1)];
xmax2 = [x(t2:indice2)];

%Velocidad maxima (con peralte)
angulo = 8;
angulor = 0.1396;
vmaxp = sqrt((g*R(t1))*((sin(angulor))+(mk*cos(angulor)))/(cos(angulor)-(mk*sin(angulor))));
vmaxp2 = sqrt((g*R(t2))*((sin(angulor))+(mk*cos(angulor)))/(cos(angulor)-(mk*sin(angulor))));

%Calculemos la distancia sin peralte
dist1 = (vmax.^2)/(2*mk*g);
dist2 = (vmax2.^2)/(2*mk*g);

%Calculemos la distancia con peralte
dist3 = (vmaxp.^2)/(2*mk*g);
dist4 = (vmaxp2.^2)/(2*mk*g);

%Preguntar al usuario para inicializar la simulacion
disp("Bienvenid@ al simulador de pista de f1")
disp("Desea hacer una simulacion con o sin peralte?")
disp("Introduzca 1 para simular sin peralte y 2 para simular con peralte")
respuesta = input("Simulacion con o sin peralte?: ");
vu = input("Introduzca una velocidad en m/s:  ");
mk = 0.9;
g = 9.81;
masa = 798;

parametro = 0;

if respuesta == 1
    if (vu < vmax) && (vu < vmax2)
        parametro = 1;

    elseif (vu >= vmax )
        parametro = 2;
        d = vmax^2/(2*mk*g);
        calor = masa*mk*g*d;

    elseif (vu < vmax) && (vu >= vmax2)
        parametro = 3;
        d = vmax2^2/(2*mk*g);
        calor = masa*mk*g*d;
    end

elseif respuesta == 2
      if (vu < vmaxp) && (vu <vmaxp2)
        parametro = 1;

    elseif (vu >= vmaxp )
        parametro = 2;
        d = vmaxp^2/(2*mk*g);
        calor = masa*mk*g*d;

    elseif (vu < vmaxp) && (vu >= vmaxp2)
        parametro = 3;
        d = vmaxp2^2/(2*mk*g);
        calor = masa*mk*g*d;
      end
end

%Plot dinamico
pause_time = 0.1;
dim = length(x);

figure(5),clf
hold on

%Dibujar la pista y las gradas.
plot(x,fcurva,'LineWidth',15,'Color','k')
plot(x,fcurva,'LineWidth',2,'Color','y')
plot(coorx1,tangente1,'--','LineWidth',2,'Color','red')
patch('Faces',[1,2,3,4],'Vertices',[gradax3 graday3;gradax4 graday4;gradax2 graday2;gradax1 graday1])
plot(coorx2,tangente2, '--','LineWidth',2,'Color','red')
patch('Faces',[1,2,3,4],'Vertices',[gradax32 graday32;gradax42 graday42;gradax22 graday22;gradax12 graday12])
xlabel('x(m)')
ylabel('y(m)')

txt1 = "Sobrepasó el limite de velocidad";
txt2 = "Calor disipado: ";
txt3 = "kJ";

point = 0;

%Movimiento del objeto
for mov = 1:dim
    if parametro == 1
        point = scatter(x(mov),fcurva(mov),50,'cyan',"filled");
        pause = pause_time;
        drawnow
        delete(point)

    elseif parametro == 2
        if x(mov)>= x(t1)
            text(100,20,txt1)
            text(100,40,[txt2 num2str(calor/1000, '%g') txt3])
            point = scatter(x(mov),m1*x(mov) + b1,40,'cyan','filled');
            pause = pause_time;
            drawnow
            delete(point)

            if sqrt((x(mov)-x(t1))^2 + (m1*(x(mov)-x(t1)))^2)>=d
                break
            else
                continue
            end
        else
            point = scatter(x(mov),fcurva(mov),50,'cyan',"filled");
            pause = pause_time;
            drawnow
            delete(point)

        end

    elseif parametro == 3
        if x(mov)>= x(t2)
            text(100,20,txt1)
            text(100,40,[txt2 num2str(calor/1000, '%.5g') txt3])
            point = scatter(x(mov),m2*x(mov) + b2,40,'cyan','filled');
            pause = pause_time;
            drawnow
            delete(point)
            

            if sqrt((x(mov)-x(t2))^2 + (m2*(x(mov)-x(t2)))^2)>=d
                break
            else
                continue  
            end

        else
            point = scatter(x(mov),fcurva(mov),50,'cyan',"filled");
            pause = pause_time;
            drawnow
            delete(point)
        end
    end
end
