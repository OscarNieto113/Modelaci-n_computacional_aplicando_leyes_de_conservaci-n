%Proyecto 
%Modelación computacional aplicando leyes de conservación (Gpo 1)
%Simulación 2D de un auto de carreras en una pista
%Integrantes del equipo 5:
%A01703936 Samantha Daniela Guanipa Ugas
%A01703556 Ariann Fernando Arriaga Alcántara
%A01703523 Salvador Santana Blanco
%A01705090 Oscar Eduardo Nieto Espitia

format long
xin=300;
yin=1100;
xfin=2800;
yfin=1200;
xin2cl=190;
yin2cl=630;
xfin2cl=800;
yfin2cl=1000;
r1=130;
r2=100;
c1=[xin+xin2cl, yin+yin2cl];
c2=[xfin-xfin2cl, yfin-yfin2cl];
p1=[xin,yin];
p2=[xin+xin2cl,yin+yin2cl+r1];
p3=[((xin+xin2cl)+(xfin-xfin2cl))/2,((yin+yin2cl+r1)+(yfin-yfin2cl-r2))/2];
p4=[xfin-xfin2cl,yfin-yfin2cl-r2];
p5=[xfin,yfin];
x=[p1(1),p2(1),p3(1),p4(1),p5(1)];
y=[p1(2),p2(2),p3(2),p4(2),p5(2)];

p=polyfit(x,y,4);
x2=linspace(300,2800,25000);
y2=polyval(p,x2);


fprintf("Simulación de cinemática de vehículo \n");
deltaX=10;
i=1;
Ts=.1;
fprintf("La velocidad máxima del vehículo es de 83 m/s \n");
v=input("Coloca la velocidad inicial del vehículo (m/s) \n");
xpos=xin;
ypos=polyval(p,xpos);
t=0;



Dife=polyder(p);
Dife2=polyder(Dife);
cf=.5;

x = [450 530 530 450];
y = [2067.16 2067.16 2077.16 2077.16];

x3 = [1841 1921 1921 1841];
y3 = [47 47 37 37];

tnew=0;
Masa=600;
g=9.81;
N=Masa*g;
Fr=N*cf;
Et=0;


fprintf("1.Simulación de vehículo sin desaceleración \n");
fprintf("2.Simulación de vehículo con perfiles de aceleración y desaceleración \n");
fprintf("3.Salir del programa \n");
Rs=input("Que simulación deseas hacer? (1-3) \n");

if Rs==1
Ap= [.3,0,-.6];
A=0;
vmaxcarro=83;

while (xpos<=xfin)  
xfut=xpos+deltaX;
yfut=polyval(p,xfut);
Angulo=atan2((yfut-ypos),(xfut-xpos));

Vx=v*cos(Angulo);
Vy=v*sin(Angulo);    
Vxanterior=Vx;
Vyanterior=Vy;
xposant=xpos;
yposant=ypos;
xpos = xpos + (Ts * Vxanterior) ;
ypos = ypos + (Ts * Vyanterior) ;

longi=sqrt((xpos-xposant)^2+(ypos-yposant)^2);
Et=Et+Fr*longi;

primer=polyval(Dife,xpos);
segundo=polyval(Dife2,xpos);
R = double(abs(((1+(primer^2))^(3/2))/segundo));
v=sqrt((Vx^2)+(Vy^2));
vmax=double(sqrt(cf*R*9.81));
deltaV=vmax-v;

if v >= vmaxcarro
    As=Ap(2);
    v=v+As;
else
    As=Ap(1);
    v=v+As;
end

if v>vmax
    while (tnew<5)
       newang=atan(primer);
       Vx=v*cos(newang);
       Vy=v*sin(newang);    
       Vxanterior=Vx;
       Vyanterior=Vy;
       xposant=xpos;
       yposant=ypos;
       xpos = xpos + (Ts * Vxanterior) ;
       ypos = ypos + (Ts * Vyanterior) ;
       v=sqrt((Vx^2)+(Vy^2));
       
       longi=sqrt((xpos-xposant)^2+(ypos-yposant)^2);
       Et=Et+Fr*longi;
       
       plot(x2,y2,'black')
       grid
       xlim([0,3500])
       ylim([0,2500])

       t2p= text(1300,2200,['Perdida de calor acumulado (J): ' num2str(Et)]);
       tp1=text(1300,2100,['Posición en x (m) : ' num2str(xpos)]);
       tp2=text(1300,2000,['Posición en y (m) : ' num2str(ypos)]);
       tp3=text(1300,1900,['Velocidad (m/s): ' num2str(v)]);
       tp4=text(1300,1800,['Tiempo actual (s): ' num2str(t)]);
       tp5=text(1300,1700,'Error critico, curso fuera de pista');
       
       hold on
       patch(x,y,'blue')
       patch(x3,y3,'red')
       plot(xpos,ypos,'r*')
       drawnow
  
       hold off
       tnew=tnew+Ts;
       t=t+Ts;
    end
    break
end
t=t+Ts;

plot(x2,y2,'black')
grid
xlim([0,3500])
ylim([0,2500])
tp= text(1300,2200,['Perdida de calor acumulado (J): ' num2str(Et)]);
tp1=text(1300,2100,['Posición en x (m) : ' num2str(xpos)]);
tp2=text(1300,2000,['Posición en y (m) : ' num2str(ypos)]);
tp3=text(1300,1900,['Velocidad (m/s): ' num2str(v)]);
tp4=text(1300,1800,['Tiempo actual (s): ' num2str(t)]);
tp5=text(1300,1700,['Aceleración actual (m/s^2) : ' num2str(As/Ts)]);

hold on
patch(x,y,'blue')
patch(x3,y3,'red')
plot(xpos,ypos,'r*')
drawnow
hold off
i=i+1;

end 
title('Simulación de vehículo')
xlabel('Desplazamiento x, (m)')
ylabel('Desplazamiento y, (m)')
elseif Rs==2
Ap= [.3,0,-.65];
A=0;
vmaxcarro=83;

while (xpos<=xfin)  
xfut=xpos+deltaX;
yfut=polyval(p,xfut);
Angulo=atan2((yfut-ypos),(xfut-xpos));

Vx=v*cos(Angulo);
Vy=v*sin(Angulo);    
Vxanterior=Vx;
Vyanterior=Vy;
xposant=xpos;
yposant=ypos;
xpos = xpos + (Ts * Vxanterior) ;
ypos = ypos + (Ts * Vyanterior) ;

longi=sqrt((xpos-xposant)^2+(ypos-yposant)^2);
Et=Et+Fr*longi;

primer=polyval(Dife,xpos);
segundo=polyval(Dife2,xpos);
R = double(abs(((1+(primer^2))^(3/2))/segundo));
v=sqrt((Vx^2)+(Vy^2));
vmax=double(sqrt(cf*R*9.81));
deltaV=vmax-v;

if deltaV<=50 && v>20
    As=Ap(3);
    v=v+As;
elseif v<10
    As=Ap(1);
    v=v+As;
elseif v >= vmaxcarro
    As=Ap(2);
    v=v+As;
else
    As=Ap(1);
    v=v+As;
end

if v>vmax
    while (tnew<5)
       newang=atan(primer);
       Vx=v*cos(newang);
       Vy=v*sin(newang);    
       Vxanterior=Vx;
       Vyanterior=Vy;
       xposant=xpos;
       yposant=ypos;
       xpos = xpos + (Ts * Vxanterior) ;
       ypos = ypos + (Ts * Vyanterior) ;
       v=sqrt((Vx^2)+(Vy^2));
       
       longi=sqrt((xpos-xposant)^2+(ypos-yposant)^2);
       Et=Et+Fr*longi;
       
       plot(x2,y2,'black')
       grid
       xlim([0,3500])
       ylim([0,2500])

       t2p= text(1300,2200,['Perdida de calor acumulado (J): ' num2str(Et)]);
       tp1=text(1300,2100,['Posición en x (m) : ' num2str(xpos)]);
       tp2=text(1300,2000,['Posición en y (m) : ' num2str(ypos)]);
       tp3=text(1300,1900,['Velocidad (m/s): ' num2str(v)]);
       tp4=text(1300,1800,['Tiempo actual (s): ' num2str(t)]);
       tp5=text(1300,1700,'Error critico, curso fuera de pista');
       
       hold on
       patch(x,y,'blue')
       patch(x3,y3,'red')
       plot(xpos,ypos,'r*')
       drawnow
  
       hold off
       tnew=tnew+Ts;
       t=t+Ts;
    end
    break
end
t=t+Ts;

plot(x2,y2,'black')
grid
xlim([0,3500])
ylim([0,2500])
tp= text(1300,2200,['Perdida de calor acumulado (J): ' num2str(Et)]);
tp1=text(1300,2100,['Posición en x (m) : ' num2str(xpos)]);
tp2=text(1300,2000,['Posición en y (m) : ' num2str(ypos)]);
tp3=text(1300,1900,['Velocidad (m/s): ' num2str(v)]);
tp4=text(1300,1800,['Tiempo actual (s): ' num2str(t)]);
tp5=text(1300,1700,['Aceleración actual (m/s^2) : ' num2str(As/Ts)]);
hold on
patch(x,y,'blue')
patch(x3,y3,'red')
plot(xpos,ypos,'r*')
drawnow
hold off
i=i+1;

end 
title('Simulación de vehículo')
xlabel('Desplazamiento x, (m)')
ylabel('Desplazamiento y, (m)')
else
    fprintf("Opción no valida")
end    