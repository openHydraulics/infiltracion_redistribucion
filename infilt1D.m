clear;clc
close all
%Cálculo de la infiltración vertical. Integración numérica de la Ec. de Richards mediante diferencias finitas
hum=[];%proporción humedad suelo
K=[];%conductividad hidráulica en subsaturación
tiempo=[];
uz=[];%componete vertical de la velocidad

%Datos
porosidad=0.2;
coef_curva_ret=2;
z0=-1.0;%Posición comienzo horizonte poros gruesos
Ks=4e-5;%Conductividad hidráulica en saturación
coef_curva_conduc=5;

Dz=-0.01;
z=0:Dz:z0;%A efectos de estabilidad, si se reduce el intervalo de z, hay que reducir también el de t
Dt=0.5;
tfinal=3600;
nintervtiempo=floor(tfinal/Dt);


%% variables(t,z). Las variables son matrices en la que cada fila se corresponde con un instante y cada columna con una profundidad

%Condición inicial
H(1,:)=z0.*ones(1,numel(z));
h(1,:)=H(1,:)-z;%Equilibrio hidrostático
hum(1,:)=porosidad.*exp(coef_curva_ret.*h(1,:));%Perfil de humedad correspondiente

%Condición de contorno (ver en la iteración késima para el tiempo)
H(1,1)=0;
h(1,1)=H(1,1)-z(1);

K(1,:)=Ks.*exp(coef_curva_conduc.*h(1,:));%Perfil de conductividad hidráulica correspondiente

%Cálculo de la velocidad
uz(1,:)=-(K(1,1:end-1)+K(1,2:end))./2-(K(1,1:end-1)+K(1,2:end))./2.*(h(1,2:end)-h(1,1:end-1))./(hum(1,2:end)-hum(1,1:end-1)).*(hum(1,2:end)-hum(1,1:end-1))./(z(2:end)-z(1:end-1));

subplot(2,2,1)
plot(hum(1,:),z)
hold on

tiempo(1)=0;
%Cálculo de la evolución
for k=2:1:nintervtiempo
  tiempo(k)=tiempo(k-1)+Dt;
  %Condición de contorno
  %if tiempo(k)<2000
    hum(k,1)=porosidad;
  %else
  %  ;
  %endif
  hum(k,end)=porosidad;
  
  %Cálculo de la humedad a partir de la divergencia de la velocidad
  hum(k,2:end-1)=hum(k-1,2:end-1)-Dt./Dz.*(uz(k-1,2:end)-uz(k-1,1:end-1));
  h(k,:)=1/coef_curva_ret.*log(hum(k,:)./porosidad); 
  K(k,:)=Ks.*exp(coef_curva_conduc.*h(k,:));
  uz(k,:)=-(K(k,1:end-1)+K(k,2:end))./2-(K(k,1:end-1)+K(k,2:end))./2.*(h(k,2:end)-h(k,1:end-1))./(z(2:end)-z(1:end-1));
  
  %Representación de algunas curvas de perfil de humedad
  if tiempo(k)<60 && (mod(k,50))==0
    plot(hum(k,:),z)
    disp(tiempo(k));
  elseif tiempo(k)>60 && (mod(k,500))==0
    plot(hum(k,:),z)
    disp(tiempo(k));
  endif
endfor


xlabel('humedad')
ylabel('z(m)')
hold off

ia=-cumsum(uz(:,1));
percol=-cumsum(uz(:,end));

subplot(2,2,2)
[ax,lines1,lines2]=plotyy(tiempo,-uz(:,1),tiempo,ia);
xlabel('t(s)')
ylabel(ax(1), 'infilt.instant.(m/s)','color','b')
ylabel(ax(2), 'infilt.acum.(m)','color','r')
set(ax(1),'Xlim',[0 tfinal])
set(ax(1),'Ylim',[0 0.0002])
set(ax(2),'Xlim',[0 tfinal])
set(ax(2),'Ylim',[0 0.4])
hold on
plot([tiempo(1) tiempo(end)],[Ks Ks],'b--')
hold off

subplot(2,2,3)
plot(tiempo,hum(:,1:5:50))
xlabel('t(s)')
ylabel('humedad')

subplot(2,2,4)
plot(tiempo,ia)
xlabel('t(s)')
ylabel('infilt.acum., percolac.(m)')
hold on
plot(tiempo,percol)