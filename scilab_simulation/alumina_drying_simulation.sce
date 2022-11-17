clc 
clear 

function [dY]=upward_drying(t, Y) 
global Rhoo Cpf Cpv Cps Cpl lambda N hy W0 Tf0 Gf Rhof M0 Tf0ac Epsilon 
dY = zeros(3*N,1); 
for i=1:N 
    M(i)=Y(i); 
    W(i)=Y(N+i); 
    Tm(i)=Y(2*N+i); 
end 
for i=1:N 
    if i==1 
        dWdy(i)=(W(i)-W0)/hy; 
        dTmdy(i)=(Tm(i)-Tf0ac)/hy; 
    else 
        dWdy(i)=(W(i)-W(i-1))/hy; 
        dTmdy(i)=(Tm(i)-Tm(i-1))/hy; 
end 

    Pab(i)=(28.97/18*W(i))/(1+28.97/18*W(i))*(695.1/760); 
    Psat(i)=exp(18.3036-3816.44/(Tm(i)+227.02))/760; 
    UR(i)= Pab(i)/Psat(i); 
    // Equilibrium isotherm for alumina at 50°C
    c_iso=1.39796669; 
    k_iso=0.14982732; 
    mo_iso=1.148646602; 
    Meq_60(i)=0.0358*UR(i)/(1-0.835*UR(i))^2;
    Meq_40(i)=mo_iso*c_iso*k_iso*UR(i)/((1-k_iso*UR(i))*(1-k_iso*UR(i)+k_iso*c_iso*UR(i))); 
    Meq_50(i) = (Meq_40(i)+Meq_60(i))/2; 
    Meq(i)=Meq_60(i); 
    Rho(i)= Rhoo;
    Ep(i)=Epsilon;
    Ko=1.124342559; 
    K(i)=Ko*exp(-2129.52871/(Tm(i)+273.15));
    dMdt(i)=-K(i)*(M(i)-Meq(i)); 
    f(i)=-(1-Ep(i))*Rho(i)*dMdt(i); 
    dY(i)= dMdt(i); // Mass balance for the solid phase
    dY(i+N)=(f(i)-(Gf*dWdy(i))/(Ep(i)*Rhof)); // Mass balance for the fluid phase
    dY(i+2*N)=(-f(i)*(lambda)- (Gf*(Cpf+W(i)*Cpv)*dTmdy(i)))/(((1-Ep(i))*Rho(i)*(Cps+M(i)*Cpl))+(Ep(i)*Rhof*(Cpf+W(i)))); // Energy balance for the mixture
end 
endfunction 

function [dY]=downward_drying(t, Y) 
global Rhoo Cpf Cpv Cps Cpl lambda N hy W0 Tf0 Gf Rhof M0 Tf0des Epsilon 
dY = zeros(3*N,1); 
for i=N:-1:1 
    M(i)=Y(i); 
    W(i)=Y(N+i); 
    Tm(i)=Y(2*N+i); 
end 
for i=N:-1:1 
    if i==N 
        dWdy(i)=(W(i)-W0)/hy; 
        dTmdy(i)=(Tm(i)-Tf0des)/hy; 
    else 
        dWdy(i)=(W(i)-W(i+1))/hy; 
        dTmdy(i)=(Tm(i)-Tm(i+1))/hy; 
end 

    Pab(i)=(28.97/18*W(i))/(1+28.97/18*W(i))*(695.1/760);
    Psat(i)=exp(18.3036-3816.44/(Tm(i)+227.02))/760; 
    UR(i)= Pab(i)/Psat(i); 
    // Equilibrium isotherm for alumina at 50°C
    c_iso=1.39796669; 
    k_iso=0.14982732; 
    mo_iso=1.148646602; 
    Meq_60(i)=0.0358*UR(i)/(1-0.835*UR(i))^2; 
    Meq_40(i)=mo_iso*c_iso*k_iso*UR(i)/((1-k_iso*UR(i))*(1-k_iso*UR(i)+k_iso*c_iso*UR(i))); 
    Meq_50(i) = (Meq_40(i)+Meq_60(i))/2; 
    Meq(i)=Meq_60(i); 
    Rho(i)= Rhoo;
    Ep(i)=Epsilon;
    Ko=1.124342559;
    K(i)=Ko*exp(-2129.52871/(Tm(i)+273.15));
    dMdt(i)=-K(i)*(M(i)-Meq(i)); 
    f(i)=-(1-Ep(i))*Rho(i)*dMdt(i); 
    dY(i)= dMdt(i); // Mass balance for the solid phase
    dY(i+N)=(f(i)-(Gf*dWdy(i))/(Ep(i)*Rhof)); // Mass balance for the fluid phase
    dY(i+2*N)=(-f(i)*(lambda)- (Gf*(Cpf+W(i)*Cpv)*dTmdy(i)))/(((1-Ep(i))*Rho(i)*(Cps+M(i)*Cpl))+(Ep(i)*Rhof*(Cpf+W(i)))); // Energy balance for the mixture
end 
endfunction 

global Rhoo Cpf Cpv Cps Cpl lambda N hy W0 Tf0 Gf Rhof M0 Tf0ac Tf0des Epsilon 

// Model parameters 
Epsilon = 0.4; // Porosity of the bed
Rhoo = 1.69; // Alumina specific mass
Cpf = 0.25; // Dry air specific heat
Cpv = 0.28; // Water vapor specific heat
Cps = 0.199914; // Solid specific heat
Cpl = 1.0; // Liquid water specific heat
lambda = 573; // Latent heat of water vaporization
altura = 10.0; // Bed height
m=0.390; // Mass flow rate of air
D=10.0; // Bed diameter

// Air reversal parameters
t0rev=10*60; // Time when the first airflow reversal happens
deltrev=10*60; // Interval between airflow reversals
nrev=11; // Number of airflow reversals

// Numerical method parameters
N=11; // Number of reversions 
hyaux=linspace(0,altura,N); 
hy=hyaux(2)-hyaux(1); 

// Experiment conditions
Ufo = 0.016; // Initial air humidity
Uso = 0.45; // Initial alumina humidity
Tfo = 60.0; // Initial bed temperature
Tf0ac=Tfo; // Initial inlet temperature of the upward flow
Tf0des= Tfo; // Initial inlet temperature of the downward flow 
Tmo = 20.8; // Initial mixture of fluid and solid phases temperature

W0=Ufo; 
M0=Uso; 
Rhof =1.2e-3*293.15/(Tfo+273.15)*(1/(1+Ufo)); // Fluid density
Gf=m*1000/(60*%pi*(D^2)/4); // Mass flux of the fluid phase
t0=0.001; // Initial drying time
delt=20; // Simulation time step
tf=t0rev; // Final drying time
y0=[Uso*ones(1,N) Ufo*ones(1,N) Tmo*ones(1,N)]; // Simulation initial state
Xd=y0; 

B=[0]; 

for e=1:nrev 
    if pmodulo(e,2)==1 // Upward flow
        tempo=[t0:delt:tf];
        [Y] = ode('stiff',y0',t0,tempo,ones(1,3*N)*1e-10,ones(1,3*N)*1e-10,upward_drying); 
        y0=Y'($,:); 
        t0=tf; 
        tf=tf+deltrev; 
    end
    if pmodulo(e,2)==0 // Downward flow
        tempo=[t0:delt:tf]; 
        [Y] = ode('stiff',y0',t0,tempo,ones(1,3*N)*1e-10,ones(1,3*N)*1e-10,downward_drying); 
        y0=Y'($,:); 
        t0=tf; 
        tf=tf+deltrev; 
    end
    Xd = [Xd;Y'(:,:)]; 
    B=[B;tempo']; 
end 

B=B/60; // Time steps values in minutes

figure(1) 
mtlb_hold 
plot(B,Xd(:,2),'k.-',B,Xd(:,6),'g.-',B,Xd(:,11),'b.-.') 
xlabel('Tempo [minutos]') 
ylabel('Umidade da particula em base seca[kg água/kg sólido seco]') 
title('Umidade da alumina em função do tempo') 
legend('X(1cm)','X(5cm)','X(10cm)') 
mtlb_hold 

figure(2) 
mtlb_hold 
plot(B,Xd(:, 2*N+2),'k.-',B,Xd(:, 2*N+6),'g.-',B,Xd(:, 2*N+11),'b.-.'); 
xlabel('Tempo [min]') 
ylabel('Temperatura da mistura [°C]'); 
title('Temperatura da mistura em função do tempo') 
legend('1 cm','5 cm','10 cm') 

mtlb_hold 
figure(3) 
mtlb_hold 
plot(B,Xd(:, N+2),'k.-',B,Xd(:, N+6),'g.-',B,Xd(:, N+11),'b.-.'); 
xlabel('Tempo [min]') 
ylabel('Umidade absoluta do ar [kg água/ kg ar]'); 
title('Umidade absoluta do ar em função do tempo') 
legend('1 cm','5 cm','10 cm') 
mtlb_hold 

// Uncomment the `disp` commands below if you want to display the results
// disp('tempo X_1cm X_5cm X_10cm') 
// disp([B' Xd'(2,:) Xd'(6,:) Xd'(11,:)]) 

// disp('tempo T_1cm T_5cm T_10cm') 
// disp([B' Xd'(2*N+2,:) Xd'(2*N+6,:) Xd'(2*N+11,:)])
 
csvWrite([B/60 Xd(:,:)], 'simulation_results_scilab.csv')