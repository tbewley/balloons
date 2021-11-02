clear; global L0 m g trial h C q a DeltaL clip Lold   % define parameters as globals
L0=1; m=1; g=9.8; a=1; T0=2*pi*sqrt(L0/g); q=2;
tol=0.01; relax=0.93; L=L0; Rmax=0;
for trial=5:5
  if trial==1,      C=0.0306,  DeltaL=0.01;   ls='k-';  % T/T0=99.6667
  elseif trial==2,  C=0.076,   DeltaL=0.025;  ls='r-';  % T/T0=40.2722
  elseif trial==3,  C=0.153,   DeltaL=0.05;   ls='b-';  % T/T0=20.1386
  elseif trial==4,  C=0.305,   DeltaL=0.1;    ls='m-';  % T/T0=10.0743
  elseif trial==5,  C=0.61,    DeltaL=0.2;    ls='k-';  % T/T0=5.0421
  else
     disp('Case not yet implemented.'); break;
  end
  figure(1); clf;
  x(1)=1; x(2)=1e-3; h=0.01; Tmax=300; t=0; L_save(1)=L0  % x(1)=phi, x(2)=phidot
  for k=1:Tmax/h+1                      % RK4 simulation of pendulum motion
    t_save(k)=t/T0; x_save(k,:)=x; Lold=L;
    [phidotdot,L,Ldot]=calulate_evolution(x(1),x(2));  L_save(k)=L;
    if k>1 & x_save(k,2)*x_save(k-1,2)<0 & abs(x(1))<0.01, break, end
    if k>1 & x_save(k,1)*x_save(k-1,1)<0, Rmax=max(Rmax,(x_save(k,1)-x_save(k-1,1))/h); end
    f1=RHS(x); f2=RHS(x+h*f1/2); f3=RHS(x+h*f2/2); f4=RHS(x+h*f3);
    x=x+h*(f1/6+(f2+f3)/3+f4/6); t=t+h;
  end
  Tmax_over_T0=t/T0, Rmax_over_R0=Rmax/(L0/T0)
  subplot(2,1,1), plot(t_save(:),x_save(:,1),ls); hold on; plot([0 t/T0],[0 0],'k-.');   ylabel('p'), title('phi(t)')
  subplot(2,1,2), plot(t_save(:),L_save(:),ls);   hold on; plot([0 t/T0],[L0 L0],'k-.'); ylabel('l'), title('L(t)')
pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=RHS(x)
ans(1)=x(2);
ans(2)=calulate_evolution(x(1),x(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phidotdot,L,Ldot]=calulate_evolution(phi,phidot)
global L0 m g C q h a DeltaL clip Lold

Vs0=0.5*m*L0^2*phidot^2 + m*g*L0*(1-cos(phi));
L=L0*(1+phi*phidot*C/Vs0^(q/2));
Ldot = (2*C*(L0^3*m*phidot^6 + 2*C*L0*phi*phidot^5 + 4*g^3*m*phi*sin(2*phi) - g^3*m*phi*sin(3*phi) + 4*C*g*phi*phidot^3 + 6*L0*g^2*m*phidot^2 + 4*L0^2*g*m*phidot^4 - 5*g^3*m*phi*sin(phi) - 4*C*g*phi*phidot^3*cos(phi) - 8*L0*g^2*m*phidot^2*cos(phi) - 4*L0^2*g*m*phidot^4*cos(phi) - 4*C*g*phi^2*phidot^3*sin(phi) + 2*L0*g^2*m*phidot^2*cos(2*phi) - 4*L0*g^2*m*phi*phidot^2*sin(phi) - L0^2*g*m*phi*phidot^4*sin(phi) + 2*L0*g^2*m*phi*phidot^2*sin(2*phi)))/(m*(L0*phidot^2 + 2*g - 2*g*cos(phi))*(L0^3*m*phidot^4 + 4*L0*g^2*m + 4*L0*g^2*m*cos(phi)^2 - 2*C*L0*phi*phidot^3 + 4*L0^2*g*m*phidot^2 - 8*L0*g^2*m*cos(phi) + 12*C*g*phi*phidot - 4*L0^2*g*m*phidot^2*cos(phi) - 12*C*g*phi*phidot*cos(phi)));

if     L>L0+DeltaL, clip=true, L=L0+DeltaL, Ldot=(L-Lold)/h
elseif L<L0-DeltaL, clip=true, L=L0-DeltaL, Ldot=(L-Lold)/h
else   clip=false;
end

phidotdot = -(g*sin(phi)+2*Ldot*phidot)/L;
end
