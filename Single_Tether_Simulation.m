clear all; close all;
x(1)=1; x(2)=0; h=0.02; Tmax=20;      % x(1)=phi, x(2)=phidot
p.L0=1; p.m=1; p.g=9.8; p.C=0.2; t=0;  % C<=0.2 converges, C>0.2 leads to divergence
for k=1:Tmax/h+1
  figure(1);  plot(t,x(1),'k*');   hold on; plot(t,Vs(x,p),'gx'); plot([0 Tmax],[0 0],'k--');       title('black is phi(t), green is Vs(t)')
  figure(2);  plot(t,L(x,p),'ro'); hold on;                       plot([0 Tmax],[p.L0 p.L0],'k--'); title('L(t) [horizontal line is at L0')
  f1=RHS(x,p); f2=RHS(x+h*f1/2,p); f3=RHS(x+h*f2/2,p); f4=RHS(x+h*f3,p);
  x=x+h*(f1/6+(f2+f3)/3+f4/6); t=t+h;
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=RHS(x,p)
L0=p.L0; m=p.m; g=p.g; C=p.C; phi=x(1); phidot=x(2);
ans(1)=phidot;
% ans(2)=-(5*g^3*m*sin(phi) + 4*C*L0*phidot^5 - 4*g^3*m*sin(2*phi) + g^3*m*sin(3*phi) + 8*C*g*phidot^3 - 8*C*g*phidot^3*cos(phi) - 8*C*g*phi*phidot^3*sin(phi) + 4*L0*g^2*m*phidot^2*sin(phi) + L0^2*g*m*phidot^4*sin(phi) - 2*L0*g^2*m*phidot^2*sin(2*phi))/(L0^3*m*phidot^4 + 4*L0*g^2*m + 4*L0*g^2*m*cos(phi)^2 - 2*C*L0*phi*phidot^3 + 4*L0^2*g*m*phidot^2 - 8*L0*g^2*m*cos(phi) + 12*C*g*phi*phidot - 4*L0^2*g*m*phidot^2*cos(phi) - 12*C*g*phi*phidot*cos(phi));
ans(2)=-(g*sin(phi)+2*Ldot(x,p)*phidot)/L(x,p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=L(x,p)                   % MUHAN: this is the new control rule I propose
L0=p.L0; m=p.m; g=p.g; C=p.C; phi=x(1); phidot=x(2);
Vs0=0.5*m*L0^2*phidot^2 + m*g*L0*(1-cos(phi));
ans=p.L0*(1+phi*phidot*C/Vs0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=Vs(x,p)
m=p.m; g=p.g; phi=x(1); phidot=x(2);
L1=L(x,p);
ans = 0.5*m*(L1^2*phidot^2)+m*g*L1*(1-cos(phi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=Ldot(x,p)
L0=p.L0; m=p.m; g=p.g; C=p.C; phi=x(1); phidot=x(2);
ans = (2*C*(L0^3*m*phidot^6 + 2*C*L0*phi*phidot^5 + 4*g^3*m*phi*sin(2*phi) - g^3*m*phi*sin(3*phi) + 4*C*g*phi*phidot^3 + 6*L0*g^2*m*phidot^2 + 4*L0^2*g*m*phidot^4 - 5*g^3*m*phi*sin(phi) - 4*C*g*phi*phidot^3*cos(phi) - 8*L0*g^2*m*phidot^2*cos(phi) - 4*L0^2*g*m*phidot^4*cos(phi) - 4*C*g*phi^2*phidot^3*sin(phi) + 2*L0*g^2*m*phidot^2*cos(2*phi) - 4*L0*g^2*m*phi*phidot^2*sin(phi) - L0^2*g*m*phi*phidot^4*sin(phi) + 2*L0*g^2*m*phi*phidot^2*sin(2*phi)))/(m*(L0*phidot^2 + 2*g - 2*g*cos(phi))*(L0^3*m*phidot^4 + 4*L0*g^2*m + 4*L0*g^2*m*cos(phi)^2 - 2*C*L0*phi*phidot^3 + 4*L0^2*g*m*phidot^2 - 8*L0*g^2*m*cos(phi) + 12*C*g*phi*phidot - 4*L0^2*g*m*phidot^2*cos(phi) - 12*C*g*phi*phidot*cos(phi)));
end
