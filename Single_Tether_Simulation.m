clear; global L0 m g trial fb C q a       % define parameters as globals   
L0=1; m=1; g=9.8; a=1; T0=2*pi*sqrt(L0/g);
for trial=0:5
  fb=false;
  if trial==0
     figure(1); clf; ls='k-.';  disp('Control off (delta=0)');
  elseif trial==1
     C=0.1;   q=0;   ls='k-';   disp(sprintf('Control on (after delay), delta = %g',C))     
  elseif trial==2
     C=0.05;  q=1;   ls='r-';   disp(sprintf('Control on (after delay), delta = %g/a',C))
  elseif trial==3
     C=0.05;  q=2;   ls='b-';   disp(sprintf('Control on (after delay), delta = %g/a^2',C))
  elseif trial==4
     figure(2); clf;
     C=0.14; q=1;  ls='r-';   disp(sprintf('Control on (after delay), delta = %g/[V_(s0)]^(1/2)',C))
  elseif trial==5   % take C<=0.2
     C=0.2;  q=2;   ls='b-';    disp(sprintf('Control on (after delay), delta = %g/[V_(s0)]',C))
  else
     disp('Case not yet implemented.'); break;
  end

  x(1)=1; x(2)=1e-3; h=0.01; Tmax=60; t=0; % x(1)=phi, x(2)=phidot
  for k=1:Tmax/h+1                      % RK4 simulation of pendulum motion
    t_save(k)=t/T0; x_save(k,:)=x;
    if ~fb & trial>0  & k>40 & x_save(k,1)*x_save(k-1,1)<0, fb=true; end
    if (trial==2 | trial==3) & k>1  & x_save(k,2)*x_save(k-1,2)<0, a=abs(x_save(k,1)); end
    if (trial<=3), a_save(k)=a; end
    [phidotdot,L,Ldot]        =calulate_evolution(x(1),x(2)); L_save(k)=L;
    [Vs_save(k),Vsdot_save(k)]=perform_analysis(x(1),x(2),phidotdot,L,Ldot);
    f1=RHS(x); f2=RHS(x+h*f1/2); f3=RHS(x+h*f2/2); f4=RHS(x+h*f3);
    x=x+h*(f1/6+(f2+f3)/3+f4/6); t=t+h;
  end
  if trial<=3
    subplot(4,1,1), plot(t_save(:),x_save(:,1),ls);    hold on; plot([0 Tmax],[0 0],'k-.');   axis([0 15 -1 1]);
      ylabel('p') % title('(dashed) phi(t), (solid) phidot(t)')
    subplot(4,1,2), plot(t_save(:),x_save(:,2),ls);    hold on; plot([0 Tmax],[0 0],'k-.');   axis([0 15 -3 3]);
      ylabel('d') % title('(dashed) phi(t), (solid) phidot(t)')
    subplot(4,1,3), plot(t_save(:),a_save(:),ls);      hold on; plot([0 Tmax],[0 0],'k--');   axis([0 15 0 1.1]);
      ylabel('a') % title('Vsdot(t)')
    subplot(4,1,4), plot(t_save(:),L_save(:),ls);      hold on; plot([0 Tmax],[L0 L0],'k-.'); axis([0 15 0.9 1.1]);
      ylabel('l') % title('L(t) [horizontal line is at L0')
  else    
    subplot(4,1,1), plot(t_save(:),x_save(:,2),ls);    hold on; plot([0 Tmax],[0 0],'k-.');   axis([0 28 -3 3]);
      ylabel('p') % title('(dashed) phi(t), (solid) phidot(t)')
    subplot(4,1,2), plot(t_save(:),Vsdot_save(:),ls);  hold on; plot([0 Tmax],[0 0],'k--');   axis([0 8 -4 0.2]);
      ylabel('d') % title('Vsdot(t)')
    subplot(4,1,3), semilogy(t_save(:),Vs_save(:),ls); hold on; axis([0 28 0 10]);
      ylabel('v') % title('Vs(t)')
    subplot(4,1,4), plot(t_save(:),L_save(:),ls);      hold on; plot([0 Tmax],[L0 L0],'k-.'); axis([0 28 0.9 1.1]);
      ylabel('l') % title('L(t) [horizontal line is at L0')
  end
end
figure(1); print -depsc -painters varpend_a_traces.eps
figure(2); print -depsc -painters varpend_Vs0_traces.eps
%  figure(2); hold on; plot3([x_old(1) x(1)],[x_old(2) x(2)],[Vs_old Vs_new]+0.01,'m-','linewidth',3)
%  figure(3); hold on; plot3([x_old(1) x(1)],[x_old(2) x(2)],[Vsdot_old Vsdot_new]+0.01,'m-','linewidth',3)
% figure(2); view([-9.9 67.9]), print -depsc -painters Vs_iso_traj.eps
% figure(3); view([-9.9 67.9]), print -depsc -painters Vsdot_iso_traj.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans=RHS(x)
ans(1)=x(2);
ans(2)=calulate_evolution(x(1),x(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phidotdot,L,Ldot]=calulate_evolution(phi,phidot)
global L0 m g C q fb trial a
if ~fb
    L=L0; Ldot = 0;
else
  if trial==1 | trial==2 | trial==3
    L=L0*(1+phi*phidot*C/a^q);
    Ldot=(C*L0*(L*phidot^2 - g*phi*sin(phi)))/(L*a^q + 2*C*L0*phi*phidot);
  elseif trial==4 | trial==5
    Vs0=0.5*m*L0^2*phidot^2 + m*g*L0*(1-cos(phi));
    L=L0*(1+phi*phidot*C/Vs0^(q/2));
    if trial==4
      Ldot = -(2*C*L0^2*m*(2*g^2*phi*sin(phi) - L*L0*phidot^4 - 2*L*g*phidot^2 - 2*g^2*phi*cos(phi)*sin(phi) + 2*L*g*phidot^2*cos(phi) + L*g*phi*phidot^2*sin(phi)))/(2^(1/2)*L*(L0*m*(L0*phidot^2 + 2*g - 2*g*cos(phi)))^(3/2) + 8*C*L0^2*g*m*phi*phidot - 8*C*L0^2*g*m*phi*phidot*cos(phi));
    elseif trial==5
      Ldot = (2*C*(L0^3*m*phidot^6 + 2*C*L0*phi*phidot^5 + 4*g^3*m*phi*sin(2*phi) - g^3*m*phi*sin(3*phi) + 4*C*g*phi*phidot^3 + 6*L0*g^2*m*phidot^2 + 4*L0^2*g*m*phidot^4 - 5*g^3*m*phi*sin(phi) - 4*C*g*phi*phidot^3*cos(phi) - 8*L0*g^2*m*phidot^2*cos(phi) - 4*L0^2*g*m*phidot^4*cos(phi) - 4*C*g*phi^2*phidot^3*sin(phi) + 2*L0*g^2*m*phidot^2*cos(2*phi) - 4*L0*g^2*m*phi*phidot^2*sin(phi) - L0^2*g*m*phi*phidot^4*sin(phi) + 2*L0*g^2*m*phi*phidot^2*sin(2*phi)))/(m*(L0*phidot^2 + 2*g - 2*g*cos(phi))*(L0^3*m*phidot^4 + 4*L0*g^2*m + 4*L0*g^2*m*cos(phi)^2 - 2*C*L0*phi*phidot^3 + 4*L0^2*g*m*phidot^2 - 8*L0*g^2*m*cos(phi) + 12*C*g*phi*phidot - 4*L0^2*g*m*phidot^2*cos(phi) - 12*C*g*phi*phidot*cos(phi)));
    end
  end
end
phidotdot = -(g*sin(phi)+2*Ldot*phidot)/L;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vs,Vsdot]=perform_analysis(phi,phidot,phidotdot,L,Ldot)
global m g
Vs    = 0.5*m*(L^2*phidot^2)+m*g*L*(1-cos(phi));
Vsdot = m*(L*Ldot*phidot^2+L^2*phidot*phidotdot) + m*g*(Ldot*(1-cos(phi)) + L*phidot*sin(phi));
end