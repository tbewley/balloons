% script varpend_table2.m, by Thomas Bewley
% Simulate variable length pendulum with Algorithm 1 with fillets.

% Physical constants
clear all; p.g=9.8; L0=10; omega0=sqrt(p.g/L0); T0=2*pi/omega0; c=4;
R0=L0/T0; fig.Renderer='Painters'; figure(1); clf

% Simulation parameters,
h=0.005; tol=0.01;

% Control parameters
fillets=false;  ls='k-';

for run=1:8
    switch run
        case 1, DeltaL_over_L0=0.01, R_over_R0=20
        case 2, DeltaL_over_L0=0.01, R_over_R0=10
        case 3, DeltaL_over_L0=0.01, R_over_R0=5
        case 4, DeltaL_over_L0=0.01, R_over_R0=2
        case 5, DeltaL_over_L0=0.01, R_over_R0=1
        case 6, DeltaL_over_L0=0.01, R_over_R0=.5
        case 7, DeltaL_over_L0=0.01, R_over_R0=.2
    end
    DeltaL=DeltaL_over_L0*L0; R=R_over_R0*R0; Ru=R; Rl=R;
    % initial conditions
    x=[0; 1; L0]; t=0; % phi=x(1), phidot=x(2), L=x(3)
    control=1; T=1000; t_on=0; t_off=0; clf;
    
    % Initialize saved variables (*.s), useful for plotting
    t_s=t; x_s=x; clear psi_s;
    
    for k=1:T/h
        psi_s(k)=x(1)*x(2);
        if control ~= 1
            L_new=L0;
        else
            if k==1
                L_new=min(x(3)+Ru*h,L0+DeltaL);
            elseif x(1)*x_s(1,k-1)<=0            % Algorithm 1
                delta=Rl/(x(2)^2);               % with or w/o fillets
                L_new=min(x(3)+Ru*h,L0+DeltaL);
                if fillets, [b,t1,t2,a]=fillet(Ru,x(3)-Ru*t,0,L0+DeltaL,-c); end
            elseif x(2)*x_s(2,k-1)<=0
                delta=abs(Ru/(x(1)*f(2)));
                L_new=max(x(3)-Rl*h,L0-DeltaL);
                if fillets, [b,t1,t2,a]=fillet(-Rl,x(3)+Rl*t,0,L0-DeltaL,c); end
            else
                if psi_s(k)>0
                    if psi_s(k)>=psi_s(k-1)
                        if fillets
                            if t<t1
                                L_new=x(3)+Ru*h;
                            elseif t>t2
                                L_new=L0+DeltaL;
                            else
                                L_new=-c*t^2/2+b*t+a;
                            end
                        else
                            L_new=min(x(3)+Ru*h,L0+DeltaL);
                        end
                    else
                        if fillets
                            s=delta*(x(2)^2+x(1)*f(2));
                            [b,t1,t2,a]=fillet(0,L0+DeltaL,s,L0+delta*psi_s(k)-s*t,-c);
                            if t<t1
                                L_new=L0+DeltaL;
                            elseif t>t2
                                L_new=L0+delta*psi_s(k);
                            else
                                L_new=-c*t^2/2+b*t+a;
                            end
                        else
                            L_new=min(L0+DeltaL,L0+delta*psi_s(k));
                        end
                    end
                else
                    if psi_s(k)<=psi_s(k-1)
                        if fillets
                            if t<t1
                                L_new=x(3)-Rl*h;
                            elseif t>t2
                                L_new=L0-DeltaL;
                            else
                                L_new=c*t^2/2+b*t+a;
                            end
                        else
                            L_new=max(x(3)-Rl*h,L0-DeltaL);
                        end
                    else
                        if fillets
                            s=delta*(x(2)^2+x(1)*f(2));
                            [b,t1,t2,a]=fillet(0,L0-DeltaL,s,L0+delta*psi_s(k)-s*t,c);
                            if t<t1
                                L_new=L0-DeltaL;
                            elseif t>t2
                                L_new=L0+delta*psi_s(k);
                            else
                                L_new=c*t^2/2+b*t+a;
                            end
                        else
                            L_new=max(L0-DeltaL,L0+delta*psi_s(k));
                        end
                    end
                end
            end
            if k>10 & x(2)*x_s(2,end-1)<=0 & abs(x(1))<tol  % if done, set L=L0
                L_new=L0; control=2; t_off=t; T=t_off+1*T0;
            end
        end
        p.Ldot=(L_new-x(3))/h;
        f1=RHS(x,p); f2=RHS(x+h*f1/2,p); f3=RHS(x+h*f2/2,p); f4=RHS(x+h*f3,p); % RK4
        f=f1/6+(f2+f3)/3+f4/6; x=x+h*f; t=t+h; x_s=[x_s x]; t_s=[t_s t];
        if t>T, break, end
    end
    psi_s=[psi_s x(1)*x(2)]; t_s=t_s/T0; x_s(3,:)=x_s(3,:)/L0;
    
    yyaxis left; Tm=5;
    plot(t_s,x_s(1,:),'b-'); hold on;
    % plot(t_s,x_s(2,:),'b--');
    axis([0 Tm -1.1 1.1]); plot([0 t_s(end)],[0 0],'b-');
    plot(t_s,psi_s,'b--');
    ylabel('h'); xlabel('t')
    yyaxis right; plot(t_s,x_s(3,:),ls); ylabel('L');
    axis([0 Tm 1-1.1*DeltaL/L0 1+1.1*DeltaL/L0]); hold on;
    if t_off>=1000, t_off=t_on; end
    T1_over_T0=(t_off-t_on)/T0
    disp(' '); pause
end
% end main script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(x,p)
phi=x(1); phidot=x(2); L=x(3);
% RHS of ODE
f=[phidot; -(p.g*sin(phi)+2*p.Ldot*phidot)/L; p.Ldot];
end % function RHS

function [b,x1,x2,a]=fillet(m1,b1,m2,b2,c)
b =(m1^2-m2^2+2*(b1-b2)*c)/(2*m1-2*m2);
x1=(m1-b)/c; x2=(m2-b)/c; a=m1*x1+b1-c*x1^2/2-b*x1;
% res=(m2*x2+b2 -c*x2^2/2 - b*x2 -a);
end % function fillet

