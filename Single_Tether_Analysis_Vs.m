disp('We first calculate, by symbolic differentiation and simplification, the (long)')
disp('expression for the symbolic function Ldot(t), which we will use later.')
clear; syms t phi(t) phidot(t) m g L0 C
Vs0(t) = 0.5*m*L0^2*phidot(t)^2 + m*g*L0*(1-cos(phi(t)));
L(t)  = L0*(1+phi(t)*phidot(t)*C/Vs0(t));                    % This is the control rule
Ldot(t) = simplify(diff(L(t),t))

disp('We have already inserted this expression for Ldot(t) in the analysis that follows, after:')
disp(' (a) manually stripping out the (t) notation [indicating Matlab "symbolic functions"],')
disp(' (b) replacing diff(phi,t) with phidot, and')
disp(' (c) replacing both diff(phidot,t) and diff(phi,t,t) with phidotdot.')
disp('We complete the rest of the analysis below using symbolic variables [without (t)],')
disp('instead of using symbolic functions [with (t)].'), disp(' ')

disp('We first rewrite the equation for L(phi,phidot), leveraging Vs [that is, a simplified V function]:')
clear, syms phi phidot L Ldot phidotdot m g L0 C
Vs0 = 0.5*m*L0^2*phidot^2 + m*g*L0*(1-cos(phi))
L   = L0*(1+phi*phidot*C/Vs0), disp(' ')

disp('Ldot(phi,phidot) and phidotdot(phi,phidot) are now determined by solving (symbolically)')
disp('a set of two simultaneous equations:')
eqn1 = Ldot == L0*((2*C*phi*phidotdot)/(L0^2*m*phidot^2 - 2*L0*g*m*(cos(phi) - 1)) + (2*C*phidot*phidot)/(L0^2*m*phidot^2 - 2*L0*g*m*(cos(phi) - 1)) - (4*C*L0*m*phi*phidot*(L0*phidot*phidotdot + g*sin(phi)*phidot))/(L0^2*m*phidot^2 - 2*L0*g*m*(cos(phi) - 1))^2)
eqn2 = phidotdot == -(g*sin(phi)+2*Ldot*phidot)/L;
SOL=solve(eqn1,eqn2,Ldot,phidotdot); Ldot=simplify(SOL.Ldot), phidotdot=simplify(SOL.phidotdot), disp(' ')

disp('Ldotdot(phi,phidot), Vs(phi,phidot) and Vdot(phi,phidot) may now be written, leveraging the above expressions:')
Vs    = simplify(0.5*m*(L^2*phidot^2) + m*g*L*(1-cos(phi)))
Vsdot = simplify(m*(L*Ldot*phidot^2+L^2*phidot*phidotdot)  + m*g*(Ldot*(1-cos(phi)) + L*phidot*sin(phi)))
disp(' ')

disp('We then assign {g,m,L0,C}={9.8,1,1,0.05}, perform some calculations over a range of {phi,phidot}, and plot.')
disp('[Give it some time; for Nphi=Nphidot=21, this calculation takes a few minutes.]')
g=9.8; m=1; L0=1; C=0.2; 
Lphi=1; Lphidot=2; N=40; eps=1e-9; Deltaphi=2*Lphi/(N-1); Deltaphidot=2*Lphidot/(N-1);

[phi_grid,phidot_grid] = meshgrid(-Lphi:Deltaphi:Lphi,-Lphidot:Deltaphidot:Lphidot);
for i=1:N; for j=1:N;
   phi=phi_grid(i,j); phidot=phidot_grid(i,j);
   phidotdot_grid(i,j)=eval(phidotdot);
   Vs_grid(i,j)=eval(Vs);
   Vsdot_grid(i,j)=eval(Vsdot);
end, end

figure(1); clf; surf(phi_grid,phidot_grid,phidotdot_grid), title('phidotdot(phi,phidot)'),  xlabel('phi'), ylabel('phidot')
figure(2); clf; surf(phi_grid,phidot_grid,Vs_grid),        title('Vs(phi,phidot)'),         xlabel('phi'), ylabel('phidot')
figure(3); clf; surf(phi_grid,phidot_grid,Vsdot_grid),     title('Vsdot(phi,phidot)'),      xlabel('phi'), ylabel('phidot')
Vsdot_grid
Vsdot_max_over_this_range=max(max(Vsdot_grid))
