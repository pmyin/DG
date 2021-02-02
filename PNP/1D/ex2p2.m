function [sol, tall, energyall]=ex2p2(h,beta0, beta1,ts)
%k=2 in P_k
%legendre P: L0=1 L1=xi L2=(3xi^2-1)/2;  xi=2(x-x_j)/h

a=0;
b=1;


x=(a:h:(b-h))'+h/2;
N=length(x);
%beta0=5;beta1=5;
t=0; dt=0.01*h^2/beta0;
sigmaa=0;
sigmab=0;
q1=1;
q2=-1;

tall=[];
energyall=[];

xi3=sqrt(3/7-2*sqrt(6/5)/7); omega3=0.5+sqrt(30)/36;
xi2=-sqrt(3/7-2*sqrt(6/5)/7); omega2=0.5+sqrt(30)/36;
xi4=sqrt(3/7+2*sqrt(6/5)/7); omega4=0.5-sqrt(30)/36;
xi1=-sqrt(3/7+2*sqrt(6/5)/7); omega1=0.5-sqrt(30)/36;




xi=[xi1,xi2, xi3, xi4]; omega=[omega1 omega2 omega3 omega4];
Lxi=[1 1 1 1;xi1 xi2 xi3 xi4 ; 0.5*(3*xi1^2-1) 0.5*(3*xi2^2-1) 0.5*(3*xi3^2-1) 0.5*(3*xi4^2-1)]';
Lxi1=[1;xi1;0.5*(3*xi1^2-1)];
Lxi2=[1;xi2;0.5*(3*xi2^2-1)];
Lxi3=[1;xi3;0.5*(3*xi3^2-1)];
Lxi4=[1;xi4;0.5*(3*xi4^2-1)];
% Legenger polynomial evaluated at xi1--xi4 (time independent)
L_xixi1=[0;1;3*xi1];
L_xixi2=[0;1;3*xi2];
L_xixi3=[0;1;3*xi3];
L_xixi4=[0;1;3*xi4];    %L_xi(xi1) etc.
%define L(1), L(-1), L_xi(1), L_xi(-1),L_{xi xi}(1), L_{xi xi}(-1),D, E
L1=[1;1;1];  LN1=[1;-1;1];
L_xi1=[0;1;3]; L_xiN1=[0;1;-3];
L_xi_xi1=[0; 0; 3];L_xi_xiN1=[0 ;0 ;3];
D=beta0*L1-L_xi1+4*beta1*L_xi_xi1;
E=beta0*LN1+L_xiN1+4*beta1*L_xi_xiN1;

Dp=beta0*L1-L_xi1;
Ep=beta0*LN1+L_xiN1;

M=zeros(3*N,3*N);c1rhs=zeros(3*N,1); c2rhs=zeros(3*N,1);
MM=[1 0 0; 0 1/3 0; 0 0 1/5]*h;
for j=1:N
    M(3*(j-1)+1:3*j,3*(j-1)+1:3*j)=[1 0 0; 0 1/3 0; 0 0 1/5]*h;
    xtemp=h*(xi)/2+x(j);
    c1=1+pi*sin(pi*xtemp);
    c2=4-2*xtemp;
    
    c1rhs(3*(j-1)+1:3*j)=(c1.*omega *Lxi)';
    c2rhs(3*(j-1)+1:3*j)=(c2.*omega *Lxi)';
   
    %Crhs(3*(j-1)+1:3*j) =( C.*omega *Lxi)';
    
end
c1rhs=c1rhs*h/2; c2rhs=c2rhs*h/2;  




c1=M\c1rhs; c2=M\c2rhs; 


G=q1*c1+q2*c2;


% 
%             figure(1), clf; hold on
%             subplot(1,3,1),hold on
%             for i=1:length(x)
%                 xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                 xitemp=(xtemp-x(i))*2/h;
%                 ytemp=c1(3*(i-1)+1)+ c1(3*(i-1)+2)*xitemp +c1(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 plot(xtemp,ytemp,'r-')
%             end
%               plot(x,c1(1:3:end),'b-*')
%             title(['c1 at t=',num2str(t)])
%             subplot(1,3,2),hold on
%             for i=1:length(x)
%                 xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                 xitemp=(xtemp-x(i))*2/h;
%                 ytemp=c2(3*(i-1)+1)+ c2(3*(i-1)+2)*xitemp +c2(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 plot(xtemp,ytemp,'r-')
%             end
%             plot(x,c2(1:3:end),'b-*')
%             title(['c2 at t=',num2str(t)])
%             pause








% construct P * \psi = b

A=-LN1*Dp'-L_xiN1*L1';
B=2*[0 0 0; 0 2 0; 0 0 6]+LN1*Ep'+ L1*Dp'+L_xiN1*LN1'-L_xi1*L1';
C=-L1*Ep'+L_xi1*LN1';

P=zeros(3*N,3*N);b=zeros(3*N,1);
for j=2:N-1
    P(3*(j-1)+1:3*(j-1)+3,3*(j-2)+1:3*(j-2)+3)=A;
    P(3*(j-1)+1:3*(j-1)+3,3*(j-1)+1:3*(j-1)+3)=B;
    P(3*(j-1)+1:3*(j-1)+3,3*(j  )+1:3*(j  )+3)=C;
    b(3*(j-1)+1:3*(j-1)+3,1)=h*MM*G(3*j-2:3*j);
end



%boundary conditions for potential
j=1; %V_a=-1 V_b=1

%set psi(0)=0
P(3*(j-1)+1:3*(j-1)+3,3*(j-1)+1:3*(j-1)+3)=B;
P(3*(j-1)+1:3*(j-1)+3,3*(j-0)+1:3*(j-0)+3)=C;
b(3*(j-1)+1:3*(j-1)+3,1)= h*MM*G(3*j-2:3*j)-LN1*h*sigmaa/2;

%b(3*(j-1)+1:3*(j-1)+3,1)= b(3*(j-1)+1:3*(j-1)+3,1) +
%(beta0*LN1+L_xiN1)*(-1);
%b(1,1)=-1;


j=N; 
P(3*(j-1)+1:3*(j-1)+3,3*(j-2)+1:3*(j-2)+3)=A;
P(3*(j-1)+1:3*(j-1)+3,3*(j-1)+1:3*(j-1)+3)=B-beta0*(L1*L1')+L_xi1*L1';
b(3*(j-1)+1:3*(j-1)+3,1)=h*MM*G(3*j-2:3*j) +L1*h*sigmab/2; 
%b(3*(j-1)+1:3*(j-1)+3,1)=b(3*(j-1)+1:3*(j-1)+3,1)- (-beta0*L1+L_xi1)*(1); 
%b(3*N,1)=1;
%P(3*(j-1)+1:3*(j-1)+3,3*(j-0)+1:3*(j-0)+3)=C;

%solving for V= psi

det(P)

pause

psi=P\b;

pause

%             K=3;
%          figure(1), clf; hold on
%                 subplot(1,3,1); hold on
%                 xx=[];
%                 for i=1:length(x)
%                     xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                     xx=[xx,xtemp];
%                     xitemp=(xtemp-x(i))*2/h;
%                     ytemp=c1(K*(i-1)+1)+ c1(K*(i-1)+2)*xitemp +c1(K*(i-1)+3)*(3*xitemp.^2-1)/2;
%                     plot(xtemp,ytemp,'r:o')
%                 end
%     
%                 c1exact=xx.^2.*(1-xx).^2*exp(-t);
%                   plot(xx,c1exact,'b-')
%                 title(['c1 at t=',num2str(t)])
%                 subplot(1,3,2),hold on
%                 for i=1:length(x)
%                     xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                     xitemp=(xtemp-x(i))*2/h;
%                     ytemp=c2(K*(i-1)+1)+ c2(K*(i-1)+2)*xitemp +c2(K*(i-1)+3)*(3*xitemp.^2-1)/2;
%                     plot(xtemp,ytemp,'r-o')
%                 end
%                 c2exact=xx.^2.*(1-xx).^3*exp(-t);
%                 plot(xx,c2exact,'b-')
%                 title(['c2 at t=',num2str(t)])
%                 subplot(1,3,3),hold on
%                 for i=1:length(x)
%                     xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                     xitemp=(xtemp-x(i))*2/h;
%                     ytemp=psi(K*(i-1)+1)+ psi(K*(i-1)+2)*xitemp +psi(K*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 
%                     plot(xtemp,ytemp,'r-o')
%                 end
%                 psiexact=-(10*xx.^7-28*xx.^6+21*xx.^5)*exp(-t)/420;
%                 plot(xx,psiexact,'b-',x,psi(1:K:end),'k:*')
%                 title(['psi at t=',num2str(t)])
%                 pause






% solving for g in c=gM, M=exp(-q\psi)
g1=zeros(3*N,1);g2=zeros(3*N,1);


for j=1:N
    % 4-point Gauss-Legendre quadrature
    cj1=c1(3*(j-1)+1:3*(j-1)+3); cj2=c2(3*(j-1)+1:3*(j-1)+3);
    psij=psi(3*(j-1)+1:3*(j-1)+3);
    % define matrix in front of g in RHS
    
     K1=omega1*(Lxi1*Lxi1')*exp(-q1*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q1*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q1*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q1*Lxi4'*psij);
      K1=K1*h/2;
     K2=omega1*(Lxi1*Lxi1')*exp(-q2*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q2*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q2*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q2*Lxi4'*psij);
            K2=K2*h/2;
    

    g1(3*(j-1)+1:3*(j-1)+3)=K1\(MM*cj1);
    g2(3*(j-1)+1:3*(j-1)+3)=K2\(MM*cj2);
end

%make g positive
for j=1:N
    g1j=g1(3*(j-1)+1:3*(j-1)+3);
    g2j=g2(3*(j-1)+1:3*(j-1)+3);
    psij=psi(3*(j-1)+1:3*(j-1)+3);
%     uval1=Lxi*g1j; uval2=Lxi*g2j;
    psijval=Lxi*psij;
    % ubar1=omega*uval1/2;     ubar2=omega*uval2/2;
    a1j= sum(omega'.*(exp(-q1*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1-xi'));
    a2j= sum(omega'.*(exp(-q2*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1-xi'));
    b1j= sum(omega'.*(exp(-q1*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1+xi'));
    b2j= sum(omega'.*(exp(-q2*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1+xi'));
    gamma1=(a1j+b1j)/2; gamma2=(a2j+b2j)/2;
    g1gamma1val=g1j'*[1;gamma1;(3*gamma1^2-1)/2];
    g2gamma2val=g2j'*[1;gamma2;(3*gamma2^2-1)/2];
        
    umin = min([g1gamma1val;L1'*g1j;LN1'*g1j]);
    if (umin < 0)
        
%         gbartt=sum(omega'.*(exp(-q1*psijval).*uval1))/sum(omega'.*exp(-q1*psijval));
         gbar=2*c1(3*j-2)/sum(omega'.*exp(-q1*psijval));
%         [gbartt, gbar]
%         if abs(gbartt-gbar)>10^-3
%            error('gbartt not the same as gbar\n');
%         end
           
        theta=min(1,gbar/(gbar-umin));
        g1(3*j-2:3*j) = theta*(g1(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        %fprintf('c1:j=%d\n',j);
    end
    umin = min([g2gamma2val;L1'*g2j;LN1'*g2j]);
    
    if (umin < 0)
        gbar=2*c2(3*j-2)/sum(omega'.*exp(-q2*psijval));
        theta=min(1,gbar/(gbar-umin));
        g2(3*j-2:3*j) = theta*(g2(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        % fprintf('c2:j=%d\n',j);
    end
    
end


%

%starting evolution 2nd oder R-K method (Heun's to be exact)

% constructing matrix T  to solve for c in tme, i.e., T * (c)_t = b
% note that T is time independent, so it is predefined before evolution




%sol=[c(1:3:3*N) c(2:3:3*N) c(3:3:3*N) psi(1:3:3*N) psi(2:3:3*N) psi(3:3:3*N)  q(1:3:3*N) q(2:3:3*N) q(3:3:3*N) ]
%pause

LL=1;
%save energyintimet00

while t<ts

    if t+dt>ts
        dt=ts-t;
    end
    
    % constructing rhs in T * c_t =rhs
    c1R1=zeros(3*N,1); c1R2p=c1R1; c1R2n=c1R1; c1R3p=c1R1; c1R3n=c1R1;
    c2R1=zeros(3*N,1); c2R2p=c2R1; c2R2n=c2R1; c2R3p=c2R1; c2R3n=c2R1;
    for j=1:N
        % 4-point Gauss-Legendre quadrature for R1
         g1j=g1(3*(j-1)+1:3*(j-1)+3);
         g2j=g2(3*(j-1)+1:3*(j-1)+3);
        psij=psi(3*(j-1)+1:3*(j-1)+3);
        
        c1R1(3*(j-1)+1:3*(j-1)+3)=-omega1*exp(-q1*Lxi1'*psij)*L_xixi1'*g1j*L_xixi1...
                                  -omega2*exp(-q1*Lxi2'*psij)*L_xixi2'*g1j*L_xixi2...
                                  -omega3*exp(-q1*Lxi3'*psij)*L_xixi3'*g1j*L_xixi3...
                                  -omega4*exp(-q1*Lxi4'*psij)*L_xixi4'*g1j*L_xixi4;
                              
        c2R1(3*(j-1)+1:3*(j-1)+3)=-omega1*exp(-q2*Lxi1'*psij)*L_xixi1'*g2j*L_xixi1...
                                  -omega2*exp(-q2*Lxi2'*psij)*L_xixi2'*g2j*L_xixi2...
                                  -omega3*exp(-q2*Lxi3'*psij)*L_xixi3'*g2j*L_xixi3...
                                  -omega4*exp(-q2*Lxi4'*psij)*L_xixi4'*g2j*L_xixi4;
                              
         %R2 and R3
        if j<N
               g1jp1=g1(3*(j)+1:3*(j)+3);
               g2jp1=g2(3*(j)+1:3*(j)+3);
              psijp1=psi(3*(j)+1:3*(j)+3);
            c1R2p(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psij)+exp(-q1*LN1'*psijp1))*(-D'*g1j+E'*g1jp1)*L1;
            c1R3p(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psij)+exp(-q1*LN1'*psijp1))*(L1'*g1j-LN1'*g1jp1)*L_xi1;
            c2R2p(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psij)+exp(-q2*LN1'*psijp1))*(-D'*g2j+E'*g2jp1)*L1;
            c2R3p(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psij)+exp(-q2*LN1'*psijp1))*(L1'*g2j-LN1'*g2jp1)*L_xi1;
            
%         else %j=N
%             q1squarebracket=1-L1'*q1j;
%             q2squarebracket=1-L1'*q2j;
%             c1R2p(3*(j-1)+1:3*(j-1)+3)=(1+L1'*c1j)*(beta0*q1squarebracket + L_xi1'*q1j*2)*L1;
%             c1R3p(3*(j-1)+1:3*(j-1)+3)=(1+L1'*c1j)*(-1)* q1squarebracket*L_xi1;
%             c2R2p(3*(j-1)+1:3*(j-1)+3)=(0+L1'*c2j)*(beta0*q2squarebracket + L_xi1'*q2j*2)*L1;
%             c2R3p(3*(j-1)+1:3*(j-1)+3)=(0+L1'*c1j)*(-1)* q2squarebracket*L_xi1;
        end
        
        if j>1
            
             g1jn1=g1(3*(j-2)+1:3*(j-2)+3);
             g2jn1=g2(3*(j-2)+1:3*(j-2)+3);
            psijn1=psi(3*(j-2)+1:3*(j-2)+3);
            c1R2n(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psijn1)+exp(-q1*LN1'*psij))*(-D'*g1jn1+E'*g1j)*LN1;
            c1R3n(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psijn1)+exp(-q1*LN1'*psij))*(L1'*g1jn1-LN1'*g1j)*L_xiN1;
            c2R2n(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psijn1)+exp(-q2*LN1'*psij))*(-D'*g2jn1+E'*g2j)*LN1;
            c2R3n(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psijn1)+exp(-q2*LN1'*psij))*(L1'*g2jn1-LN1'*g2j)*L_xiN1;
%         else %j=1
%             c1R2n(3*(j-1)+1:3*(j-1)+3)=(0+LN1'*c1j)*(beta0*q1squarebracket + L_xiN1'*q1j*2)*LN1;
%             c1R3n(3*(j-1)+1:3*(j-1)+3)=-(0+LN1'*c2j)*q1squarebracket*L_xiN1;
%             c2R2n(3*(j-1)+1:3*(j-1)+3)= (1+LN1'*c2j)*(beta0*q2squarebracket + L_xiN1'*q2j*2)*LN1;
%             c2R3n(3*(j-1)+1:3*(j-1)+3)=-(1+LN1'*c2j)*q2squarebracket*L_xiN1;
        end
    end
    
    c1R2=c1R2p-c1R2n; c2R2=c2R2p-c2R2n;
    c1R3=c1R3p+c1R3n; c2R3=c2R3p+c2R3n;
    
    
    
    b1temp=M\(2*c1R1/h+(c1R2+c1R3)/(2*h));
    b2temp=M\(2*c2R1/h+(c2R2+c2R3)/(2*h));
    c1temp=c1+dt*b1temp;
    c2temp=c2+dt*b2temp;
  
    G=q1*c1temp+q2*c2temp;
    
    for j=1:N
    btemp(3*(j-1)+1:3*(j-1)+3,1)=h*MM*G(3*j-2:3*j);
    end
    j=1;
    btemp(3*(j-1)+1:3*(j-1)+3,1)= btemp(3*(j-1)+1:3*(j-1)+3,1)-LN1*h*sigmaa/2;
   
    j=N;
    btemp(3*(j-1)+1:3*(j-1)+3,1)= btemp(3*(j-1)+1:3*(j-1)+3,1) +L1*h*sigmab/2;   
 
    
    
    %btemp(1:3,1)=btemp(1:3,1)-L_xiN1/h;
    %btemp(3*N-2:3*N,1)=btemp(3*N-2:3*N,1)-L_xi1/h;
    psitemp= P\btemp;
    
    
    % solving for gtemp
g1temp=zeros(3*N,1);g2temp=zeros(3*N,1);
  
for j=1:N
    % 4-point Gauss-Legendre quadrature
    cj1=c1temp(3*(j-1)+1:3*(j-1)+3); cj2=c2temp(3*(j-1)+1:3*(j-1)+3);
    psij=psitemp(3*(j-1)+1:3*(j-1)+3);
    % define matrix in front of g in RHS
    
     K1=omega1*(Lxi1*Lxi1')*exp(-q1*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q1*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q1*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q1*Lxi4'*psij);
      K1=K1*h/2;
     K2=omega1*(Lxi1*Lxi1')*exp(-q2*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q2*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q2*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q2*Lxi4'*psij);
            K2=K2*h/2;
    

    g1temp(3*(j-1)+1:3*(j-1)+3)=K1\(MM*cj1);
    g2temp(3*(j-1)+1:3*(j-1)+3)=K2\(MM*cj2);
end
    %make gtemp positive
for j=1:N
    g1j=g1temp(3*(j-1)+1:3*(j-1)+3);
    g2j=g2temp(3*(j-1)+1:3*(j-1)+3);
    psij=psitemp(3*(j-1)+1:3*(j-1)+3);
    %uval1=Lxi*g1j; uval2=Lxi*g2j;
    psijval=Lxi*psij;
    % ubar1=omega*uval1/2;     ubar2=omega*uval2/2;
    a1j= sum(omega'.*(exp(-q1*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1-xi'));
    a2j= sum(omega'.*(exp(-q2*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1-xi'));
    b1j= sum(omega'.*(exp(-q1*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1+xi'));
    b2j= sum(omega'.*(exp(-q2*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1+xi'));
    gamma1=(a1j+b1j)/2; gamma2=(a2j+b2j)/2;
    g1gamma1val=g1j'*[1;gamma1;(3*gamma1^2-1)/2];
    g2gamma2val=g2j'*[1;gamma2;(3*gamma2^2-1)/2];
        
    umin = min([g1gamma1val;L1'*g1j;LN1'*g1j]);
    if (umin < 0)
        
        %gbar=sum(omega'.*(exp(-q1*psijval).*uval1))/sum(omega'.*exp(-q1*psijval));
        gbar=2*c1temp(3*j-2)/sum(omega'.*exp(-q1*psijval));
        theta=min(1,gbar/(gbar-umin));
        g1temp(3*j-2:3*j) = theta*(g1temp(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        %fprintf('c1:j=%d\n',j);
    end
    umin = min([g2gamma2val;L1'*g2j;LN1'*g2j]);
    
    if (umin < 0)
        gbar=2*c2temp(3*j-2)/sum(omega'.*exp(-q2*psijval));
        theta=min(1,gbar/(gbar-umin));
        g2temp(3*j-2:3*j) = theta*(g2temp(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        % fprintf('c2:j=%d\n',j);
    end
    
end

    
    
    %second step in RK2
    
    
    c1R1=zeros(3*N,1); c1R2p=c1R1; c1R2n=c1R1; c1R3p=c1R1; c1R3n=c1R1;
    c2R1=zeros(3*N,1); c2R2p=c2R1; c2R2n=c2R1; c2R3p=c2R1; c2R3n=c2R1;
   for j=1:N
        % 4-point Gauss-Legendre quadrature for R1
         g1j=g1temp(3*(j-1)+1:3*(j-1)+3);
         g2j=g2temp(3*(j-1)+1:3*(j-1)+3);
        psij=psitemp(3*(j-1)+1:3*(j-1)+3);
        
        c1R1(3*(j-1)+1:3*(j-1)+3)=-omega1*exp(-q1*Lxi1'*psij)*L_xixi1'*g1j*L_xixi1...
                                  -omega2*exp(-q1*Lxi2'*psij)*L_xixi2'*g1j*L_xixi2...
                                  -omega3*exp(-q1*Lxi3'*psij)*L_xixi3'*g1j*L_xixi3...
                                  -omega4*exp(-q1*Lxi4'*psij)*L_xixi4'*g1j*L_xixi4;
                              
        c2R1(3*(j-1)+1:3*(j-1)+3)=-omega1*exp(-q2*Lxi1'*psij)*L_xixi1'*g2j*L_xixi1...
                                  -omega2*exp(-q2*Lxi2'*psij)*L_xixi2'*g2j*L_xixi2...
                                  -omega3*exp(-q2*Lxi3'*psij)*L_xixi3'*g2j*L_xixi3...
                                  -omega4*exp(-q2*Lxi4'*psij)*L_xixi4'*g2j*L_xixi4;
                              
         %R2 and R3
        if j<N
               g1jp1=g1temp(3*(j)+1:3*(j)+3);
               g2jp1=g2temp(3*(j)+1:3*(j)+3);
              psijp1=psitemp(3*(j)+1:3*(j)+3);
            c1R2p(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psij)+exp(-q1*LN1'*psijp1))*(-D'*g1j+E'*g1jp1)*L1;
            c1R3p(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psij)+exp(-q1*LN1'*psijp1))*(L1'*g1j-LN1'*g1jp1)*L_xi1;
            c2R2p(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psij)+exp(-q2*LN1'*psijp1))*(-D'*g2j+E'*g2jp1)*L1;
            c2R3p(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psij)+exp(-q2*LN1'*psijp1))*(L1'*g2j-LN1'*g2jp1)*L_xi1;
            
%         else %j=N
%             q1squarebracket=1-L1'*q1j;
%             q2squarebracket=1-L1'*q2j;
%             c1R2p(3*(j-1)+1:3*(j-1)+3)=(1+L1'*c1j)*(beta0*q1squarebracket + L_xi1'*q1j*2)*L1;
%             c1R3p(3*(j-1)+1:3*(j-1)+3)=(1+L1'*c1j)*(-1)* q1squarebracket*L_xi1;
%             c2R2p(3*(j-1)+1:3*(j-1)+3)=(0+L1'*c2j)*(beta0*q2squarebracket + L_xi1'*q2j*2)*L1;
%             c2R3p(3*(j-1)+1:3*(j-1)+3)=(0+L1'*c1j)*(-1)* q2squarebracket*L_xi1;
        end
        
        if j>1
            
             g1jn1=g1temp(3*(j-2)+1:3*(j-2)+3);
             g2jn1=g2temp(3*(j-2)+1:3*(j-2)+3);
            psijn1=psitemp(3*(j-2)+1:3*(j-2)+3);
            c1R2n(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psijn1)+exp(-q1*LN1'*psij))*(-D'*g1jn1+E'*g1j)*LN1;
            c1R3n(3*(j-1)+1:3*(j-1)+3)=(exp(-q1*L1'*psijn1)+exp(-q1*LN1'*psij))*(L1'*g1jn1-LN1'*g1j)*L_xiN1;
            c2R2n(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psijn1)+exp(-q2*LN1'*psij))*(-D'*g2jn1+E'*g2j)*LN1;
            c2R3n(3*(j-1)+1:3*(j-1)+3)=(exp(-q2*L1'*psijn1)+exp(-q2*LN1'*psij))*(L1'*g2jn1-LN1'*g2j)*L_xiN1;
%         else %j=1
%             c1R2n(3*(j-1)+1:3*(j-1)+3)=(0+LN1'*c1j)*(beta0*q1squarebracket + L_xiN1'*q1j*2)*LN1;
%             c1R3n(3*(j-1)+1:3*(j-1)+3)=-(0+LN1'*c2j)*q1squarebracket*L_xiN1;
%             c2R2n(3*(j-1)+1:3*(j-1)+3)= (1+LN1'*c2j)*(beta0*q2squarebracket + L_xiN1'*q2j*2)*LN1;
%             c2R3n(3*(j-1)+1:3*(j-1)+3)=-(1+LN1'*c2j)*q2squarebracket*L_xiN1;
        end
    end
    
    c1R2=c1R2p-c1R2n; c2R2=c2R2p-c2R2n;
    c1R3=c1R3p+c1R3n; c2R3=c2R3p+c2R3n;
    
    

    
    b1=M\(2*c1R1/h+(c1R2+c1R3)/(2*h));
    b2=M\(2*c2R1/h+(c2R2+c2R3)/(2*h));
    c1=c1+0.5*dt*(b1+b1temp);
    c2=c2+0.5*dt*(b2+b2temp);
    
    
      G=q1*c1+q2*c2;
   
     for j=1:N
    b(3*(j-1)+1:3*(j-1)+3,1)=h*MM*G(3*j-2:3*j);
    end
    j=1;
    b(3*(j-1)+1:3*(j-1)+3,1)= b(3*(j-1)+1:3*(j-1)+3,1)-LN1*h*sigmaa/2;
   
    j=N;
    b(3*(j-1)+1:3*(j-1)+3,1)= b(3*(j-1)+1:3*(j-1)+3,1)+ L1*h*sigmab/2;   
 
    
    
    psi= P\b;
    
   % solving for q
g1=zeros(3*N,1);g2=zeros(3*N,1);

for j=1:N
    % 4-point Gauss-Legendre quadrature
    cj1=c1(3*(j-1)+1:3*(j-1)+3); cj2=c2(3*(j-1)+1:3*(j-1)+3);
    psij=psi(3*(j-1)+1:3*(j-1)+3);
    % define matrix in front of g in RHS
    
     K1=omega1*(Lxi1*Lxi1')*exp(-q1*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q1*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q1*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q1*Lxi4'*psij);
      K1=K1*h/2;
     K2=omega1*(Lxi1*Lxi1')*exp(-q2*Lxi1'*psij)...
        +omega2*(Lxi2*Lxi2')*exp(-q2*Lxi2'*psij)...
         +omega3*(Lxi3*Lxi3')*exp(-q2*Lxi3'*psij)...
          +omega4*(Lxi4*Lxi4')*exp(-q2*Lxi4'*psij);
            K2=K2*h/2;
    

    g1(3*(j-1)+1:3*(j-1)+3)=K1\(MM*cj1);
    g2(3*(j-1)+1:3*(j-1)+3)=K2\(MM*cj2);
end
    
        
%make g positive
  for j=1:N
    g1j=g1(3*(j-1)+1:3*(j-1)+3);
    g2j=g2(3*(j-1)+1:3*(j-1)+3);
    psij=psi(3*(j-1)+1:3*(j-1)+3);
    %uval1=Lxi*g1j; uval2=Lxi*g2j;
    psijval=Lxi*psij;
    % ubar1=omega*uval1/2;     ubar2=omega*uval2/2;
    a1j= sum(omega'.*(exp(-q1*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1-xi'));
    a2j= sum(omega'.*(exp(-q2*psijval).*(xi'-xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1-xi'));
    b1j= sum(omega'.*(exp(-q1*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q1*psijval).*(1+xi'));
    b2j= sum(omega'.*(exp(-q2*psijval).*(xi'+xi'.^2)))/sum(omega'.*exp(-q2*psijval).*(1+xi'));
    gamma1=(a1j+b1j)/2; gamma2=(a2j+b2j)/2;
    g1gamma1val=g1j'*[1;gamma1;(3*gamma1^2-1)/2];
    g2gamma2val=g2j'*[1;gamma2;(3*gamma2^2-1)/2];
        
    umin = min([g1gamma1val;L1'*g1j;LN1'*g1j]);
    if (umin < 0)
         
        %gbar=sum(omega'.*(exp(-q1*psijval).*uval1))/sum(omega'.*exp(-q1*psijval));
        gbar=2*c1(3*j-2)/sum(omega'.*exp(-q1*psijval));
        theta=min(1,gbar/(gbar-umin));
        g1(3*j-2:3*j) = theta*(g1(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        %fprintf('c1:j=%d\n',j);
    end
    umin = min([g2gamma2val;L1'*g2j;LN1'*g2j]);
    
    if (umin < 0)
        gbar=2*c2(3*j-2)/sum(omega'.*exp(-q2*psijval));
        theta=min(1,gbar/(gbar-umin));
        g2(3*j-2:3*j) = theta*(g2(3*j-2:3*j)-[gbar;0;0])+[gbar;0;0];
        % fprintf('c2:j=%d\n',j);
    end
    
   end
        t=t+dt;
    
    
%         if abs(t-0.01*LL)<dt
%             LL=LL+1;
%     
%             figure(1), clf; hold on
%             subplot(1,3,1); hold on
%             xx=[];
%             for i=1:length(x)
%                 xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%         
%                 xitemp=(xtemp-x(i))*2/h;
%                 ytemp=c1(3*(i-1)+1)+ c1(3*(i-1)+2)*xitemp +c1(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 plot(xtemp,ytemp,'r:o')
%             end
%            
% 
%             title(['c1 at t=',num2str(t)])
%             subplot(1,3,2),hold on
%             for i=1:length(x)
%                 xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                 xitemp=(xtemp-x(i))*2/h;
%                 ytemp=c2(3*(i-1)+1)+ c2(3*(i-1)+2)*xitemp +c2(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 plot(xtemp,ytemp,'r-o')
%             end
% 
%             title(['c2 at t=',num2str(t)])
%             subplot(1,3,3),hold on
%             for i=1:length(x)
%                 xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%                 xitemp=(xtemp-x(i))*2/h;
%                 ytemp=psi(3*(i-1)+1)+ psi(3*(i-1)+2)*xitemp +psi(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%                 plot(xtemp,ytemp,'r-o')
%             end
% 
%             title(['psi at t=',num2str(t)])
%             pause(0.001)
%     

%         end
      
    
    

    
    
%     energy=0;
%     for j=1:N
%         c1val=Lxi*c1(3*j-2:3*j);
%         c2val=Lxi*c2(3*j-2:3*j);
%         psival=Lxi*psi(3*j-2:3*j);
%         
%         energy=energy+omega*( c1val.^2+c2val.^2 + psival.*(c2val-c1val) );
%     end
%     energy=energy*h/2;
%     
%    
%     tall=[tall, t];
%     energyall = [energyall,energy];
%     
%     
    if abs(t-0.01)<=dt/2
        save energyintimet01
    end
    if abs(t-0.005)<=dt/2
        save energyintimet005
    end
    if abs(t-0.02)<=dt/2
        save energyintimet02
    end
    if abs(t-0.03)<=dt/2
        save energyintimet03
    end
    if abs(t-0.04)<=dt/2
        save energyintimet04
    end
    
    if abs(t-0.05)<=dt/2
        save energyintimet05
    end
    if abs(t-0.06)<=dt/2
        save energyintimet06
    end
    if abs(t-0.07)<=dt/2
        save energyintimet07
    end
    if abs(t-0.08)<=dt/2
        save energyintimet08
    end
    
    
    if abs(t-0.09)<=dt/2
        save energyintimet09
    end
    
    
%             figure(1), clf; hold on
%             for i=1:length(x)
%             xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%             xitemp=(xtemp-x(i))*2/h;
%             ytemp=c(3*(i-1)+1)+ c(3*(i-1)+2)*xitemp +c(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%             plot(xtemp,ytemp,'r-o')
%             end
%             title('c')
%         %     subplot(1,3,1), plot(x,c(1:3:3*N),'ro-'); title(['c at t=',num2str(t)]);%ylim([ymin,ymax])
%         %     subplot(1,3,2), plot(x,psi(1:3:3*N),'ro-'); title(['\psi at t=',num2str(t)]);%ylim([ymin,ymax])
%         %     subplot(1,3,3), plot(x,q(1:3:3*N),'ro-'); title(['q at t=',num2str(t)]);%ylim([ymin+log(ymin),ymax+log(ymax)])
%             figure(2), clf; hold on
%             for i=1:length(x)
%             xtemp=linspace(x(i)-h/2, x(i)+h/2, 11);
%             xitemp=(xtemp-x(i))*2/h;
%             ytemp=psi(3*(i-1)+1)+ psi(3*(i-1)+2)*xitemp +psi(3*(i-1)+3)*(3*xitemp.^2-1)/2;
%             plot(xtemp,ytemp,'r-o')
%             end
%             title('psi')
%             pause(0.001)
    
    
 energy=0;
    for j=1:N
        c1val=Lxi*c1(3*j-2:3*j);
        c2val=Lxi*c2(3*j-2:3*j);
        psival=Lxi*psi(3*j-2:3*j);
        
        energy=energy+omega*(  c1val.*log(c1val)+ c2val.*log(c2val)+ 0.5*psival.*(q1*c1val+q2*c2val) );
    end
    energy=energy*h/2;
    
   
    tall=[tall, t];
    energyall = [energyall,energy];
    
    
end
t

sol=[c1 c2 g1 g2 psi];


