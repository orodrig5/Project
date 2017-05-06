%implicit testing
clear all
clc
D=1;            %Given difrusivity       
ax=-3.14;       %Lower x Limit
ay=-3.14;       %Lower y Limit  
bx=3.14;        %Upper x limit
by=3.14;        %Upper y Limit

Nxx=39;                          %Number of steps in space(x)
Nyy=39;                          %Number of steps in space(y)       

nt=0.5;                         %Number of time steps 
delta_t=.065;                    %size time step
h1=6.28/(Nxx);                  %Size of space step(x)

u=zeros(Nxx+1,Nyy+1);          %Creating the grid for temperature        
un=zeros(Nxx+1,Nyy+1);         %initioal grid for temperature, Initial time condition is also represented by this line
x=ax:h1:bx;                      %defining the range of x
y=ay:h1:by;                      %defining the range of y



fa=y.*(y-ay).^2 ;                   %defining f(a)
ga=((y-ay).^2).*cos(3.14.*y/(ay)); %defining g(a)
Neuman_ay=0;               %defning the neumand condition given by the problem

u(1,:)=ga; %detriclet non linear  condition in the  right side boundary
u(Nxx,:)=fa; %detriclet condition non linear condition left side boudary 
u(:,Nyy)=fa(39)+((x+3.14)/(6.28)*(ga(39)-fa(39))); %detrichlet condition upper boundary 

lambda=((2*D*delta_t)/(h1^2))
a1(2:Nxx)=1;
b1(1:Nxx+1)=-2-2/lambda;
c1(1:Nxx+1)=1;
d1(1:Nxx+1)=0;
a1(1)=0;b1(1)=1;c1(1)=0;
a1(Nxx+1)=0;b1(Nxx+1)=1;c1(Nxx+1)=0;
a2(1:Nyy+1)=1;
b2(1:Nyy+1)=-2-2/lambda;
c2(1:Nyy+1)=1;  
a2(1)=0;b2(1)=1;c2(1)=0;
a2(Nyy+1)=0;b2(Nyy+1)=1;c2(Nyy+1)=0;
d2(1:Nyy+1)=0;
uhalf=zeros(Nxx+1,Nyy+1);
t=0;

for it=0:nt
    %un=u;                  %definite initial temperature (initial time condition)
    hei=surf(x,y,u');      %Initiation of 3 dimensional graph capable of showing time dependece of problem
    shading interp% color gradient for graph
    xlabel('(x)') %lable axixs
    ylabel('(y)') %lable axis
    zlabel('Temperature') %lable axis
   
    title({['2-D Diffusion with D = ',num2str(D)];['time (\itt) = ',num2str(it*delta_t)]})
      
    drawnow;
    refreshdata(hei) %frefresh of graph per time interval
for j=2:Nyy
% call TDMA
d1(1)=u(1,j);
d1(Nxx+1)=u(Nxx,j);
for i=2:Nxx
d1(i)=-1*u(i,j-1)-1*u(i,j+1)+(2*1-2/lambda)*u(i,j);
end; 
uhalf(:,j)=tridiag(Nxx+1,a1,b1,c1,d1);
end;
for i=2:Nxx
    % call TDMA
d2(1)=u(i,1);
d2(Nyy+1)=u(i,Nyy+1);
for j=2:Nyy
d2(j)=-uhalf(i-1,j)-uhalf(i+1,j)+(2-2/lambda)*uhalf(i,j);
end; 
u(:,1)=u(:,2)-Neuman_ay*h1; %Neuman condition bottom boundary
u(i,:)=tridiag(Nyy+1,a2,b2,c2,d2);
end;
end
