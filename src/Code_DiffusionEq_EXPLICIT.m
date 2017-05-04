%MECE 5397 PROJECT
%Oscar Rodriguez
%10717156   
%
%Project B- Diffusion Equation
%Bc2-4

%Solving Equation Using Explicit Method

%Defining Parameters

clear all
clc

D=1;            %Given difrusivity       
ax=-3.14;       %Lower x Limit
ay=-3.14;       %Lower y Limit  
bx=3.14;        %Upper x limit
by=3.14;        %Upper y Limit

Nxx=10;                           %Number of steps in space(x)
Nyy=10;                           %Number of steps in space(y)       


nt=100;                           %Number of time steps 
delta_t=.1;                      %size time step
h1=6.28/(Nxx-1);                 %Size of space step(x)

x=ax:h1:bx;                        %Range of x(0,2) and specifying the grid points
y=ay:h1:by;                        %Range of y(0,2) and specifying the grid points

lamnda=((2*D*delta_t)/(h1^2));

if lamnda<=.5
    
else
    disp('Given parameters will make calculation blow up')
    return
end

                    
fa=y.*(y-ay).^2;
ga=((y-ay).^2).*cos(3.14.*y/(ay));

UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/n=UnN)

u=zeros(Nxx,Nyy);                  %Preallocating u
un=zeros(Nxx,Nyy);                 %Preallocating un
for it=0:nt
    un=u;
    h=surf(x,y,u);       %plotting the field variable
    shading interp
    axis ([-3.14 3.14 -3.14 3.14 -3.14 3.14])
    title({['2-D Diffusion with {\nu} = ',num2str(D)];['time (\itt) = ',num2str(it*delta_t)]})
    xlabel('(x)')
    ylabel('(y)')
    zlabel('Temperature')
    drawnow; 
    refreshdata(h)
       %Explicit method:
   for j=2:Nyy-1 
    for i=2:Nxx-1
    u(i,j)=un(i,j)+(D*delta_t*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(h1*h1))+(D*delta_t*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(h1*h1));
    %Boundary conditions
    %Dirichlet:
    u(1,:)=ga;
    u(Nxx,:)=fa;
    u(:,Nyy)=fa(10)+((x+3.14)/(6.28)*(ga(10)-fa(10)));
    %Neumann:
    
   u(:,1)=u(:,2)-UnS*h1;
    
   end
   end
   
end