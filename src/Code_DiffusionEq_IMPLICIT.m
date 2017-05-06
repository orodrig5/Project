%MECE 5397 PROJECT
%Oscar Rodriguez
%10717156   
%
%Project B- Diffusion Equation
%Bc2-4

%Solving Equation Using Explicit Method

%Defining Parameters

bx=6.28; by=6.28; lb=2; 
u_initial=0; 
D=1; %Deffusive 

Xx=60; Yy=60; %Mesh size
nx=Xx+1; ny=Yy+1; %Finding the correct size for iteration of the tridiagonals
dx=bx/Xx; dy=by/Yy; %Space steps
del=(dx*dx)/(dy*dy); %dummy variable used in tridiagonal calulation
dt=0.99; tn=20; %time step size and total time of calculation
lambda=D*dt/(dx*dx); %lamda needed for the coefficients of the tridiagonal matrix
x=linspace(0,bx,nx); y=linspace(0,by,ny); %size of x and y coordinates


u=zeros(nx,ny); u(:,:)=u_initial;

% boundary conditions
%left b.c.
fa=y.*(y).^2 ;                   %defining f(a)
ga=((y).^2).*cos(-y);            %defining g(a)

u(1,:)=ga; %left boundary condition

u(nx,:)=fa; % right boundary condition

% top boundary
u(:,ny)=fa(ny)+((x+3.14)/(6.28)*(ga(ny)-fa(ny))); %top boundary condition
Neuman_ay=0;%dummy for variable for the neuman boundary


surf(x,y,u'); shading interp;


%defining coefficents for the varilables in the tridiagonal matrix
g1(2:Xx)=1;
b1(1:Xx+1)=-2-2/lambda;
c1(1:Xx+1)=1;
f1(1:Xx+1)=0;
g1(1)=0;b1(1)=1;c1(1)=0;
g1(Xx+1)=0;b1(Xx+1)=1;c1(Xx+1)=0;
a2(1:Yy+1)=del;
b2(1:Yy+1)=-2*del-2/lambda;
c2(1:Yy+1)=del;
a2(1)=0;b2(1)=1;c2(1)=0;
a2(Yy+1)=0;b2(Yy+1)=1;c2(Yy+1)=0;
f2(1:Yy+1)=0;
uhalf=zeros(nx,ny);
t=0;
%building up the matrix for the tridiagonal
while (t<tn)
for j=2:Yy%visiting the columns to create left side matrix
f1(1)=u(1,j);
f1(Xx+1)=u(Xx+1,j);
for i=2:Xx%visiting the rows tp create left side matrix
f1(i)=-del*u(i,j-1)-del*u(i,j+1)+(2*del-2/lambda)*u(i,j);
end;
uhalf(:,j)=tridiag(Xx+1,g1,b1,c1,f1);
end;
for i=2:Xx
f2(1)=u(i,1);
f2(Yy+1)=u(i,Yy+1);
for j=2:Yy %gitvisiting the columns to create right side matrix
f2(j)=-uhalf(i-1,j)-uhalf(i+1,j)+(2-2/lambda)*uhalf(i,j);
end;
u(i,:)=tridiag(Yy+1,a2,b2,c2,f2);
end;
t=t+dt;

fprintf('\n');
set(gca,'zlim',[0 70]);
hei=surf(x,y,u');
shading interp;
 u(:,1)=u(:,2)-Neuman_ay*dx;%neuman boundary
end;
  