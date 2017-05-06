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

Nxx=30;                          %Number of steps in space(x)
Nyy=30;                          %Number of steps in space(y)    

h1=6.28/(Nxx-1);
delta_t=.999*.25*(h1^2);  
nt=1000;                         %Number of time steps 
t=0:delta_t:nt;
                

x=ax:h1:bx;                      %defining the range of x
y=ay:h1:by;                      %defining the range of y

                  
fa=y.*(y-ay).^2 ;                   %defining f(a)
ga=((y-ay).^2).*cos(-y); %defining g(a)
Neuman_ay=0;               %defning the neumand condition given by the problem             
u=zeros(Nxx,Nyy);          %Creating the grid for temperature        
un=zeros(Nxx,Nyy);         %initioal grid for temperature, Initial time condition is also represented by this line

u(1,:)=ga; %detriclet non linear  condition in the  right side boundary
u(Nxx,:)=fa; %detriclet condition non linear condition left side boudary 
u(:,Nyy)=fa(Nxx)+((x+3.14)/(6.28)*(ga(Nxx)-fa(Nxx))); %detrichlet condition upper boundary 
for it=0:nt                %Begin time loop
    un=u; %definite initial temperature (initial time condition)
   
    %% Explicit Method Iteration
     
   for j=2:Nyy-1 %initiation of iterative process to calculate temperature (y axis)
    for i=2:Nxx-1%initiation of iterative process to calculate temperature (y axis)
    u(i,j)=un(i,j)+(D*delta_t*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(h1*h1))+(D*delta_t*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(h1*h1));
        
    u(:,1)=u(:,2)-Neuman_ay*h1; %Neuman condition bottom boundary
    
    end %end for loop for iteration process of x
   end  %end for loop for iteration process of y 
    hei=surf(x,y,u);      %Initiation of 3 dimensional graph capable of showing time dependece of problem
    shading interp% color gradient for graph
    xlabel('(x)') %lable axixs
    ylabel('(y)') %lable axis
    zlabel('Temperature') %lable axis
   
    title({['2-D Diffusion with D = ',num2str(D)];['time (\itt) = ',num2str(it*delta_t)]})
      
    drawnow;
    
    
    refreshdata(hei) %frefresh of graph per time interval
end %end time loop