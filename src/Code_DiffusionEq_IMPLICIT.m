%MECE 5397 PROJECT
%Oscar Rodriguez
%10717156   
%
%Project B- Diffusion Equation
%Bc2-4

%Solving Equation Using Implicit Method

%Defining Parameters

D=1;            %Given difrusivity       
ax=-3.14;       %Lower x Limit
ay=-3.14;       %Lower y Limit  
bx=3.14;        %Upper x limit
by=3.14;        %Upper y Limit

Nxx=39;                          %Number of steps in space(x)
Nyy=39;                          %Number of steps in space(y)  

nt=100;                         %Number of time steps 
delta_t=.00065;                    %size time step
h1=6.28/(Nxx);                  %Size of space step(x)

u=zeros(Nxx+1,Nyy+1);          %Creating the grid for temperature        
un=zeros(Nxx+1,Nyy+1);         %initioal grid for temperature, Initial time condition is also represented by this line
x=ax:h1:bx;                      %defining the range of x
y=ay:h1:by;                      %defining the range of y



