clc;
clear all;

%size of the problem
n=9;

%levels of refinement
L=round(log2(n/2));

%original plus L levels of coefficients to store
m=L+1;

%maximum additional length of zeros before and after the relevant coefficients
sd_max=2^L;

%analytical function for rhs
x=-1:2/(n-1):1; %1D domain discretization
y=f(x);

%start and end of relevant indices 
i_start=sd_max+1;
i_end=i_start+n-1;

CL=zeros(n+2*sd_max,m);%additional zeros befor and after the entries
CL(i_start+1:i_end,:)=-1;

CC=zeros(n+2*sd_max,m);
CC(i_start:i_end,:)=2;
for i=i_end+1:n+2*sd_max
    CC(i,:)=1;
end

CU=zeros(n+2*sd_max,m);
CU(i_start:i_end-1,:)=-1;

d=zeros(n+2*sd_max,m);
d(i_start:i_end,1)=y;

alpha=matlabsolve(i_start,i_end,CL(:,1),CC(:,1),CU(:,1),d(:,1));
fprintf('alpha_mat=%f\n',alpha);
fprintf('\n');

alpha=CyclicReduction(i_start,i_end,CL,CC,CU,d);
fprintf('alpha=%f\n',alpha(i_start:i_end));

function [x]= CyclicReduction( i_start,i_end,CL,CC,CU,d )
%Cyclic reduction for parallel Thomas algorithm

%number of unknowns
n=i_end-i_start+1;

%level of refinement
L=round(log2(n/2));
m=L+1;
sd_max=2^L;

%Indices of interest are between 2 and n+1.
x=zeros(n+2*sd_max,1);

%stride in arrays
sd=1;

%eliminate equations and reduce the system of unknowons to two in L levels
for ll=1:L            
    %Update coefficients    
    for i=(i_start-1)+2*sd:2*sd:i_end
        fprintf('forward reduction:\tl=%d,\ti=%d\n',ll,i);
        
        k1=CL(i,ll)/CC(i-sd,ll);
        k2=CU(i,ll)/CC(i+sd,ll);
    
        CL(i,ll+1)=-k1*CL(i-sd,ll);
        CC(i,ll+1)=CC(i,ll)-k1*CU(i-sd,ll)-k2*CL(i+sd,ll);
        CU(i,ll+1)=-k2*CU(i+sd,ll);
        d(i,ll+1)=d(i,ll)-k1*d(i-sd,ll)-k2*d(i+sd,ll);
    end    
    sd=sd*2;
    fprintf('\n');
end

i_back=i;

%Thomas algorithm
j=floor(n/2);
if(mod(j,2)==0)
    CU(i_back-sd,m)=CU(i_back-sd,m)/CC(i_back-sd,m);
    CU(i_back,m)=CU(i_back,m)/(CC(i_back,m)-CU(i_back-sd,m)*CL(i_back,m));
    
    d(i_back-sd,m)=d(i_back-sd,m)/CC(i_back-sd,m);
    d(i_back,m)=(d(i_back,m)-d(i_back-sd,m)*CL(i_back,m))/(CC(i_back,m)-CU(i_back-sd,m)*CL(i_back,m));
    
    x(i_back)=d(i_back,m);
    x(i_back-sd)=d(i_back-sd,m)-CU(i_back-sd,m)*x(i_back);
else
    x(i_back)=d(i_back,m)/CC(i_back,m);
end
 
for ll=L:-1:2        
    sd=sd/2;
    
    %last equation in level ll
    if i_end >= i_back+sd
        i_back=i_back+sd;
    else
        i_back=i_back-sd;
    end    
    
    %all unknowns with odd indices in the reduction levels are found 
     for i=i_back:-2*sd:i_start   
        fprintf('backward substitution:\tl=%d,\ti=%d\n',ll,i)
        x(i)=(d(i,ll)-CL(i,ll)*x(i-sd)-CU(i,ll)*x(i+sd))/CC(i,ll);    
     end
     fprintf('\n');
end  

%all other unknowns with stride of 2 arw found
for i=i_start:2:i_end
    x(i)=(d(i)-CL(i)*x(i-1)-CU(i)*x(i+1))/CC(i);
end
       
end

function [alpha]=matlabsolve(i_start,i_end,CL,CC,CU,d)

A=diag(CL(i_start+1:i_end),-1)+diag(CC(i_start:i_end))+diag(CU(i_start:i_end-1),1);
alpha=A\d(i_start:i_end);

end

function [y]=f(x)
y=-x.^2+1;
end


