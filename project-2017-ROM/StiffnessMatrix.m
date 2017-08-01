function [KK]= StiffnessMatrix(method,num_of_points,time_step,dx,D)
% stiffness matrix
  KK = zeros(num_of_points,num_of_points);
  n=num_of_points;
    
  % introduce factor for compact notation
  if strcmp(method,'ImplicitEuler')      
      fact = time_step/dx^2 * D;      
  elseif strcmp(method,'CN')      
      fact = 0.5*time_step/dx^2 * D; 
  end
  
  for i=2:n-1
      KK(i,i-1) = -fact;
      KK(i,i) =  (1+2*fact);
      KK(i,i+1) = -fact;
  end
  
  % homogeneous Neumann B.C. - no flow
  % right side
    KK(n,n) = 1;
    KK(n,n-1) = -1;
  % left side
    KK(1,1) = 1;
    KK(1,2) = -1;
  
  %dlmwrite('Matrix_2D.csv',KK,'delimiter','\t', 'precision',3);
  
end