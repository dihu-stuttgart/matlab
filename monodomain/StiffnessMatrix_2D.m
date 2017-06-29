function [KK]= StiffnessMatrix_2D(method,num_of_points,time_step,dx,D)
% stiffness matrix
  KK = zeros(num_of_points,num_of_points);
  n=sqrt(num_of_points);
  
  % introduce factor for compact notation
  if strcmp(method,'ImplicitEuler')      
      fact = time_step/dx^2 * D;      
  elseif strcmp(method,'CN')      
      fact = 0.5*time_step/dx^2 * D; 
  end
  
  for i=n+1:n*n-n
    if i>n
      KK(i,i-n) = -fact;
    end
    KK(i,i-1) = -fact;
    KK(i,i) =  (1+4*fact);
    KK(i,i+1) = -fact;
    if i<n*n-n+1
      KK(i,i+n) = -fact;
    end
  end
  
  % homogeneous Neumann B.C. - no flow
  % bottom side
  for i=1:n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i+n) = -1;
  end
  % top side
  for i=n*n-n+1:n*n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i-n) = -1;
  end
  % right side
  for i=2*n:n:n*n-n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i-1) = -1;
  end
  % left side
  for i=n+1:n:n*n-2*n+1
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i+1) = -1;
  end
  
end