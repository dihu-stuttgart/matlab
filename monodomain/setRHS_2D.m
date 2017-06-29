function [bb]=setRHS_2D(method,V_m,time_step,dx,D)
    num_of_points=length(V_m);
    n=sqrt(num_of_points);
  
    % right hand side vector
    bb  = zeros(num_of_points,1);

    if strcmp(method,'ImplicitEuler')
        % set the right hand side vector
        bb = V_m;      
    elseif strcmp(method,'CN')
        fact = 0.5*time_step/dx^2 * D;
        for i=n+1:n*n-n
            bb(i) =  (1-4*fact)*V_m(i)+fact*(V_m(i-1)+V_m(i+1)+V_m(i-n)+V_m(i+n));   
        end
    end
           
    % homogeneous Neumann B.C., no flow 
    % bottom side
    for i=1:n
      bb(i) = 0;
    end
    % top side
    for i=n*n-n+1:n*n
      bb(i) = 0;
    end
    % right side
    for i=2*n:n:n*n-n
      bb(i) = 0;
    end
    % left side
    for i=n+1:n:n*n-2*n+1
      bb(i) = 0;
    end
        
end