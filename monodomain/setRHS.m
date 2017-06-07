function [bb]=setRHS(method,V_m,time_step,dx,D)

    n=length(V_m);
  
    % right hand side vector
    bb  = zeros(n,1);

    if strcmp(method,'ImplicitEuler')
        % set the right hand side vector
        bb = V_m;      
    elseif strcmp(method,'CN')
        fact = 0.5*time_step/dx^2 * D;
        for j=2:n-1
            bb(j)=(1-2*fact)*V_m(j)+fact*(V_m(j-1)+V_m(j+1));
        end                        
    end
    
    % homogeneous Neumann B.C., no flow         
        bb(1) = 0;% left side       
        bb(n) = 0;% right side
        
end