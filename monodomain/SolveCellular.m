function [vS_hh,ALL_STATES]= SolveCellular(order,time,num_of_points,ALL_STATES,CONSTANTS,i_Stim,tspan)
vS_hh = zeros(num_of_points,1);

% Integrate the cellular model at each discretisation point
    for point = 1:num_of_points
      
      % load the last state as initial state for the next integration
      LAST_STATES = ALL_STATES(point,:);
      
      %---------------------------------------------------------------------
      % INTEGRATE CELLULAR MODEL WITH ODE/DAE SOLVER
      % NOTE: no timestep for the integration of the cellular model has to
      % be provided. MATLAB chooses an appropriate time step
      if(order==1)
      [STATES] = ode1(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_Stim(point,time)), tspan, LAST_STATES); 
      else
          if(order==2)
      [STATES] = ode2(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_Stim(point,time)), tspan, LAST_STATES);
          else
              fprintf('orders 1 and 2 are supported');
          end
      end
      
      % update ALL_STATES
      ALL_STATES(point,:) = STATES(end,:);
      % update transmembrane voltage of cellular model
      vS_hh(point) = STATES(end,1);
      
    end % points

end