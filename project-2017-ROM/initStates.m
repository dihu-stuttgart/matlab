function [STATES] = initStates(model)
    if(model=='HODGKIN_HUXLEY')
        STATES(:,1) = -75;
        STATES(:,2) = 0.05;
        STATES(:,3) = 0.6;
        STATES(:,4) = 0.325;
        if (isempty(STATES)), warning('Initial values for states not set'); end
    end
end
