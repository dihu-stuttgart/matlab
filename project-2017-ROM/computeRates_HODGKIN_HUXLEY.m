function [RATES, ALGEBRAIC] = computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_Stim)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    ALGEBRAIC(:,2) = (  - 0.100000.*(STATES(:,1)+50.0000))./((exp(( - (STATES(:,1)+50.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,6) =  4.00000.*(exp(( - (STATES(:,1)+75.0000)./18.0000)));
    RATES(:,2) =  ALGEBRAIC(:,2).*(1.00000 - STATES(:,2)) -  ALGEBRAIC(:,6).*STATES(:,2);
    ALGEBRAIC(:,3) =  0.0700000.*(exp(( - (STATES(:,1)+75.0000)./20.0000)));
    ALGEBRAIC(:,7) = 1.00000./((exp(( - (STATES(:,1)+45.0000)./10.0000)))+1.00000);
    RATES(:,3) =  ALGEBRAIC(:,3).*(1.00000 - STATES(:,3)) -  ALGEBRAIC(:,7).*STATES(:,3);
    ALGEBRAIC(:,4) = (  - 0.0100000.*(STATES(:,1)+65.0000))./((exp(( - (STATES(:,1)+65.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,8) =  0.125000.*(exp(((STATES(:,1)+75.0000)./80.0000)));
    RATES(:,4) =  ALGEBRAIC(:,4).*(1.00000 - STATES(:,4)) -  ALGEBRAIC(:,8).*STATES(:,4);
    ALGEBRAIC(:,5) =  CONSTANTS(:,3).*(STATES(:,2) .^ 3.00000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6));
    ALGEBRAIC(:,9) =  CONSTANTS(:,4).*(STATES(:,4) .^ 4.00000).*(STATES(:,1) - CONSTANTS(:,7));
    ALGEBRAIC(:,10) =  CONSTANTS(:,5).*(STATES(:,1) - CONSTANTS(:,8));
% tomo
%    ALGEBRAIC(:,1) = piecewise({VOI>=10.0000&VOI<=10.5000, 20.0000 }, 0.00000);
    ALGEBRAIC(:,1) = i_Stim;
% tomo end
    RATES(:,1) =  - ( - ALGEBRAIC(:,1)+ALGEBRAIC(:,5)+ALGEBRAIC(:,9)+ALGEBRAIC(:,10))./CONSTANTS(:,2);
   RATES = RATES';
end