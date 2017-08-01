% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic_HODGKIN_HUXLEY(ALGEBRAIC, CONSTANTS, STATES, VOI, i_Stim)
    ALGEBRAIC(:,2) = (  - 0.100000.*(STATES(:,1)+50.0000))./((exp(( - (STATES(:,1)+50.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,6) =  4.00000.*(exp(( - (STATES(:,1)+75.0000)./18.0000)));
    ALGEBRAIC(:,3) =  0.0700000.*(exp(( - (STATES(:,1)+75.0000)./20.0000)));
    ALGEBRAIC(:,7) = 1.00000./((exp(( - (STATES(:,1)+45.0000)./10.0000)))+1.00000);
    ALGEBRAIC(:,4) = (  - 0.0100000.*(STATES(:,1)+65.0000))./((exp(( - (STATES(:,1)+65.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,8) =  0.125000.*(exp(((STATES(:,1)+75.0000)./80.0000)));
    ALGEBRAIC(:,5) =  CONSTANTS(:,3).*(STATES(:,2) .^ 3.00000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6));
    ALGEBRAIC(:,9) =  CONSTANTS(:,4).*(STATES(:,4) .^ 4.00000).*(STATES(:,1) - CONSTANTS(:,7));
    ALGEBRAIC(:,10) =  CONSTANTS(:,5).*(STATES(:,1) - CONSTANTS(:,8));
% tomo;
    ALGEBRAIC(:,1) = i_Stim;
% tomo end
end