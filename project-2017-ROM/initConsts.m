function [CONSTANTS] = initConsts(model)
    if model=='HODGKIN_HUXLEY'
        VOI = 0; CONSTANTS = [];
        CONSTANTS(:,1) = -75;
        CONSTANTS(:,2) = 1;
        CONSTANTS(:,3) = 120;
        CONSTANTS(:,4) = 36;
        CONSTANTS(:,5) = 0.3;
        CONSTANTS(:,6) = CONSTANTS(:,1)+115.000;
        CONSTANTS(:,7) = CONSTANTS(:,1) - 12.0000;
        CONSTANTS(:,8) = CONSTANTS(:,1)+10.6130;
    end
end