function [algebraicVariableCount] = getAlgebraicVariableCount(model) 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".
    if(model=='HODGKIN_HUXLEY')
        % There are a total of 4 entries in each of the rate and state variable arrays.
        % There are a total of 8 entries in the constant variable array.
        algebraicVariableCount =10;
    end
end