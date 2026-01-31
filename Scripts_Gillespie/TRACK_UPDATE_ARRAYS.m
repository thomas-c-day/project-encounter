function [Conc_output, ix_ind, Bank_output] = TRACK_UPDATE_ARRAYS(Bank, Rates, Conc_obj, Water_Vol)

    % Update the concentration of items, and assign the object to a
    % particular individual. Then, return the bank of objects encountered
    % by the individual.

    % Update the concentration array:
    num_obj = Conc_obj * Water_Vol;
    num_obj = num_obj - 1;
    Conc_output = num_obj / Water_Vol;

    % Choose which individual caught the object:
    probs = Rates./sum(Rates);
    ix_ind = randsample(length(Rates), 1, true, probs);

    % Update the bank:
    Bank_output = Bank;
    Bank_output(ix_ind) = Bank(ix_ind) + 1;

end