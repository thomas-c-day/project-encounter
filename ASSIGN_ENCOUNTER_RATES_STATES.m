function GAMMA = ASSIGN_ENCOUNTER_RATES_STATES(n, b, G0, alpha, phi, Dct)
    % Assign a mean encounter rate according to the current size of the
    % microbial aggregate.

    th_current = zeros(1,length(n));
    DiffConstant = th_current;
    for ii = 1:length(n)
        if b(ii) <= 2 % this means they are still single cells
            th_current(ii) = 0.5 * b(ii).^(1/3);            % this returns the current cell radius
            % DiffConstant(ii) = Dct;
        else % this means they are more than a single cell
            th_current(ii) = 0.5 * (b(ii)/(phi)).^(1/3);    % this returns the current aggregate radius
            % DiffConstant(ii) = 0;
        end

        % Single cell state swims, multi-celled state doesn't:
        if n(ii) == 1
            DiffConstant(ii) = Dct;
        else
            DiffConstant(ii) = 0;
        end
    end
    GAMMA = G0 * (th_current + 0.5).^alpha + 4*pi * DiffConstant .* (th_current + 0.5);

end