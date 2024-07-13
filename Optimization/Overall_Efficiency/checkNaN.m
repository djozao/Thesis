function checkNaN(value, variableName, Temperature, Pressure)
    if ~isnan(value)
        return
    end
    conditionStr = '';
    if ~isempty(Temperature)
        conditionStr = sprintf('at Temperature %.2f K', Temperature);
    end
    if ~isempty(Pressure)
        if ~isempty(conditionStr)
            conditionStr = [conditionStr, ' and '];
        end
        conditionStr = [conditionStr, sprintf('at Pressure %.2f Pa', Pressure)];
    end
    
    if isnan(value)
        if isempty(conditionStr)
            error('Interpolated value for %s is NaN. Check input range.', variableName);
        else
            error('Interpolated value for %s %s is NaN. Check input range.', variableName, conditionStr);
        end
    end
end