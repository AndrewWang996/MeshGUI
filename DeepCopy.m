function new = DeepCopy(object)
    % Instantiate new object of the same class.
    
    new = feval(class(object));

    % Copy all non-hidden properties.
    p = properties(object);
    for i = 1:length(p)
        try 
            new.(p{i}) = object.(p{i});
        catch
            
        end
    end

end