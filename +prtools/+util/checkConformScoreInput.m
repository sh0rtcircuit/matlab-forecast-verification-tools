function [probDists,measurements] = checkConformScoreInput(probDists,measurements)
%checkConformScoreInput Function that validates the input for score
%computation.

    % conform dimensions
    probDists = reshape(probDists,numel(probDists),1);
    measurements = reshape(measurements,numel(measurements),1);
    
    % check for dimension equality
    if length(probDists) ~= length(measurements)
        error(['Error: Unequal number of pdfs and measurements: handles: ' ...
            num2str(size(probDists,1)) ', Measurements: ' num2str(size(measurements,1))]);
    end
    
    % check that all elements are objects
    if any(~isobject(probDists{1}))
        error('Error: ProbDists are no classes, has to be probability distribution class');
    end
    
    % check that classes are indeed probability distributions
    % (simplification)
    cName = class(probDists{1});
    if ~strcmp(cName(1:5),'prob.')
       error('Class has to be a probability distribution (e.g., created using ''makedist'' '); 
    end
end

