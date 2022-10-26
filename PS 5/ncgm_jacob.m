function jacob = ncgm_jacob(x0,  params)
% This function computes a numerical approximation of a jacobian matrix at
% a specific point (combination of values)
% x0: vector of values to evaluate the function at
% params: set of parameters, to be fed into the function


% Set matrix of appropriate size
jacob = zeros(length(x0), length(x0));

% Loop over each input variable
for i = 1:length(x0)
    dx = zeros(length(x0),1); % prepare vector of variable changes
    dx(i) = x0(i)*1e-3; % set variable-specific change
    x1 = x0 + dx; % compute new input value
    %Calculate numerical approximation to jacobian by forward difference
    jacob(:,i) = (rbc_obj_start(x1, params)- rbc_obj_start(x0, params))/ dx(i);
end
