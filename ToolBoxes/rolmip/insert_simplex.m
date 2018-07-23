function [varargout] = insert_simplex(varargin)
% Procedure to insert a new multisimplex variable in an already existent
% polynomial variable
%
% [poly] = insert_simplex(poly,vertices)
% Inserts, in the variable poly, the multisimplex with the given number of
% vertices, which corresponds to a row vector.
%
% [poly] = insert_simplex(label,vertices)
% Inserts the new multisimplex in the polynomial with the given label.

if (nargin ~= 2)
    error('Input error. Type ''help insert_simplex'' for more details');
    return
end

if ischar(varargin{1})
    label = varargin{1};
    poly = rolmip('getvar',label);
    if isempty(poly)
        error('The informed label does not correspond to any polynomial');
        return
    end
else
    poly = varargin{1};
    label = poly.label;
end

vertices = varargin{2};

poly.vertices = [poly.vertices vertices];

numsimplex = length(poly.data(1).exponent);
for cont=1:length(poly.data)
    for continsert=1:length(vertices)
        poly.data(cont).exponent{numsimplex+continsert} = [0];
    end
end


varargout{1} = poly;



return