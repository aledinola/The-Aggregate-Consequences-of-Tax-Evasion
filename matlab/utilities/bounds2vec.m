function bounds_vec = bounds2vec(bounds,names)

%{
Purpose: convert a structure with K parameter bounds into a (K,2) array

INPUTS
- bounds: name of the structure. Each field of the structure MUST be a
1*ncol array

- (OPTIONAL) pnames: cell of characters, with K names of parameters to extract

OUTPUT

bounds_vec: it is a (K,ncol) array, typically (K,2)

EXAMPLE

bounds.beta  = [0.92 0.95];
bounds.alpha = [1 2];

vec = bounds2vec(bounds,{'beta';'alpha'})


%}

narginchk(1,2)

if ~isstruct(bounds)
    error('First argument must be a structure.')
end

if nargin==1
    names = fieldnames(bounds);
end


bounds_vec = nan(numel(names),size(bounds.(names{1}),2));

for ii=1:numel(names)
   
   bounds_vec(ii,:) = bounds.(names{ii});
  
end
