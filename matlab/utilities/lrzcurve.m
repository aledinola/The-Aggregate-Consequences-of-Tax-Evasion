function [pw, meanw, stdvw, giniw] = lrzcurve(p, w, color)

% INPUTS:
% p is the probability distrib (sum to 1)
% w is a vector with the values of the variable of interest
% (for example, the policy function for bequests/wealth)

% OUTPUTS:
% pw is a [n,2] vector with (x,y) to plot Lorenz curve


%% Eliminate elements with zero probability
p1=p; p(p1==0)  = []; w(p1==0)  = [];

%% Standard Deviation
meanw     = w'*p;
stdvw  = sqrt((w'-meanw).^2*p);

%% compute relative value for each "wealth" level
wi=p.*w; wt=sum(wi); wr=wi/wt;

%% Keep the smallest population, needed to normalize the Gini coefficient
minpop = min(p);

%% Store in a single array
pw = [p, wr, w];

%% Sort with respect to "Wealth" Per Capita
pw = sortrows(pw, 3);
pw(:, 3) = [];
pw = [zeros(1, 2); pw];

%% Cumulative p & w
pw = cumsum(pw);


%% Average bases and height for right trapezoids
height = diff(pw(:, 1));
base = (pw(1:(end - 1), 2) + pw(2:end, 2)) / 2;

% The Gini Coefficient is normalized with respect to its highest possible
% value which is obtained if the smallest population owns all the existing
% wealth.
giniw = (1 - 2 * height' * base) ./ (1 - minpop);

% plot(pw(:, 1), pw(:, 2), 'LineWidth', 4, 'Color', color);            % Lorenz Curve

