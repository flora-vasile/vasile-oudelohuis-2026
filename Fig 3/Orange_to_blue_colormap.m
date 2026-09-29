n = 256;  % number of colors
blue = [0 0.2 0.7];
white = [1 1 1];
orange = [0.9 0.4 0];

% Interpolate from blue ? white ? orange
cmap = [interp1([1 n/2], [blue; white], 1:n/2); ...
        interp1([1 n/2], [white; orange], 1:n/2)];