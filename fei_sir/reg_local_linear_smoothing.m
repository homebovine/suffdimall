% Local linear regression
% Author: Xiaohui Chen (xiaohuic@ece.ubc.ca)
% Version: 2012-Feb

function [m, m_0] = reg_local_linear_smoothing(x, y, kerfun, h, x_0),

% Local Linear Smoothing (LLS) for 1-dim curve fitting.
%
% Input:
% x -- independent variable
% y -- dependent variable
% kerfun -- kernel function with options: 
% kerfun = 0: Uniform
% kerfun = 1: Gaussian
% kerfun = 2: Epanechnikov
% kerfun = 3: Triangular
% kerfun = 4: Quartic (biweight)
% kerfun = 5: Triweight (tricube) kernel
% kerfun = 6: Cosine
% h -- bandwidth of kernel function
% x_0 -- points-of-interests (optional)
%
% Output:
% m -- smoothed mean vector of size n*1
% m_0 -- fitted values at x_0
%
% Ref:
% Fan (1993) Local linear regression smoothers and their mimimax
% efficiencies. Annals of Statistics, 21(1): 196-216.

n = length(y);
if nargin < 5, x_0 = 1; end     % Default value is the end of window
if nargin < 4, h = 1; end      % Default bandwidth = 1
n0 = length(x_0);

% Compute weights with the kernel function
K = zeros(n, n);
W = zeros(n, n);
s_1 = zeros(n, 1);
s_2 = zeros(n, 1);

K_0 = zeros(n, n0);
W0 = zeros(n, n0);
s0_1 = zeros(n0, 1);
s0_2 = zeros(n0, 1);

if 0 == kerfun,     % Uniform kernel
    for i = 1:n,
        K(:,i) = .5 * ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = .5 * ( abs( x_0(i0) - x ) <= h );
    end
elseif 1 == kerfun,     % Gaussian kernel
    for i = 1:n,
        K(:,i) = normpdf( ( x(i) - x ) / h, 0, 1 );
    end
    for i0 = 1:n0,
        K_0(:,i0) = normpdf( ( x_0(i0) - x ) / h, 0, 1 );
    end
elseif 2 == kerfun,           % Epanechnikov kernel
    for i = 1:n,
        K(:,i) = 3/4 * ( 1 - ( ( x(i) - x ) / h ).^2 ) .* ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = 3/4 * ( 1 - ( ( x_0(i0) - x ) / h ).^2 ) .* ( abs( x_0(i0) - x ) <= h );
    end
elseif 3 == kerfun,         % Triangular kernel
    for i = 1:n,
        K(:,i) = ( 1 - abs( ( x(i) - x ) / h ) ) .* ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = ( 1 - abs( ( x_0(i0) - x ) / h ) ) .* ( abs( x_0(i0) - x ) <= h );
    end
elseif 4 == kerfun,         % Quartic (biweight) kernel
    for i = 1:n,
        K(:,i) = 15/16 * ( 1 - ( ( x(i) - x ) / h ).^2 ).^2 .* ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = 15/16 * ( 1 - ( ( x_0(i0) - x ) / h ).^2 ).^2 .* ( abs( x_0(i0) - x ) <= h );
    end
elseif 5 == kerfun,         % Triweight (tricube) kernel
    for i = 1:n,
        K(:,i) = 35/32 * ( 1 - ( ( x(i) - x ) / h ).^2 ).^3 .* ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = 35/32 * ( 1 - ( ( x_0(i0) - x ) / h ).^2 ).^3 .* ( abs( x_0(i0) - x ) <= h );
    end
elseif 6 == kerfun,         % Cosine kernel
    for i = 1:n,
        K(:,i) = pi/4 * ( cos(pi / 2 * ( x(i) - x ) / h ) ) .* ( abs( x(i) - x ) <= h );
    end
    for i0 = 1:n0,
        K_0(:,i0) = pi/4 * ( cos(pi / 2 * ( x_0(i0) - x ) / h ) ) .* ( abs( x_0(i0) - x ) <= h );
    end
end
for i = 1:n,
    s_1(i) = sum( ( x(i) - x ) .* K(:,i) );
    s_2(i) = sum( ( x(i) - x ).^2 .* K(:,i) );
    W(:,i) = K(:,i) .* ( s_2(i) - ( x(i) - x ) * s_1(i) );
end
for i0 = 1:n0,
s0_1(i0) = sum( ( x_0(i0) - x ) .* K_0(:,i0) );
s0_2(i0) = sum( ( x_0(i0) - x ).^2 .* K_0(:,i0) );
W0(:,i0) = K_0(:,i0) .* ( s0_2(i0) - ( x_0(i0) - x ) * s0_1(i0) );
end

% Smoothing
if any( 0 == sum(W) ),
    m = ( y' * W ./ ( sum(W) + n^-2 ) )';       % Avoid dividing zero
else
    m = ( y' * W ./ sum(W) )';
end

if any( 0 == sum(W0) ),
    m_0 = ( y' * W0 ./ ( sum(W0) + n^-2 ) )';
else
    m_0 = ( y' * W0 ./ sum(W0) )';
end

end
