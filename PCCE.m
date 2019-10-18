function [y, conv_y1, counter] = PCCE(m, h, beta, gamma, ismono, ini_y1)

%%%%%%%%%%%%%%%%%%%
% Input Variables
%%%%%%%%%%%%%%%%%%%
% m : modified histogram or guide histogram
% h : input original histogram

%%%%%%%%%%%%%%%%%%%
% Input Variables
%%%%%%%%%%%%%%%%%%%
% y 		: optimized result
% conv_y1 	: first element of y when the iteration is converged
% counter	: log for the number of secant iterations

%%%%%%%%%%%%%%%%%%%
% Input Parameters
%%%%%%%%%%%%%%%%%%%
% beta  	: regulation parameter between power term and contrast term. See Eq. (28)
% gamma		: gamma value for the power modeling in Section 3-A.
% ismono	: should the transformation fuction be monotonic? Option for monotonic constraint

%%%%%%%%%%%%%%%%%%%
% Other parameters 
%%%%%%%%%%%%%%%%%%%
% L		: original input length, ex) 256
% niter : the number of iterations for the secant analysis
% toll	: tollerance for escaping the loop
% ips	: very small number for the secant method

if nargin < 4
    gamma = 2.2;
end
if nargin < 5
    ismono = 1;
end
if nargin < 6
    y1_tmp = -100.0;
else
    y1_tmp = ini_y1;
end

% default parameter setting
niter = 100;
toll = 5e-11;
ips = 1e-8;
L = length(h);

% for final y
y = zeros(L,1);

% how many iterations are there - log
counter = zeros(L-1,1);

m = m/sum(m) * (L-1);
alpha = beta/([0:L-1]*h);

for k = 1:L-1
    
    y(k,1) = 0;
    
    h_tmp = h(k+1:end,1);
    m_bar_tmp = m(k+1:end,1);   % m_bar_tmp = (L-1) * m_bar_tmp/sum(m_bar_tmp);
    
    y1 = zeros(niter, 1);
    f = zeros(niter, 1);
    
    if y1_tmp == -100.0
        y1(1,1) = m_bar_tmp(1,1);
        y1(2,1) = m_bar_tmp(1,1) - ips;
    else
        y1(1,1) = y1_tmp;
        y1(2,1) = y1_tmp - ips;
    end
    
    for n=3:niter
        f(n-2,1) = sum(calc_const(y1(n-2,1), h_tmp, m_bar_tmp, alpha, gamma, ismono)) - (L-1);
        f(n-1,1) = sum(calc_const(y1(n-1,1), h_tmp, m_bar_tmp, alpha, gamma, ismono)) - (L-1);
        
        y1(n,1) = y1(n-1,1) - ( (y1(n-1,1) - y1(n-2,1))/(f(n-1,1) - f(n-2,1)) * f(n-1,1) );
        
        % counter increment
        counter(k,1) = counter(k,1) + 1;
        
        if abs(y1(n,1) - y1(n-1,1)) < toll
            break
        end
        
    end
    
    if (y1(n,1) >= 0) || (ismono == 0)
        y(k+1:end,1) = calc_const(y1(n,1), h_tmp, m_bar_tmp, alpha, gamma, ismono);
        
        conv_y1 = y1(n,1);
        break
    end    
end

end

function y = calc_const(y1, h, m, alpha, gamma, ismono)


L = length(h);
y_tmp = zeros(L,1); y_tmp(1,1) = y1;
y = zeros(L,1); y(1,1) = y1;

h_tmp = h;
m_tmp = m;
    
for i=1:L-1
    y_tmp(i+1,1) = m_tmp(i+1,1) - m_tmp(i,1) + y_tmp(i,1) + (alpha*gamma/2)*h_tmp(i,1)*real(sum(y(1:i))^(gamma-1));
    y(i+1) = y_tmp(i+1,1);
    
    if ismono && y_tmp(i+1,1) <=0
        y(i+1) = 0;
%         y_tmp(i+1) = y_tmp(i);
%         h_tmp(i+1) = h_tmp(i) + h_tmp(i+1);
%         m_tmp(i+1) = m_tmp(i);
    end
end


end

