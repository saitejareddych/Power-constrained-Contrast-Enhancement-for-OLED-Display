function m = LHM( h, mu )

% Implementation of the Log-based Histogram Modification.
% 
% Parameters: (see the paper for details)
% h		: input histogram vector whose length is 255.
% mu	: determines strength of the modification. Typical values are from 5.0 to 5.5.
% m		: output modified histogram. 

temp = max(h) * 10^(-mu);
denominator = max(h)*temp + 1;
numerator = h*temp + 1;

m = log(numerator) / log(denominator);
m = m/sum(m);		% m is normalized such that the sum of which is to be 1.


end

