function [a,g,R] = getlpc(data,order)
% function [a,g,R] = getlpc(data,order)
% finds the inverse filter coefficients, a=[1, a(1), ..., a(N)], of an N:th
% order linear predictor.
%
% Inputs:
%       data:       signal frame
%       order:      order of the LPC (default = 20)
% Outputs:
%       a:          LPC inverse filter coefficients
%       g:          gain 
%       R:          autocorrelation matrix

if nargin<2 
    order = 20;
end

N = length(data);

% Calculate autocorrelation for lags [0,1, ..., p] 
% See Eq. 4.10 in instructions.
% (Note: pay attention to indexing in equations vs. in MATLAB vectors that
% always start from index 1).
p = order;

% doing it for a window

% r = ?
r = zeros(p+1, 1);

for i = 0:p
    for n = i:N-1
        r(i+1,:) = r(i+1,:) + data(n+1)*data(n-i+1);
    end
end
r;

% Create autocorrelation matrix from autocorrelation vector (Eq. 4.11)
R = toeplitz(r(1:end-1)); 
R;

% Solve Ra = r (Eq. 4.12) for LPC prediction coefficients (Eq. 4.13)
% a = ? 
a = R\r(2:end);

% Gain calculation (Eq. 4.14)
% g = ?
temp_sum = 0;
for k=1:p
    temp_sum = temp_sum + a(k)*r(k+1);
end 
g = sqrt(r(1)-temp_sum);

% Convert prediction coefficients to synthesis coefficients
% (Eq. 4.15)
a = [1;-a(1:order)];


