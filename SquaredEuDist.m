function D = SquaredEuDist(X,Y)
%EUDIST  Compute the squared Euclidean distance matrix.
% Problem description: D_ij = ||x_i-y_j||^2
% 
% Usgae: 
%     D = SquaredEuDist(X,Y)
% Input:
%     - X: d*nX data matrix with nX samples
%     - Y: d*nY data matrix with nY samples
% Output:
%     - D: nX*nY distance matrix of X and Y
% 
% 2021/9/19 



% Main
X_sum_row = (sum(X.^2, 1));
Y_sum_row = (sum(Y.^2, 1));
D = bsxfun(@minus, X_sum_row',bsxfun(@minus,2*X'*Y, Y_sum_row));

end

