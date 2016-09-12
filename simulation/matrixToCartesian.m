% Function for Starting Simulation March 2013
% Elena Carano
% April 15, 2013

% Function to transfer a matrix to reflect cartesian coordinates
% Changes values of matrix so that (1,1) information is located in [(0,0)
% (0,delta) (delta,0) (delta,delta)] segment of a cartesian map:
% Translation:
%   - matrix origin located at top left corner
%   - map origin located at bottom left, sometimes negative numbers
%   - matrix furthest point located at bottom right corner
%   - map furthest point located at top right corner
%   - matrix columns are i, rows are j, while map coords are x,y

function [updatedMatrix] = matrixToCartesian(matrix);

updatedMatrix = zeros(size(matrix'));
for i = 1:length(matrix(1,:)) % x coordinates
    for j = 1:length(matrix(:,1)) % y coordinates
        updatedMatrix(end-j+1,i) = matrix(i,j);
    end
end

end