%%% Hankel Matric Construction
%%% Allen Lee

% Each column of y must be a data from the same time instant!
function H = Hankel_Matrix(y,H_row,H_col,Output_dim,Input_dim,index)
[m,data_length] = size(y);

if((H_col+H_row+index-1)>data_length)
    warning("Extend data samples or decrease matrix dimension!")
else
    H = [];
    % Collect data from index 
    for i=index:H_col+index-1 % for each column
        h = [y(:,Input_dim*(i-1)+1:i*Input_dim)];%h=[y(:,i)];
        for j = 1:H_row-1 % for each row
            h = [h;y(:,Input_dim*(i+j-1)+1:(i+j)*Input_dim)];
        end
        H = [H h];
    end

end