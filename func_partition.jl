# all credit goes to Mark Toth!

function partition(matrix, sub_rows, sub_cols)
    # sub_rows = number of rows in the desired submatrices
    # sub_cols = number of columns in the desired submatrices
    if (size(matrix)[1] % sub_rows != 0) | (size(matrix)[2] % sub_cols != 0)
        error("Please insert feasible numbers for rows and columns of submatrices.")
    end
    n_rows = size(matrix)[1]
    n_cols = size(matrix)[2] 
    parts_collected = Array{Array{Any,2}}(undef, n_rows ÷ sub_rows, n_cols ÷ sub_cols) 
    for i in 1:size(parts_collected)[1]
        for j in 1:size(parts_collected)[2]
            current_row_min = ((i-1) * sub_rows) + 1
            current_row_max = i * sub_rows
            current_col_min = ((j-1) * sub_cols) + 1
            current_col_max = j * sub_cols
            parts_collected[i,j] = matrix[current_row_min:current_row_max, 
                current_col_min:current_col_max]
        end
    end
    parts_collected
end