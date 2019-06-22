%{
Objective: lower the rank of given matrix by minimizing the nulcear norm
- input: a matrix parsed from a csv file, which contains many NaNs(unknown entries); desired accuracy;
- output: parse result back to a csv file
- NaN entries changes in every loop; Given entries remains unchanged after every loop
%}
function optimize(in_mat, accuracy)
    curr = zeros(50,200);
    new = fillWithZero(in_mat);
    indices = find(isnan(in_mat)); %find all the NaNs; the entries needed to be filled in
    
    while norm(new - curr) > accuracy
        curr = new;
        [U, D, V] = svd(curr);
        minv = min(svds(curr)); %get the minimum singular value
        disp(minv);
        
        %soft D
        for i = 1:size(D) 
            sigma = D(i,i) - minv; %subtract the minimum singular value off each diagonal entry
            if (sigma < 0)
                sigma = 0;
            end
            D(i,i) = sigma;
        end       
        disp(D(3,3)); %check if diagonal entries are gradually decreasing

        new_unaligned = U*D*V';
        new = combineMat(curr,new_unaligned,indices); %Given entries remains unchanged after every loop
             
    end 
    disp(new);
    parse_back('res1',new);
end 


%Objective: find all the Nans and set them to zero
function in_mat = fillWithZero(in_mat)
    indices = find(isnan(in_mat)); %find all the NaNs
    
    for i = 1 : size(indices)
        [row, col] = ind2sub(size(in_mat), indices(i)); %convert an index to [row,col] tuple
        in_mat(row, col) = 0;
    end
end



%{
- Objective: 
  Fill the unknown entries with new data; 
  while given entries remains unchanged
- indices: the indices of unknown entries. 
- out_mat: output matrix, already contains known entries
- un_mat: unaligned matrix--known entries not yet been aligned
%}
function out_mat = combineMat(out_mat, un_mat, indices)
    for i = 1 : size(indices)
        [row, col] = ind2sub(size(out_mat), indices(i)); % %convert an index to [row,col] tuple
        out_mat(row, col) = un_mat(row, col);
    end
end
