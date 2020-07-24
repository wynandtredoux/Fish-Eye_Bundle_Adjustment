function col = printCell(fileID,Mcell,prefix,padding)
    % find longest string in 1st column
    max = 0;
    for i = 1:size(Mcell,1)
        if length(Mcell{i,1}) > max
            max = length(Mcell{i,1});
        end
    end
    
    col = max+padding+3;
    for i = 1:size(Mcell,1)
        % print line ------
        if strcmp(Mcell{i,1},'\line')
            for j = 1:(max+padding+6)
                fprintf(fileID,'-');
            end
        elseif strcmp(Mcell{i,1},'\n') % line break \n
            % do nothing
        else        
            % print prefix
            fprintf(fileID, prefix);        
            % print 1st col
            fprintf(fileID, [Mcell{i,1} ' ']);
            % print padding
            for j = 1:(max+padding-length(Mcell{i,1}))
                fprintf(fileID,'.');
            end
            fprintf(fileID,' ');
            % 2nd col
            if ischar(Mcell{i,2})
                fprintf(fileID, Mcell{i,2});
            elseif isa(Mcell{i,2},'double')
                fprintf(fileID, num2str(Mcell{i,2}));
            else
                disp('Error printCell(): invalid datatype in Mcell. Must be char or double')
            end
        end
        % new line
        fprintf(fileID,'\n');
    end
end