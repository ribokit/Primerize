function lottery_96_well (total_num, pick_num)

%
% LOTTERY_96_WELL (total_num, pick_num)
%
% Generates non-repeat random pick_num wells for NanoDrop (or other) 
%   sampling from total_num wells. Output format "plate/row/column", e.g. 
%   "1A3" means plate 1, row A, column 3.
%
% Input
% =====
%   total_num           Number of total available wells on the plate(s).
%                           Multiple plates allowed.
%   pick_num            Number of sampling selection.
%
% by T47, Oct 2013.
%


row_char = 'ABCDEFGH';
no = randperm(total_num);
no = sort(no(1:pick_num));

for i = 1:pick_num
    plate_ind = ceil(no(i)/96);
    row_ind =  mod(no(i)-1, 8)+1;
    column_ind = ceil((no(i) - (plate_ind -1) * 96)/8);
    fprintf([num2str(plate_ind),row_char(row_ind),num2str(column_ind),'\t']);
    if mod(i,10) == 0; fprintf('\n'); end;
end;
fprintf('\n');