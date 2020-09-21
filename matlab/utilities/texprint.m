% texprint prints matrix in tex table format for importing into SW
% A is the matrix to be printed
% fmt is a format string. For now, same for all elements of A Typical is
% ' %6.2f ' ' %6.2G '. Dont forget spaces
% head is optional header. It's either a character array or a cell array of strings
% row is optional row lables. Again a character array or a cell array of strings

% for example, 
% head = { 'col 1'; 'col 2'; 'column 3'};
% A = [1 2 3; 4 5 6]; 
% fmt = ' %G6.2 '; 
% row = { 'row 1'; 'row 2 label'}; 
% texprint(A, fmt, head, row);

% you can omit head and row, but they must be empty, e.g. 
% texprint(A, fmt, head, []); 

function texprint(A, fmt, head, row); 

if iscellstr(head);
    head = char(head); 
end; 
if iscellstr(row);
    row = char(row);
end; 

if (size(head,1)~= size(A,2))&(~isempty(head)) ; 
    disp('texprint: rows of header not equal to colums of matrix'); 
end; 
if (size(row,1)~=size(A,1))&(~isempty(row)); 
    disp('texprint: rows of row head not equal to rows of matrix'); 
end; 

% display \begin{tabular}{llll}
els = 'l';
fmtstr = '%c';
for i=2:size(A,2); 
    els = [els 'l'];
    fmtstr = [fmtstr '%c']; 
end; 
if ~(isempty(row)); 
    els = [els 'l'];
    fmtstr = [fmtstr '%c']; 
end; 
fprintf(['\\begin{tabular}{' fmtstr '}\n'],els); 
 

% display column headers
if ~isempty(head); 
    fmtstr = '%c';
    for i = 2:size(head,2);
        fmtstr = [fmtstr '%c']; 
    end; 
    fmtstr2 = fmtstr; 
    for i = 2:size(head,1); 
        fmtstr2 = [fmtstr2 ' & ' fmtstr]; 
    end;
    if ~isempty(row);
        fmtstr2 = ['     & ' fmtstr2 ' \\\\ \n']; 
    else; 
        fmtstr2 = [ fmtstr2 ' \\\\ \n'];
    end; 
    fprintf(fmtstr2 , head'); 
end; 

% display row headers and then rows

fmtstr2 = fmt; 
for i = 2:size(A,2); 
    fmtstr2 = [fmtstr2 ' & ' fmt]; 
end; 

if ~isempty(row); 
    rowfmtstr = '%c';
    for i = 2:size(row,2);
        rowfmtstr = [rowfmtstr '%c']; 
    end; 
end; 

for j=1:size(A,1)-1;
    if ~isempty(row);
        fprintf([rowfmtstr ' & '], row(j,:)); 
    end; 
    fprintf([fmtstr2 '\\\\ \n'], A(j,:));
end; 

if ~isempty(row);
        fprintf([rowfmtstr ' & '], row(end,:)); 
end; 
fprintf([fmtstr2 ' \n'], A(end,:));

fprintf('\\end{tabular} \n'); 



