function mb_save_json_file(object, filename)

% mb_save_json_file(object, filename)
%
% Convert a matlab object into a JSON string and save it to a file

fid = fopen(filename,'wt');
if fid < 0
  error('File %s cannot be written', filename);
end
fprintf(fid, '%s', jsonencode(object));
fclose(fid);
