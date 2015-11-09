function script_file_path = create_slurm_script( data_file, output_folder, code_dir, k, Model )
%CREATE_SLURM_SCRIPT Summary of this function goes here
%   Detailed explanation goes here

matlab_cmd = ['/bin/matlab -nodesktop -nosplash -singleCompThread -logfile ' output_folder filesep 'logfile.txt'];
func_call = [' -r "cd Figure3_and_4/functions; runAlexsdata_prediction_func ' data_file ' ' output_folder ' ' code_dir ' ' num2str(k) ' ' Model '; exit;" '];

script_file_path = [output_folder filesep 'k' num2str(k) '_run' output_folder(end)];

fID = fopen(script_file_path,'w');
fprintf(fID, '#!/bin/bash\n');
hours = num2str(min(k,6)*4);
if length(hours) <2, hours = ['0' hours]; end;
fprintf(fID, ['#SBATCH --time=' hours ':00:00\n']);
fprintf(fID, '#SBATCH --qos=normal\n');
fprintf(fID, '#SBATCH --cpus-per-task=1\n');
fprintf(fID, ['cd ' code_dir '\n']);

% Find oldest matlab version on the current machine, and use that one
fprintf(fID,'cur_version=$(ls -d -t /opt/matlab-R201* | tail -1)\n');
fprintf(fID, ['$cur_version' matlab_cmd func_call '\n']);
fprintf(fID, 'echo $SLURM_JOB_ID');
fclose(fID);

end

