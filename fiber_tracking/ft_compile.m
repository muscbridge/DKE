spm_path = 'Y:\helpern_data\Programs\spm8_r4667'; 
addpath(genpath(spm_path))

fprintf('building ft (win64)\n')

fprintf('compiling ft\n')
mcc -v -m dke_ft.m

