spm_path = 'D:\spm8'; 
addpath(genpath(spm_path))

% gradientVector_path = '.\gradientVectors';
% addpath(genpath(gradientVector_path))

% mcc -v -m dke.m ...
%     -a gradient_vectors_siemens6.dat ...
%     -a gradient_vectors_siemens10.dat ...
%     -a gradient_vectors_siemens12.dat ...
%     -a gradient_vectors_siemens20.dat ...
%     -a gradient_vectors_siemens30.dat ...
%     -a gradient_vectors_siemens64.dat ...
%     -a gradient_vectors_siemens256.dat

% mcc -v -m dke.m -a ..\gradientVectors\*.dat
% 
% % dke_utils_compile
% mcc -m dke_preprocess_dicom.m -a spm_dicom_dict.mat
% mcc -m map_interpolate.m

fprintf('building DKE (win64, internal)\n')

fprintf('compiling DKE\n')
mcc -v -m dke.m -a ..\gradientVectors\*.dat

% dke_utils_compile
fprintf('compiling dke_preprocess_dicom\n')
mcc -v -m dke_preprocess_dicom.m -a spm_dicom_dict.mat

fprintf('compiling map_interpolate\n')
mcc -v -m map_interpolate.m

fprintf('completed\n\n')
