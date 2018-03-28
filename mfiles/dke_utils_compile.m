spm_path = 'D:\spm8';
addpath(genpath(spm_path))

mcc -v -m dke_preprocess_dicom.m ...
    -a spm_dicom_dict.mat

mcc -v -m map_interpolate.m
