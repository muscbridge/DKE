clear all

%move dicoms of interest to folder 'b12root/dicom_dke'

%subj_list={'ET_016'};
%{'ET_017'}
subj_list={'test'};	

for i=1:length(subj_list)
    b12root = ['/Users/lewis/prisma/' subj_list{i}];

    %convert dcm to nii
    
    eval(['!/Applications/MRIcroGL/dcm2niix -f %p ''' b12root '/dicom_dke'''])
    
    %move new niis
    
    %this part is needed when you have an additional separate b0 sequence
%     mkdir([b12root '/nifti/B0'])
%     in=fullfile(b12root, 'dicom_dke', '*B0*.nii');
%     out=fullfile(b12root, 'nifti/B0');
%     copyfile(in,out);

    mkdir([b12root '/nifti/DKI1'])
    file_in=dir(fullfile(b12root, 'dicom_dke', '*DKI*.nii'));
    in=fullfile(b12root, 'dicom_dke', file_in(1).name);
    out=fullfile(b12root, 'nifti/DKI1', '4D.nii');
    copyfile(in,out);
    
     %% Split 4D nii into 3D nii
        
        %for seperate b0
% Vdir = dir(fullfile(b12root, 'nifti/B0', '*.nii'));
% V=fullfile(b12root, 'nifti/B0', Vdir(1).name);
% Vo = spm_file_split(V,[b12root '/nifti/B0']);

Vdir = dir(fullfile(b12root, 'nifti/DKI1', '*.nii'));
V=fullfile(b12root, 'nifti/DKI1', Vdir(1).name);
Vo = spm_file_split(V,[b12root '/nifti/DKI1']);
    
    mkdir([b12root '/dke']);
    
    %% DENOISE
command=['/usr/local/mrtrix3/bin/dwidenoise ''' b12root '/nifti/DKI1/4D.nii'' ''' b12root '/nifti/DKI1/4D_DN.nii'' ' '-noise ''' b12root '/nifti/DKI1/noise.nii'''];
[status,cmdout] = system(command);

    %% UNRING

    DN=spm_read_vols(spm_vol(fullfile(b12root,'/nifti/DKI1','4D_DN.nii')));
    list=dir(fullfile(b12root,'/nifti/DKI1','*00*'));
    [dim1,dim2,dim3,dim4]=size(DN);
   parfor j=1:dim4
    img(:,:,:,j)=unring(DN(:,:,:,j));
    hdr=spm_vol(fullfile(b12root,'nifti/DKI1',[list(j).name]));
    hdr.dt=[16 0];
    int=img(:,:,:,j);
    int(isnan(int))=0;
    spm_write_vol(hdr,int);
   end
    
    %% Rename nii (3D_vol#_bval.nii)
    %may need to check indices
    
list=dir(fullfile(b12root,'/nifti/DKI1','*00*.nii'));
for k=2:65(list);
    list(k).bval = 'b1000';
end
for k=66:129(list); 
    list(k).bval = 'b2000';
end
for k=[1 130:138](list);
    list(k).bval = 'b0';
end

for k=1:length(list)
    [list(k).newname] = strrep(list(k).name,'4D', '3D');
    [list(k).newname] = strrep(list(k).newname, '.nii', ['_' list(k).bval '.nii']);
    movefile(fullfile(b12root,'/nifti/DKI1',list(k).name), fullfile(b12root,'/nifti/DKI1',list(k).newname));
end

%now move b0s to folder 'b12root/nifti/b0'
clear list
list=dir(fullfile(b12root,'/nifti/DKI1','*b0*.nii'));
mkdir([b12root '/nifti/B0']);
for l=1:length(list)
    in=fullfile(b12root, 'nifti/DKI1',list(l).name);
    out=fullfile(b12root, '/nifti/B0', list(l).name);
    movefile(in,out);
end
    
    %% make gradient file
    
    cd([b12root '/dicom_dke']);
    name=dir('*.bvec');
    A=importdata([name(1).name]);
    B=A(:,any(A));
    Gradient=B';
    Gradient1=Gradient(1:(round(end/2)),:);
    save(fullfile(b12root,'dke/gradient_dke.txt'),'Gradient1','-ASCII')
    
    %% create ft and dke params files
    
    %dke params
    fid=fopen('/Users/lewis/prisma/Scripts/dke_parameters.txt'); %Original file. CHANGE THIS
    fout=fullfile(b12root,'dke/dke_parameters.txt');% new file 

    fidout=fopen(fout,'w');

	while(~feof(fid))
    s=fgetl(fid);
    s=strrep(s,'PATH_REPLACE',b12root); %s=strrep(s,'A201', subject_list{i}) replace subject
    fprintf(fidout,'%s\n',s);
%    disp(s)
    end
    fclose(fid);
    fclose(fidout);
    
    %ft params
    fid=fopen('/Users/lewis/prisma/Scripts/ft_parameters.txt'); %Original file. CHANGE THIS
    fout=fullfile(b12root,'dke/ft_parameters.txt');% new file 

    fidout=fopen(fout,'w');

	while(~feof(fid))
    s=fgetl(fid);
    s=strrep(s,'PATH_REPLACE',b12root); %s=strrep(s,'A201', subject_list{i}) replace subject
    fprintf(fidout,'%s\n',s);
%    disp(s)
    end
    fclose(fid);
    fclose(fidout);
    

%% coregister b0s to b0
    %ONLY USE WITH INTERLEAVED B0S

% dirb0=dir(fullfile(b12root,'nifti/B0/*.nii'));
% dirdki=dir(fullfile(b12root,'nifti/DKI1/*.nii'));
% 
% fprintf('Co-registering images...\n')
% fn_source = fullfile(b12root, 'nifti/B0', dirb0(1).name); % source file is the first b = 0 image in the series returned by the operating system
% 
%             fn_target = fullfile(b12root, 'nifti/DKI1', dirdki(1).name);   
%             M=coregister(fn_target, fn_source, fullfile(b12root, 'nifti/B0'),'.nii');
% 
% delete(fullfile(b12root, 'nifti/DKI1', 'r*.nii'))
% movefile(fullfile(b12root, 'nifti/B0/r*.nii'),fullfile(b12root, 'nifti/B0_coreg')); 
% fprintf('Co-registration complete.\n')      

%% average b0's
list = dir(fullfile(b12root,'nifti/B0/*.nii'));
hdr = spm_vol(fullfile(b12root,'nifti/B0',list(1).name));

imgavg = spm_read_vols(hdr);
for j = 2:length(list)
    hdr = spm_vol(fullfile(b12root,'nifti/B0',list(j).name));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
end

imgavg = imgavg / (length(list));

mkdir(b12root, '/nifti/combined');
hdr.dt=[16 0];
hdr.fname = fullfile(b12root,'nifti/combined/B0_avg.nii');
imgavg(isnan(imgavg))=0;
spm_write_vol(hdr, imgavg);

%% move DKI & make 4D nii
copyfile([b12root '/nifti/DKI1/*00*'], [b12root '/nifti/combined']);
files = dir([b12root '/nifti/combined/*.nii']);
make_4D_nii([b12root '/nifti/combined'],{files.name},'4D.nii');
movefile([b12root '/nifti/combined/4D.nii'],[b12root '/dke/4D.nii'])

 img=spm_read_vols(spm_vol([b12root '/dke/4D.nii']));
 img(isnan(img))=0;
 make_4D_nii(spm_vol([b12root '/dke/4D.nii']),img,'4D.nii');

%% run DKE and DKE_FT

    fn_params=fullfile(b12root, 'dke/dke_parameters.txt')
    dke(fn_params)

    FT_parameters=fullfile(b12root, 'dke/ft_parameters.txt')
    dke_ft(FT_parameters)

end