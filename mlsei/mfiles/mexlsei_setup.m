% This file contains the procedure for building the mexlsei function from
% scratch. Not every step in this file may be necessary as some will likely
% have been performed already, such as steps 2 and 3. 
% 

%% 1. Change to the installation directory (where this file is)
% you may need to do this manually if mexlsei_setup is not on the path

rootdir = fullfile(fileparts(which('mexlsei_setup')));
cd(rootdir);

%% 2. Build the f2c lib which must be included

temp = pwd;
cd('..');

cd(fullfile(pwd, 'c', 'libf2c'));


% get the local machine architecture
machine = computer('arch');

if strcmp(machine,'win32')

    % To compile use Visual C++ Express Edition 2008
    % Open the visual studio 2008 Command prompt and run the
    % following
    %
    % > nmake -f makefile.vc64
    libfilename = 'vcf2c.lib';
    fprintf(1, 'Manual steps are required for win32');
    break;
    
elseif strcmp(machine, 'win64')
    
    % To compile use Visual C++ Express Edition 2008
    % Open the visual studio Cross-Tools x64 Command prompt and run the
    % following
    %
    % > nmake -f makefile.vc64
    libfilename = 'vc64f2c.lib';
    fprintf(1, 'Manual steps are required for win64');
    break;
    
elseif strcmp(machine,'glnxa64')
    
    % compile using gcc on 64 bit platform
    !make -f makefile.lx64
    libfilename = 'libf2cx64.a';
    
elseif strcmp(machine,'glnx86')
    
    % compile using gcc
    !make -f makefile.lx86
    libfilename = 'libf2cx86.a';
    
end

% Move the library file to the directory above
movefile(libfilename, fullfile(rootdir, 'c'));

cd(temp);

%% 3. Convert the fortran code files into a single C file dlsei.c using
%% f2c.exe

if ispc
    
    !type .\fortran\*.f | ".\c\f2c.exe" -A -R > .\c\dlsei.c

else
    
    !cat ./fortran/*.f | "./c/f2c.exe" -A -R > .\c\dlsei.c
    
end

%% 4. Now mex the files, including the appropriate lib

machine = computer('arch');

if strcmp(machine,'win32')

    % To compile use Visual C++ Express Edition 2008
    % I had to move the library file (vcf2c.lib) to a directory with no
    % spaces, 'C:\libraries' on my system, you will probably need to change
    % this directory for your windows machine
    mex .\c\mexlsei.c .\c\dlsei.c -outdir . -L"C:\libraries" -lvcf2c LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:LIBCMT.lib" -v
    
elseif strcmp(machine,'win64')
    
	% To compile use Visual C++ Express Edition 2008
    % I had to move the library file (vc64f2c.lib) to a directory with no
    % spaces, 'C:\libraries' on my system, you will probably need to change
    % this directory for your windows machine
    mex .\c\mexlsei.c .\c\dlsei.c -outdir . -L"C:\libraries" -lvc64f2c LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:LIBCMT.lib" -v
    
elseif strcmp(machine,'glnxa64')
    
    % compile using gcc
    mex ./c/mexlsei.c ./c/dlsei.c -outdir . -L"/home/s0237326/Postgrad_Research/fortran/mexlsei/c" -lf2cx64 -v

elseif strcmp(machine,'glnx86')
    
    % compile using gcc
    mex ./c/mexlsei.c ./c/dlsei.c -outdir . -L"/home/s0237326/Postgrad_Research/fortran/mexlsei/c" -lf2cx86 -v
    
end


%%

load Test_lsei.mat

copyfile('~/Postgrad_Research/fortran/mexlsei/mexlsei.mexw64', '/home/s0237326/Postgrad_Research/MATLAB_Scripts/subversion/matlab/Useful_Functions/mlsei/')

[x, rnorme, rnorml, mode] = mlsei(A, b, Mineq, rhsineq, Meq, rhseq);

