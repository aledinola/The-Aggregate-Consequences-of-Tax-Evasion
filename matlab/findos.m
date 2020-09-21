function [os] = findos()

%This function finds the OS of the computer on which the program runs.
%This is needed for the usage of the FORTRAN executables.

os_raw = computer('arch');

switch os_raw
     case 'win32'
		os = 'windows';
     case 'win64'
		os = 'windows';
     case 'glnx86'
		os = 'unix';
     case 'glnxa64'
		os = 'unix64';
     case 'maci64'
		os = 'mac';
end

end %END FUNCTION "findos"
