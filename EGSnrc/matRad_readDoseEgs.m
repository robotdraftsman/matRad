function [bixelDose,bixelDoseError] = matRad_readDoseEgs(filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad binary dose import from vmc++
% 
% call
%   [bixelDose,bixelDoseError] = matRad_readDoseVmc(filename)
%
% input
%   filename:   path of input file
%
% output
%   bixelDose       = vector of imported dose values, [D]      = 10^-(10) Gy cm^2
%   bixelDoseError  = vector of imported dose errors, [deltaD] = 10^-(10) Gy cm^2
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%modify this to read the dos files created by Egs...so need that header
%information, even though the 3ddose to dos conversion only seems to put
%one thing in the header, the first one (no regions). So either sub the
%info in somehow else (basically  make it have a single-variable header and
%type up the other used variables here, i.e. header(5) -> dump_dose), or
%make the 3ddose to dos script make the full header

fid = fopen(filename,'r');

% read header (no regions, no histories, no batches, no beamlets, format specifier (dump_dose))
Header     = fread(fid,5,'int32');
no_regions = Header(1);
dump_dose  = 1; %Header(5);

% read dose array
if dump_dose == 2
    dmax            = fread(fid, 1, 'double');
    bixelDose       = fread(fid, no_regions, 'uint16');
    bixelDose       = bixelDose/65534*dmax; % conversion short integers to floating numbers
    bixelDoseError  = 0;
elseif dump_dose == 1
    bixelDose       = fread(fid, no_regions, 'float32');
    bixelDoseError  = fread(fid, no_regions, 'float32');
end

fclose(fid);
return;