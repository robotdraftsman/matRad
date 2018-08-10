%read in a binary file (the phase space file)
function [m m2 charges lastParticle numParticlesLeft] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip)
%clear;

% read the header:
 m = memmapfile(phspFile,...
     'Format',{'int8',[1,5],'mode';...
        'int32',[1,1],'NUM_PHSP_TOT';...
        'int32',[1,1],'PHOT_PHSP_TOT';...
        'single',[1,1],'EKMAX_PHSP_SHORT';...
        'single',[1,1],'EKMINE_PHSP_SHORT';...
        'single',[1,1],'NINC_PHSP_SHORT'},'Repeat',1);

%if want to read it all in at once:
%readThisMuch = double(m.Data.NUM_PHSP_TOT);

readOffset = 28 + numParticlesToSkip*28;
%when running this 1M particles at a time:
% numParticlesToSkip = (i-1)*1000000, with i = iteration we're on in the
% loop we're dealing with in this other file we're calling this from

if(readThisMuch > double(m.Data.NUM_PHSP_TOT))
   readThisMuch =  double(m.Data.NUM_PHSP_TOT);
   "a"
end
if(readThisMuch + numParticlesToSkip > double(m.Data.NUM_PHSP_TOT))
   readThisMuch = double(m.Data.NUM_PHSP_TOT)-numParticlesToSkip;
   "b"
end
numParticlesLeft = (numParticlesToSkip + readThisMuch)
if(numParticlesToSkip >= double(m.Data.NUM_PHSP_TOT))
   fprintf("You're finished with the file!\n");
   numParticlesLeft = 0;
   phspData = [0];
   lastParticle = -1;  %interpret this as meaning it's done
   "c"
   quit;
end

%read everything after the header:
%but whether we read in zlast or not depends on the mode:
if (m.Data.mode == 'MODE0')
    "Running in MODE0: no ZLAST"
    %store just the LATCH variable for easy extraction of the charge:
    m3 = memmapfile(phspFile,'Offset',readOffset,...
        'Format',{'uint8',[readThisMuch*28 1],'LATCHandEverythingElse'},'Repeat',1);
    
    %now to get the required information out of latch
    charges = zeros(length(readThisMuch),1);
    %for complement representation of signed integers:
    %(would probbaly take around half an hour for all 36 million particles)
    count = 1;
    for i = 1:28:length(m3.Data.LATCHandEverythingElse)
        getBits = dec2bin(m3.Data.LATCHandEverythingElse(i));
        abyte = '00000000';
        abyte(9-length(getBits):8) = getBits(:);
        if(abyte(2) == '1' && abyte(3) == '0')
            charge = -1;
        elseif (abyte(2) == '0' && abyte(3) == '1')
            charge = 1;
        elseif (abyte(2) == '0' && abyte(3) == '0')
            charge = 0;
        else
            charge = 42;    %just so it's easily found
            fprintf("You have a 11 in binary for charge, dude!");
        end
        charges(count) = charge;
        count = count + 1;
    end
    clear m3;
    

    fprintf("Got the charges. Now the rest...\n");

    m2 = memmapfile(phspFile,'Offset',readOffset,...
        'Format',{'single', [readThisMuch*7 1], 'allTheStuff'},'Repeat',1);
    
elseif (m.Data.mode == 'MODE2')
    fprintf("Sorry, I only take phsp files without ZLAST (only run in MODE0).\n");
else
    "Uhh it's not mode0 or mode2 dude??\n"
end

strMode = char(m.Data.mode);
lastParticle = numParticlesToSkip + readThisMuch; %later use as numParticlesToSkip

end