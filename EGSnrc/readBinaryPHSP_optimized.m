function [m m2 charges lastParticle numParticlesLeft] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% read in a binary file (the phase space file)
% 
% Output arguments:
% 
% m:                    Header information from the phase space file
% m2:                   Particle data from phsp file (in one long array)
% charges:              The charges of these particles
% lastParticle:         Last particle read in, if only part of file read
% numParticlesLeft:     Particles yet to be read in
% 
% Note on the last two: these are not used currently, but are here so that,
% if the phsp file is too big to be read in at once, it can go a bit at a
% time. I didn't finish developing this methodology
% 
% Input arguments:
% 
% phspFile:             The phase space file's name
% readThisMuch:         The number of particles to read in (see note above)
% numParticlesToSkip:   If reading in only part of the file, and part has
%                       already been read in, skip forward by this many particles
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
example_txtmsg('Large PHSP file read-in:','Read the header. Now check that user inputs make sense...');

readOffset = 28 + numParticlesToSkip*28;
% if running this 1M particles at a time:
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
   example_txtmsg('Large PHSP file read-in:','Inputs said to skip way wast all the particles...no can do, buddy. THIS ENDS HERE >:(');
   error("you're telling me to skip past all the particles, and then some...");
end

%read everything after the header:
%but whether we read in zlast or not depends on the mode:
if (m.Data.mode == 'MODE0')
    "Running in MODE0: no ZLAST"
    %store just the LATCH variable for easy extraction of the charge:
    m3 = memmapfile(phspFile,'Offset',readOffset,...
        'Format',{'uint8',[readThisMuch*28 1],'LATCHandEverythingElse'},'Repeat',1);
    
    example_txtmsg('Large PHSP file read-in:','Extracting the charges...');
    
    %now to get the required information out of latch
    charges = inf(length(readThisMuch),1);
    %for complement representation of signed integers:
    tic;
    count = 1;
    notifyIncrement = 3760000;
    for i = 4:28:length(m3.Data.LATCHandEverythingElse)
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
        if (mod(count,notifyIncrement)==0)
            didSuch = strcat('percent of charges extracted = ',num2str((double(count)/double(m.Data.NUM_PHSP_TOT))*100));
            example_txtmsg('Large PHSP file reading process:',didSuch);
        end
    end
    clear m3;
    toc;

    example_txtmsg('Large PHSP file read-in:','Got the charges. Now the rest...');
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

example_txtmsg('Large PHSP file read-in:','Finished reading in the file - all done with this function :)');

end