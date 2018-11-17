function headers = writeBinaryPHSP_optimized(phspData,charges,header)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% write phase space file for each bixel in a hypothetical beam that
% contains all possible bixels/beamlets
% 
% phspData:     An array containing everything read in from a phase space
%               file get it in the right format from readBinaryPHSP_optimized.m
% charges:      An arry of the charge of every particle in that phsp file
% header:       The header of the large phsp file
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phspPath = 'EGSnrc/beamletPHSPfiles';
filebase = 'dividedPhsp36M';

bixelSize = 5; %in mm

%generate the beamlets. If you have a bigger beam with more bixels or
%bigger bixels, you might want to change this...
xmin = -45; %these are in mm. Get these from the matRad steering file
xmax = 45;
ymin = -45;
ymax = 45;

allBeamlets = zeros(((xmax-xmin)/bixelSize + 1 )*( (ymax-ymin)/bixelSize + 1),2);
n = 1;
for x = xmin:bixelSize:xmax
    for y = ymin:bixelSize:ymax
        allBeamlets(n,:) = [x y];
        n = n + 1;
    end
end

"arranged the beamlets...\n"

%now loop over the above beamlets and cut up the phsp file:

%one row for each beamlet. Col1 = numParticles, col2 = numPhotons:
headerNumbers = zeros(((xmax-xmin)/bixelSize + 1 )*( (ymax-ymin)/bixelSize + 1),2);

%one row for each beamlet. col1 = maxEnergy, col2 = min elec energy in beamlet:
%set all entries to inf so that don't have to worry about 0s messing up min
%energy values -> any particle's energy is less. So on first iteration it
%becomes a finite number: the first particle in that beamlet's energy.
maxEnergy = zeros(((xmax-xmin)/bixelSize + 1 )*( (ymax-ymin)/bixelSize + 1),1);
minEnergy = inf(((xmax-xmin)/bixelSize + 1 )*( (ymax-ymin)/bixelSize + 1),1);


%create array of fids for all the beamlets:
example_txtmsg('PHSP file dividing process:','opening the small beamlet phsp files...');
filenum = 1;
for x = xmin:bixelSize:xmax
    for y = ymin:bixelSize:ymax
        phspFile = fullfile(phspPath, strcat(filebase,num2str(filenum),'.egsphsp1'));
        fids(filenum) = fopen(phspFile,'w+');
        filenum = filenum + 1; 
    end
end

%placeholder for the header:
for i = 1:(filenum-1)
   fwrite(fids(i),[0 0 0 0 0 0 0],'int32'); 
end

"opened the files, made placeholder headers..."



%loop over all the particles in the big phase space file array and assign each particle to a beamlet

example_txtmsg('PHSP file dividing process:','looping over the particles in the big file: populating beamlets...');
notifyIncrement = 3760000;  %roughly 10% of the particles in the phsp I was using; change if you want

tic;
for i = 1:7:length(phspData.Data.allTheStuff)
    
    %see which beamlet it's within (boundaries included):
    I = find( ( (double(phspData.Data.allTheStuff(i+2))*10 - allBeamlets(:,1)/2 >= -2.5/2) & (double(phspData.Data.allTheStuff(i+2))*10 - allBeamlets(:,1)/2 <= 2.5/2) ) & ( (double(phspData.Data.allTheStuff(i+3))*10 - allBeamlets(:,2)/2 >= -2.5/2) & (double(phspData.Data.allTheStuff(i+3))*10 - allBeamlets(:,2)/2 <= 2.5/2)) );
    %check if it shares a border, i.e. I > 1 (and must be < 4):
    if (length(I)>1)
        inBeamletNum = I(end);
    elseif (length(I) == 1)
        inBeamletNum = I;
    else
    	inBeamletNum = 0;
    end
    
    %if this particle is in one of the beamlets, do all the stuff:
    if (inBeamletNum ~= 0)
        headerNumbers(inBeamletNum,1) = headerNumbers(inBeamletNum,1)+1;
        if(charges( (i-1)/7 + 1 ) == 0)
           headerNumbers(inBeamletNum,2) = headerNumbers(inBeamletNum,2) + 1;
        end
        fwrite(fids(inBeamletNum),phspData.Data.allTheStuff(i),'single');   %write LATCH
        fwrite(fids(inBeamletNum),abs(phspData.Data.allTheStuff(i+1)),'single');    %write energy (but we don't want the minus signs so just using abs value)
        fwrite(fids(inBeamletNum),phspData.Data.allTheStuff(i+2:i+6),'single'); %write the rest (x,y,u,v,wt)

        if(maxEnergy(inBeamletNum) < abs(double(phspData.Data.allTheStuff(i+1))))
            maxEnergy(inBeamletNum) = abs(double(phspData.Data.allTheStuff(i+1)));
        end
        %check if it's an electron and its energy is smaller than min elec
        %energy for the given beamlet:
        if( (charges( (i-1)/7 + 1 ) == -1) & (minEnergy(inBeamletNum) > abs(double(phspData.Data.allTheStuff(i+1)))))
            minEnergy(inBeamletNum) = abs(double(phspData.Data.allTheStuff(i+1)));
        end
    end
    if (mod((i+6)/7,notifyIncrement)==0)
        didSuch = strcat('percent of particles assigned to beamlets = ',num2str((double((i+6)/7)/double(header.Data.NUM_PHSP_TOT))*100));
        example_txtmsg('PHSP file dividing process:',didSuch);
    end
end
toc;

example_txtmsg('PHSP file dividing process:','Done looping over the particles. :) now making the headers...');
"finished looping over all those many particles!"

clear allBeamlets;


tic;
%go back through the files and write their headers:
for(i = 1:length(fids))
    fseek(fids(i),0,'bof');
    fwrite(fids(i), header.Data.mode, 'int8');
    fwrite(fids(i), headerNumbers(i,1),'int32'); %write the number of particles
    fwrite(fids(i), headerNumbers(i,2),'int32'); %number of photons
    fwrite(fids(i), maxEnergy(i),'single');    %max kinetic energy
    %min elec kinetic energy (don't want places w 0 electrons to be inf):
    if(minEnergy(i) == inf)
        minEnergy(i) = 0;
    end
    fwrite(fids(i), minEnergy(i),'single');
    fwrite(fids(i), header.Data.NINC_PHSP_SHORT,'single');
    fwrite(fids(i),[0 0 0],'int8');
end
toc;
example_txtmsg('PHSP file dividing process:','made the headers. We are done here!')
"wrote the headers!! Alll donee!\n"

headers(:,1) = headerNumbers(:,1);
headers(:,2) = headerNumbers(:,2);
headers(:,3) = maxEnergy(:);
headers(:,4) = minEnergy(:);

fclose all;
clear fids;


% [communist geologists]
% Ivan: look at all this mica
% Igor: no comrade, it's ourca.

end
