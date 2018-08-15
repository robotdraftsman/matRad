%write phase space file for each beamlet
function headers = writeBinaryPHSP_optimized(phspData,charges,header)

%total number of particles through this beamlet:
%beamlet boundaries in x and y are (stf.ray.rayPos_bev(x or y) +/- 2.5)/2
%-> in mm. divide by 20 instead for cm. And y mentioned here is actually z
%in stf

%test out a beamlet in first beam:
%gets all the points in a given beamlet:
phspPath = 'beamletPHSPfiles';
filebase = 'dividedPhsp36M';

xmin = -45; %these are in mm
xmax = 45;
ymin = -45;
ymax = 45;

allBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);
n = 1;
for x = xmin:5:xmax
    for y = ymin:5:ymax
        allBeamlets(n,:) = [x y];
        n = n + 1;
    end
end

"arranged the beamlets...\n"

%now loop over these and cut up the phsp file
%scatter(allTheStuff(3:7:length(phspData.Data.allTheStuff))*10,phspData.Data.allTheStuff(4:7:length(phspData.Data.allTheStuff))*10)
%hold on

%one row for each beamlet. Col1 = numParticles, col2 = numPhotons:
headerNumbers = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);

%one row for each beamlet. col1 = maxEnergy, col2 = min elec energy in beamlet:
%set all entries to inf so that don't have to worry about 0s messing up min
%energy values -> any particle's energy is less. So on first iteration it
%becomes a finite number: the first particle in that beamlet's energy
maxEnergy = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),1);
minEnergy = inf(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),1);

%new structure: loop through phspData and check which beamlet each particle
%is in. Can either accumulate stuff to write to each file or just write
%directly to each file. In latter case, would have an array of fids and
%just access the relevant one based on the x and y of the particle we're
%looking at

%create array of fids for all the beamlets:
example_txtmsg('PHSP file dividing process:','opening the smol beamlet phsp files...');
filenum = 1;
for x = xmin:5:xmax
    for y = ymin:5:ymax
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
tic;
%loop over all the particles in the array:
example_txtmsg('PHSP file dividing process:','looping over the particles in the big file: populating beamlets...');
for i = 1:7:length(phspData.Data.allTheStuff)
    
    %see which beamlet it's within (boundaries included):
    I = find( ( (double(phspData.Data.allTheStuff(i+2))*10 - allBeamlets(:,1)/2 >= -2.5/2) & (double(phspData.Data.allTheStuff(i+2))*10 - allBeamlets(:,1)/2 <= 2.5/2) ) & ( (double(phspData.Data.allTheStuff(i+3))*10 - allBeamlets(:,2)/2 >= -2.5/2) & (double(phspData.Data.allTheStuff(i+3))*10 - allBeamlets(:,2)/2 <= 2.5/2)) );
    %check if it shares a border, i.e. I > 1 (and must be < 4):
    if (length(I)>1)
        inBeamletNum = I(end);
    else
        inBeamletNum = I;
    end
    
    if (inBeamletNum ~= 0)
        headerNumbers(inBeamletNum,1) = headerNumbers(inBeamletNum,1)+1;
        if(charges( (i-1)/7 + 1 ) == 0)
           headerNumbers(inBeamletNum,2) = headerNumbers(inBeamletNum,2) + 1;
        end
        fwrite(fids(inBeamletNum),phspData.Data.allTheStuff(i:i+6),'single');

        if(maxEnergy(inBeamletNum) < double(phspData.Data.allTheStuff(i+1)))
            maxEnergy(inBeamletNum) = double(phspData.Data.allTheStuff(i+1));
        end
        %check if it's an electron and its energy is smaller than min elec
        %energy for the given beamlet:
        if( (charges( (i-1)/7 + 1 ) == -1) & (minEnergy(inBeamletNum) > double(phspData.Data.allTheStuff(i+1))))
            minEnergy(inBeamletNum) = double(phspData.Data.allTheStuff(i+1));
        end
    end
    if (mod((i+6)/7,round(header.Data.NUM_PHSP_TOT/10))==0)
        didSuch = strcat('done=',num2str(((i+6)/7)/m.Data.NUM_PHSP_TOT),'% of the particles!');
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
