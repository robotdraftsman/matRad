%write phase space file for each beamlet

%% load stuff so I don't have to run matRad every time:

load TG119.mat
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [0:72:359]; % [?]
pln.propStf.couchAngles     = [0 0 0 0 0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);


% optimization settings
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% now do the actual stuff

%total number of particles through this beamlet:
%beamlet boundaries in x and y are (stf.ray.rayPos_bev(x or y) +/- 2.5)/2
%-> in mm. divide by 20 instead for cm. And y mentioned here is actually z
%in stf

%later on make it a loop over all beamlets in a beam, but for now just make
%it work with the first one:

%put things in a matrix
%tic;
phspSrc = zeros(length(m2.Data),7+3); %if in mode2 (include zlast) 7->8
%the +3 part is because I'm storing each latch byte separately (the other
%additional spaces are not used)

for i = 1:length(m2.Data)
    phspSrc(i,1:4) = m2.Data(i).LATCH;
    phspSrc(i,5) = m2.Data(i).ESHORT;
    phspSrc(i,6) = m2.Data(i).X_PHSP_SHORT;
    phspSrc(i,7) = m2.Data(i).Y_PHSP_SHORT;
    phspSrc(i,8) = m2.Data(i).U_PHSP_SHORT;
    phspSrc(i,9) = m2.Data(i).V_PHSP_SHORT;
    phspSrc(i,10) = m2.Data(i).WT_PHSP_SHORT;
end
%toc; %so for 10^6 particles would take ~ 20 minutes :/ but can't think of a way to speed it up

%test out a beamlet in first beam:
%gets all the points in a given beamlet:
%n = 156;   %when looping over beamlets n = beamlet number. rn 156 (in centre of beam so overpapped w/ mini phsp file)
phspPath = 'beamletPHSPfiles';
filebase = 'dividedPhsp';

xmin = -40; %these are in mm
xmax = 40;
ymin = -45;
ymax = 45;

allBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);

%now loop over these and cut up the phsp file as I did before
%scatter(phspSrc(:,6)*10,phspSrc(:,7)*10)
%hold on
n = 1;

numParticles = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),1);

for x = xmin:5:xmax
    for y = ymin:5:ymax
        if(x == xmax && y == ymax)
            I = find( ( (phspSrc(:,6)*10 - x/2) >= -2.5/2 & (phspSrc(:,6)*10 - x/2) <= 2.5/2 ) & ( (phspSrc(:,7)*10 - y/2) >= -2.5/2 & (phspSrc(:,7)*10 - y/2) <= 2.5/2 ) ); 
        elseif(x==xmax)
            I = find( ( (phspSrc(:,6)*10 - x/2) >= -2.5/2 & (phspSrc(:,6)*10 - x/2) <= 2.5/2 ) & ( (phspSrc(:,7)*10 - y/2) >= -2.5/2 & (phspSrc(:,7)*10 - y/2) < 2.5/2 ) ); 
        elseif(y==ymax)
            I = find( ( (phspSrc(:,6)*10 - x/2) >= -2.5/2 & (phspSrc(:,6)*10 - x/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - y/2) >= -2.5/2 & (phspSrc(:,7)*10 - y/2) <= 2.5/2 ) ); 
        else
            I = find( ( (phspSrc(:,6)*10 - x/2) >= -2.5/2 & (phspSrc(:,6)*10 - x/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - y/2) >= -2.5/2 & (phspSrc(:,7)*10 - y/2) < 2.5/2 ) ); 
        end
        
        allBeamlets(n,:) = [x y];
        
        phspFile = fullfile(phspPath, strcat(filebase,num2str(n),'.egsphsp1'));
        fid = fopen(phspFile,'w');

        % energies, charges, and electron energies
        importantStuffs = zeros(length(I),3);
        eindex = 1;
        for i = 1:length(I)
           %need charges, max energies
           importantStuffs(i,1) = phspSrc(I(i),5);  %energies
           importantStuffs(i,2) = charges(I(i)); %charges
           
           if(importantStuffs(i,2) == -1)
              importantStuffs(eindex,3) = importantStuffs(i,1); 
              eindex = eindex + 1;
           end
        end


        fwrite(fid, m.Data.mode, 'int8');
        fwrite(fid, length(I),'int32'); %write the number of particles
        fwrite(fid, length(importantStuffs(:,2) == 0),'int32'); %number of photons
        fwrite(fid, max(abs(importantStuffs(:,1))),'single');    %max kinetic energy
        
        % now get the min electron energy:
        fwrite(fid, min(abs(importantStuffs(1:eindex,3))),'single');    %min kinetic energy

        fwrite(fid, m.Data.NINC_PHSP_SHORT,'single');
        
        numParticles(n) = length(I);
        
        
        %header finished, now write those weird 3 bytes that are formatted in:
        fwrite(fid, [0 0 0],'uint8');

        %And write the rest of the stuff
        for i = 1:length(I)
            fwrite(fid, phspSrc(I(i),1), 'uint8');
            fwrite(fid, phspSrc(I(i),2), 'uint8');
            fwrite(fid, phspSrc(I(i),3),'uint8');
            fwrite(fid, phspSrc(I(i),4),'uint8');
            fwrite(fid, phspSrc(I(i),5),'single');
            fwrite(fid, phspSrc(I(i),6),'single');
            fwrite(fid, phspSrc(I(i),7),'single');
            fwrite(fid, phspSrc(I(i),8),'single');
            fwrite(fid, phspSrc(I(i),9),'single');
            fwrite(fid, phspSrc(I(i),10),'single');
        end
        scatter( [(x+2.5)/2 (x-2.5)/2 (x+2.5)/2 (x-2.5)/2],[ (y + 2.5)/2 (y - 2.5)/2 (y - 2.5)/2 (y + 2.5)/2])
        fclose(fid);
        n = n + 1;
    end
end

%in a holistic perspective, probably the best naming conventioin for the
%file names would just be *1, *2, *3,... because then it's easier to
%automate the dosxyz running stuff and the phasespace files being stitched
%together, and then here in matlab I can have some sort of checker to use
%with the beamlets, or an array that matches numbers to x and y coordinates
%- the reference array I can make even while making the beamlets phsp files

%Now for the specific beams:
%For a given beam, am looping over beamlets that it has. I'm given the
%centre of the beamlet. Can then do a find(beamlet centre == phsp array)
%which'll give the part in the phsp comparison array (need to make this !)
%that

%So loop over beamlets. For each, see which index it is in the checker
%array = number in filename. Collect all these indices in an array, and we
%know which files we want.
%So what do I do with this information? Once I move this over to the
%cluster, I can make it run dosxyz on each of these files (a loop that
%loads parameters, loads file (do the stuff from the command line, except
%write it up in here and stick commands into comand line using SYSTEM)) in
%a loop that submits one at a time to be queued to the cluster. Have the
%phsp files sitting in the cluster, so then this list of which files are
%included in a particular beam is then used only to ---
%OH WAIT - I can go and run dosxyz on all these phsp files made from
%initially (pre-beams) cutting up the phsp file. Then I have the .3ddose for
%each of these. Then I just use the indices when stitching together the
%doses for each beam - then save as *beam1,2,3,.... To run it I'll use
%SYSTEM to call the c program (add_...). I can modify it more if I'm unable
%to make matlab run concurrently w/ c -> bc I need to be able to send the
%file name base and number of files to the program while it's running. Not
%100% sure what I'll do there, but maybe I might have to make a shell
%script to do things. Or make it run in bk (background).


%Do the loop over beamlets here for now. Store info in a matrix (say each
%column is a beam, and each element in the column is --

%Oh actually two ways to do this: construct final matrix of what's what
%from 1s and 0s so that 1 in an index means that that beamlet in the
%rectangle with that index is included and 0 is excluded. Or can have a
%list of the actual numbers/indices.
%With binary, loop through column (when making final doses) and if(==1)
%then you include the file of the index you're on.
%With other, you loop through column and include the file of number
%*entryValue.3ddose.
%I think the binary one is more genralizable


%So this makes an array where each column is for a different beam, and each
%entry in the column is a binary saying whether the beam does/doesn't
%include the beamlet corresponding to its index (and hence also whether the
%phsp file for that should be used to construct the dose).
whichBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1)  ,length(stf));  %beamlets included (index in a given row) for each beam (columns)
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        indexOfPhsp = find( (allBeamlets(:,1) == stf(i).ray(j).rayPos_bev(1)) & (allBeamlets(:,2) == stf(i).ray(j).rayPos_bev(3)) );    %find index in reference vector that corresponds to this beamlet
        whichBeamlets(indexOfPhsp,i) = 1;
    end
end




% =======================================================================
% So the following can be picked at and later modified to suit our needs:
% =======================================================================


%loop over beams ==> now just need to change stf(1) to stf(j) everywhere
%tic;
%for j = 1:length(stf)
    scatter(phspSrc(:,6)*10,phspSrc(:,7)*10)
    hold on
%     j = 1;  %so I only have to look at beam 1 right now
%     tic;
%     %loop over beamlets/rays in beam
%     for n = 1:stf(j).numOfRays
% 
%         %if(n == 136 || n == 137 || n == 138 || n == 155 || n == 156 || n == 157 || n == 174 || n == 175 || n == 176)
%         
%         %for the following: how to detect if it's at an edge or corner?
%         %Seeing as we just have one long vector of rayPos points, it's not
%         %an array in the shape of the beam, must get creative...
%         %goes along y for a given x, then goes to next x and goes along y
%         %and so on. The moment x changes the we know we're at max y. The
%         %moment a given y is < previous y, we know we've changed x (but I
%         %don't think this info is useful - want end x)
%         %Maybe do a find() on it, for points that would satisfy x or y at
%         %max for a given row/column, then say if (index in loop = oneof
%         %these indices) then do the inclusion. Have three sets of indices
%         %then, for each case: x max, y max, and x,y max.
%         
%         %
%         
%         
%         %tic;
% %         if(is nothing to its right _and_ above it)
% %             I = find( ( (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) <= 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) <= 2.5/2 ) );
% %         elseif(is nothing above it)
% %             I = find( ( (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) <= 2.5/2 ) );
% %         elseif(is nothing to its right)
% %             I = find( ( (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) <= 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) < 2.5/2 ) );
% %         else
%             %is just inside the beam, not on the edges
%             I = find( ( (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(j).ray(n).rayPos_bev(1)/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(j).ray(n).rayPos_bev(3)/2) < 2.5/2 ) );
%         %end
%         %toc;    %would take ~10.5 minutes per beamlet for 10^7 particles
% 
%         %I gives the indices in phspSrc corresponding to the contained particles.
%         %Loop through this and I can extract the stuff to be written to the new
%         %binary phsp file
%         phspFile = fullfile(phspPath, strcat('beamletTestphsp',num2str(n),'.egsphsp1'));
%         fid = fopen(phspFile,'w');
% 
%         importantStuffs = zeros(length(I),2);
%         for i = 1:length(I)
%            %need charges, max energies
%            importantStuffs(i,1) = phspSrc(I(i),5);  %energies
%            importantStuffs(i,2) = charges(I(i)); %charges
%         end
% 
%         [emax,imax] = max(abs(importantStuffs(:,1)));
%         [emin,imin] = min(abs(importantStuffs(:,1)));
% 
%         fwrite(fid, m.Data.mode, 'int8');
%         fwrite(fid, length(I),'int32'); %write the number of particles
%         fwrite(fid, length(importantStuffs(:,2) == 0),'int32'); %number of photons
%         fwrite(fid, importantStuffs(imax),'single');    %max kinetic energy
%         fwrite(fid, importantStuffs(imin),'single');    %min kinetic energy
%         fwrite(fid, m.Data.NINC_PHSP_SHORT,'single');
% 
%         %header finished, now write those weird 3 bytes that are formatted in:
%         fwrite(fid, [0 0 0],'uint8');
% 
%         %And write the rest of the stuff
%         for i = 1:length(I)
%             fwrite(fid, phspSrc(I(i),1), 'uint8');
%             fwrite(fid, phspSrc(I(i),2), 'uint8');
%             fwrite(fid, phspSrc(I(i),3),'uint8');
%             fwrite(fid, phspSrc(I(i),4),'uint8');
%             fwrite(fid, phspSrc(I(i),5),'single');
%             fwrite(fid, phspSrc(I(i),6),'single');
%             fwrite(fid, phspSrc(I(i),7),'single');
%             fwrite(fid, phspSrc(I(i),8),'single');
%             fwrite(fid, phspSrc(I(i),9),'single');
%             fwrite(fid, phspSrc(I(i),10),'single');
%         end
%         fclose(fid);
%         %end
%         scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% 
%     end
%     toc;
%end
%toc;

% 
% n = 155;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 156;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 157;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 136;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 137;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 138;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 174;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 175;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
% n = 176;
% scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])


%relevant beamlets are numbers 155, 156, 157, 136, 137, 138, 174, 175, 176

%hold off

clear m;
clear m2;

% [communist geologists]
% Ivan: look at all this mica
% Igor: no comrade, it's ourca.
