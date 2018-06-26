%write phase space file for each beamlet

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
tic;
for n = 1:stf(1).numOfRays

    %if(n == 136 || n == 137 || n == 138 || n == 155 || n == 156 || n == 157 || n == 174 || n == 175 || n == 176)
    %tic;
    I = find( ( (phspSrc(:,6)*10 - stf(1).ray(n).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(1).ray(n).rayPos_bev(1)/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(1).ray(n).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(1).ray(n).rayPos_bev(3)/2) < 2.5/2 ) );
    %toc;    %would take ~10.5 minutes per beamlet for 10^7 particles

    %I gives the indices in phspSrc corresponding to the contained particles.
    %Loop through this and I can extract the stuff to be written to the new
    %binary phsp file
    phspFile = fullfile(phspPath, strcat('beamletTestphsp',num2str(n),'.egsphsp1'));
    fid = fopen(phspFile,'w');

    %put header stuff here:
    %for that I need all the stuff from prev header as it was EXCEPT num photons
    %and num particles: those will be based on what's in the square
    %or WAIT I would want to change E_kin min and max, no? According to max and
    %min in the beamlet, not in the space as a whole.

    %get min and max kinetic energy of the contained particles: want to have
    %array of just the contained particles then?

    importantStuffs = zeros(length(I),2);

    for i = 1:length(I)
       %need charges, max energies
       importantStuffs(i,1) = phspSrc(I(i),5);  %energies
       importantStuffs(i,2) = charges(I(i)); %charges
    end

    [emax,imax] = max(abs(importantStuffs(:,1)));
    [emin,imin] = min(abs(importantStuffs(:,1)));

    fwrite(fid, m.Data.mode, 'int8');
    fwrite(fid, length(I),'int32'); %write the number of particles
    fwrite(fid, length(importantStuffs(:,2) == 0),'int32'); %number of photons
    fwrite(fid, importantStuffs(imax),'single');    %max kinetic energy
    fwrite(fid, importantStuffs(imin),'single');    %min kinetic energy
    fwrite(fid, m.Data.NINC_PHSP_SHORT,'single');

    %header finished, now write those weird 3 bytes that are formatted in:
    fwrite(fid, [0 0 0],'uint8');

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
    fclose(fid);
    %end
end
toc;

scatter(phspSrc(:,6)*10,phspSrc(:,7)*10)
hold on
n = 155;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 156;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 157;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 136;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 137;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 138;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 174;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 175;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])
n = 176;
scatter( [(stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2 (stf(1).ray(n).rayPos_bev(1)+2.5)/2 (stf(1).ray(n).rayPos_bev(1)-2.5)/2],[ (stf(1).ray(n).rayPos_bev(3) + 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) - 2.5)/2 (stf(1).ray(n).rayPos_bev(3) + 2.5)/2])


%relevant beamlets are numbers 155, 156, 157, 136, 137, 138, 174, 175, 176

hold off

% %loop over beamlets in -one- beam (make whole loop later)
% for i = 1:stf(1).numOfRays
%     %get stuff from phspSrc that's in the beamlet
%     I = find( ( (phspSrc(:,6)*10 - stf(1).ray(i).rayPos_bev(1)/2) >= -2.5/2 & (phspSrc(:,6)*10 - stf(1).ray(i).rayPos_bev(1)/2) < 2.5/2 ) & ( (phspSrc(:,7)*10 - stf(1).ray(i).rayPos_bev(3)/2) >= -2.5/2 & (phspSrc(:,7)*10 - stf(1).ray(i).rayPos_bev(3)/2) < 2.5/2 ) )
%     
%     
% end

%mode = "MODE0";
% write the header:
%fwrite(fid,mode,'char*1');

% [communist geologists]
% Ivan: look at all this mica
% Igor: no comrade, it's ourca.


% the variables in the header:

% mode (MODE0 or MODE2)
% NUM_PHSP_TOT;  total number of particles
% PHOT_PHSP_TOT;  no. of photons
% EKMAX_PHSP_SHORT; max k.e. of particles
% EKMINE_PHSP_SHORT; min. k.e. of electrons
% NINC_PHSP_SHORT;  no. of particles incident from original source
% this last one should be equal to what's in the source header file
%so like with the example as of now, it's 10000
