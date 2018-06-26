function middle = matRad_exportCtEgs(ct,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ASCII CT export for EGSnrc
% 
% call
%   middle = matRad_exportCtEgs(ct,filename)
% where "middle" is an array with the midpoint of the ct cube
% -> is used in calcPhotonDoseEgs
%
% input
%   ct:             matRad ct struct
%   filename:       path where CTfile is created
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output formats
formatx = " %f";
formaty = " %f";
formatz = " %f";
formatdim = "%d   %d   %d\n";
formnum = "%d";
formx = "  %d";

fid = fopen(filename,'w');

%number of materials:
numMats = 4;
fprintf(fid, "%d\n", numMats);

%the max possible ct numbers for a material (and other characteristics):
%order (important!): lowest-to-highest ct number
matNames = ["AIR700ICRU", "LUNG700ICRU", "ICRUTISSUE700ICRU", "ICRPBONE700ICRU"];    %material names
ctMax = [-950,-700,125,3000];   %material max ct values (in same order as matNames)
matNum = [1,2,3,4]; %material number (in those semi-visual blocks)
matDens = [1.2048E-03, 2.6000E-01, 1.0, 1.85];
%doing things like this should make it easier to add materials;
%just keep adding the stuff into these arrays
%material densities are from 700 ICRU pegs4 data file

    
%write the material names in the file
fprintf(fid, "%s\n", matNames);

%ESTEPE values except they're not needed ? So just dummy inputs
dummy = zeros([1,numMats]);
fprintf(fid, "%d  ", dummy);
fprintf(fid,"\n");

% write ct dimensions
fprintf(fid,formatdim,ct.cubeDim);

% write voxel corner location in cm in physical cs with ct cube corner at [.5 .5 .5]
X = [.5:(ct.cubeDim(1)+.5)]*ct.resolution.x/10;
Y = [.5:(ct.cubeDim(2)+.5)]*ct.resolution.y/10;
Z = [.5:(ct.cubeDim(3)+.5)]*ct.resolution.z/10;

%get the coordinates of the centre of the ct cube:
%assuming the voxel corner referred to above is top left back corner
middle = [ (X(length(X)) + X(1))/2, (Y(length(Y)) + Y(1))/2, (Z(length(Z)) + Z(1))/2 ];
%fprintf("(brakelaatabasaasta feed me) centre = (%f, %f, %f)\n", middle);


fprintf(fid, formatx, X);
fprintf(fid,"\n");
fprintf(fid, formaty, Y);
fprintf(fid,"\n");
fprintf(fid, formatz, Z);
fprintf(fid,"\n");

%setting up the 3d density matrix (filled in the nested for loops)
%ctDens = zeros(ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3));

%defining red2md = matrix with information needed to convert Relative
% Electron Density to Mass Density:

%first we define properties of the elements that comprise the materials:
Zi = [1,6,7,8,11,12,15,16,17,18,19,20,30];
Ai = [1.008 12.011 14.007 15.999 22.990 24.312 30.974 32.064 35.453 39.948 39.102 40.080 65.37];

wiAir = [0 1.24000E-04 7.55200E-01 2.31800E-01 0 0 0 0 0 1.28300E-02 0 0 0];
wiLung = [1.03000E-01 1.05000E-01 3.10000E-02 7.49000E-01 2.00000E-03 0 2.00000E-03 3.00000E-03 3.00000E-03 0 2.00000E-03 0 0];
wiTissue = [1.01172E-01 1.11000E-01 2.60000E-02 7.61828E-01 0 0 0 0 0 0 0 0 0];
wiBone = [4.72340E-02 1.44330E-01 4.19900E-02 4.46096E-01 0 2.20000E-03 1.04970E-01 3.15000E-03 0 0 0 2.09930E-01 1.00000E-04];

red2md = zeros(numMats,length(wiAir));

red2md(1,:) = wiAir;
red2md(2,:) = wiLung;
red2md(3,:) = wiTissue;
red2md(4,:) = wiBone;

edens_water = 3.343E23; %in 1/cm^3
N_A = 6.0221E23;


%calculate material relative electron densities, store in an array:
eDens = zeros(numMats,2);
eDens(:,2) = matDens;
for i = 1:numMats
    summation = 0;
    for j = 1:length(wiAir)
       summation = summation + (Zi(j)/Ai(j))*red2md(i,j);
    end
    eDens(i,1) = (matDens(i)*N_A*summation)/edens_water;
end
%make sure it's sorted from lowest to highest relative electron density:
eDens = sortrows(eDens, 1);


%write voxel material numbers (e.g. 1 for air, 4 for bone...)
for z = 1:ct.cubeDim(3);
    for x = 1:ct.cubeDim(1);
        for y = 1:ct.cubeDim(2);
            
            %see which material it is and set values accordingly:
            unassigned = 1;
            i = 1;
            while(unassigned == 1)
                  ctnum = ct.cubeHU{1}(x,y,z);    %ct number of the voxel
                      if(i>length(ctMax))
                          fprintf("this material has no place in this world :(");
                          unassigned = 0;
                      end
                      if(ctnum <= ctMax(i))
                          name = matNames(i);
                          number = matNum(i);
                          unassigned = 0;
                      else
                          if(i <= length(ctMax) - 1)
                              i=i+1;
                          else
                              name = matNames(i+1);
                              number = matNum(i+1);
                              unassigned = 0;
                          end
                      end

            end            
            %print the material number in that voxel
            fprintf(fid, formnum, number);
            
            %convert electron density to mass density according to equation
            % #1 here: https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4875976
%             sum = 0;
%             for j = 1:length(wiAir)
%                sum = sum + (Zi(j)/Ai(j))*red2md(i,j);
%             end
%             physDens = (edens_water*ct.cube{1}(x,y,z))/(N_A*sum);
%             ctDens(x,y,z) = physDens;
        end
        
        fprintf(fid, "\n");
    end
    fprintf(fid, "\n");
end

fprintf("first loop done\n");

%write voxel densities:
for z = 1:ct.cubeDim(3);
    for x = 1:ct.cubeDim(1);
        for y = 1:ct.cubeDim(2);            
            %interpolate to get density of voxel:
%             if(ct.cube{1}(x,y,z) >= 0)
%                 if(ct.cube{1}(x,y,z) > eDens(numMats,1))
%                     slope = (eDens(numMats,2) - eDens(numMats - 1, 2))/(eDens(numMats,1) - eDens(numMats-1,1));
%                     intercept = eDens(numMats,2) - slope*eDens(numMats,1);
%                     physDens = slope*ct.cube{1}(x,y,z) + intercept;
%                 elseif(ct.cube{1}(x,y,z) < eDens(1,1))
%                     slope = (eDens(2,2) - eDens(1, 2))/(eDens(2,1) - eDens(1,1));
%                     intercept = eDens(2,2) - slope*eDens(2,1);
%                     physDens = slope*ct.cube{1}(x,y,z) + intercept;
%                 else
%                     i = 2;
%                     if(ct.cube{1}(x,y,z) == eDens(1,1))
%                         physDens = eDens(1,2);
%                         i = numMats + 42;
%                     end
%                     while(i <= numMats)
%                         if(ct.cube{1}(x,y,z) == eDens(i,1))
%                             physDens = eDens(i,2);
%                             i = numMats + 17;
%                         else
%                             if(ct.cube{1}(x,y,z) < eDens(i,1))
%                                 slope = (eDens(i,2) - eDens(i-1,2))/(eDens(i,1) - eDens(i-1,1));
%                                 intercept = eDens(i,2) - slope*eDens(i,1);
%                                 physDens = slope*ct.cube{1}(x,y,z) + intercept;
%                                 i = numMats + 313;
%                             else
%                                 i = i + 1;
%                             end
%                         end
%                     end
%                 end
%             else
%                 fprintf("wth you have a negative electron density ?? at (x,y,z) = (%d, %d, %d)\n", x,y,z);
%             end
            %fprintf(fid,formx,physDens);
            fprintf(fid,formx,ct.cube{1}(x,y,z));
        end  
        fprintf(fid, "\n");
    end
    fprintf(fid, "\n");
end

fprintf("wrote the file :)\n");

fclose(fid);