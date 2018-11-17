function matRad_exportCtEgs(ct,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ASCII CT export for EGSnrc
% 
% call
%   matRad_exportCtEgs(ct,filename)
%
% input
%   ct:             matRad ct struct
%   filename:       path where CTfile is created and its name
%
%
% References
% 
% A note on how this whole thing works:
% it prints out CT info like materials present and their numbers, then
% prints out two 3D CT cubes, going one slice (each slice a Z value) at a
% time.
% The first cube is the material number for each voxel
% The second cube is the density of each voxel
% BUT it's currently set to print the relative electron densities in the
% density cube. If you want the physical density, follow my instructions on
% what to comment/uncomment.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output formats
formatDims = " %f";
formatDims = " %f";
formatDims = " %f";
formatdim = "%d   %d   %d\n";
formnum = "%d";
formx = "  %d";

fid = fopen(filename,'w');

%number of materials:
numMats = 4;
fprintf(fid, "%d\n", numMats);

% OKAY a word about how material information is ordered. This is important
% for this to work, so if, say you added some materials, you need to know
% how to input this extra information:
% Each material in following matrices has is naturally in the same place 
% from one matrix to the next. This place is dictated by the material's
% (upper bound of its) CT number: lowest to highest.
% Hence lowest (air) is first, higher (e.g. bone) is last.
% And when adding new materials, be sure to update the info around line 90
% REMEMBER THIS LIKE YOU REMEMBER THE FIFTH OF NOVEMBER


%the max possible ct numbers for each material (and other characteristics):
matNames = ["AIR700ICRU", "LUNG700ICRU", "ICRUTISSUE700ICRU", "ICRPBONE700ICRU"];    %material names
ctMax = [-950,-700,125,3000];   %material max ct values
matNum = [1,2,3,4]; %material number (used in those semi-visual blocks)
matDens = [1.2048E-03, 2.6000E-01, 1.0, 1.85];  %material densities; from 700 ICRU pegs4 data file

% just stick new materials into these arrays in the proper place (refer to
% ctMax to determine where it fits in with the other materials)

%write the material names in the file
fprintf(fid, "%s\n", matNames);

%ESTEPE values except they're not needed, so just dummy inputs
dummy = zeros([1,numMats]);
fprintf(fid, "%d  ", dummy);
fprintf(fid,"\n");

% write ct dimensions
fprintf(fid,formatdim,ct.cubeDim);

% write voxel corner location in cm in physical cs with ct cube corner at [.5 .5 .5]
X = (.5:(ct.cubeDim(1)+.5))*ct.resolution.x/10;
Y = (.5:(ct.cubeDim(2)+.5))*ct.resolution.y/10;
Z = (.5:(ct.cubeDim(3)+.5))*ct.resolution.z/10;

fprintf(fid, formatDims, X);
fprintf(fid,"\n");
fprintf(fid, formatDims, Y);
fprintf(fid,"\n");
fprintf(fid, formatDims, Z);
fprintf(fid,"\n");

% the following bit is for doing physical density instead of relative
% electron density (so converting the ct cube to physical density)
% ======== {

% --------
%defining red2md = matrix with information needed to convert Relative
% Electron Density (red) to Mass Density (md):

%first we define properties of the elements that comprise the materials:
%each spot in the following arrays (Zi, Ai, the Wis) corresponds to an
%element found in at least one of the materials
Zi = [1,6,7,8,11,12,15,16,17,18,19,20,30];  %atomic number
Ai = [1.008 12.011 14.007 15.999 22.990 24.312 30.974 32.064 35.453 39.948 39.102 40.080 65.37];    %atomic mass

%weights of each element in the given materials:
wiAir = [0 1.24000E-04 7.55200E-01 2.31800E-01 0 0 0 0 0 1.28300E-02 0 0 0];
wiLung = [1.03000E-01 1.05000E-01 3.10000E-02 7.49000E-01 2.00000E-03 0 2.00000E-03 3.00000E-03 3.00000E-03 0 2.00000E-03 0 0];
wiTissue = [1.01172E-01 1.11000E-01 2.60000E-02 7.61828E-01 0 0 0 0 0 0 0 0 0];
wiBone = [4.72340E-02 1.44330E-01 4.19900E-02 4.46096E-01 0 2.20000E-03 1.04970E-01 3.15000E-03 0 0 0 2.09930E-01 1.00000E-04];

red2md = zeros(numMats,length(wiAir));

red2md(1,:) = wiAir;
red2md(2,:) = wiLung;
red2md(3,:) = wiTissue;
red2md(4,:) = wiBone;

% --------

edens_water = 3.343E23; %electron density of water, in 1/cm^3
N_A = 6.0221E23;    %avogadro's number

%calculate material relative electron densities, store in an array:
eDens = zeros(numMats,2);
eDens(:,2) = matDens;   %second column is physical density of the materials we're using
for i = 1:numMats
    summation = 0;
    for j = 1:length(wiAir)
       summation = summation + (Zi(j)/Ai(j))*red2md(i,j);
    end
    %fill the first column with the electron density of those materials
    eDens(i,1) = (matDens(i)*N_A*summation)/edens_water;
end
%make sure it's sorted from lowest to highest relative electron density:
eDens = sortrows(eDens, 1);

% ========= }

%write voxel material numbers (e.g. 1 for air, 4 for bone...)
for z = 1:ct.cubeDim(3)
    for x = 1:ct.cubeDim(1)
        for y = 1:ct.cubeDim(2)
            
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
        end
        
        fprintf(fid, "\n");
    end
    fprintf(fid, "\n");
end

fprintf("first loop done; material number cube has been written.\n");

%write voxel densities:
for z = 1:ct.cubeDim(3)
    for x = 1:ct.cubeDim(1)
        for y = 1:ct.cubeDim(2)
            
            % UNCOMMENT THE FOLLOWING BLOCK TO USE PHYSICAL DENSITY CUBE
            % INSTEAD OF RELATIVE ELECTRON DENSITY CUBE
            % (be sure to also comment out the last line in this loop)
            
%             %interpolate to get density of voxel:
%             if(ct.cube{1}(x,y,z) >= 0)
%                 if(ct.cube{1}(x,y,z) > eDens(numMats,1))
%                     %if electron density of voxel is greater than highest
%                     %electron density of the materials we're using
%                     slope = (eDens(numMats,2) - eDens(numMats - 1, 2))/(eDens(numMats,1) - eDens(numMats-1,1));
%                     intercept = eDens(numMats,2) - slope*eDens(numMats,1);
%                     physDens = slope*ct.cube{1}(x,y,z) + intercept;
%                 elseif(ct.cube{1}(x,y,z) < eDens(1,1))
%                     %if electron density of the voxel is less than the
%                     %lowest electron density of the materials we're using
%                     slope = (eDens(2,2) - eDens(1, 2))/(eDens(2,1) - eDens(1,1));
%                     intercept = eDens(2,2) - slope*eDens(2,1);
%                     physDens = slope*ct.cube{1}(x,y,z) + intercept;
%                 else
%                     %if it's between the upper and lowe bounds, then
%                     interpolate between the points with electron
%                     densities directly adjacent this voxel's electron
%                     density
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
% 
%             %finally, print the voxel density
%             fprintf(fid,formx,physDens);
            
            fprintf(fid,formx,ct.cube{1}(x,y,z));
        end  
        fprintf(fid, "\n");
    end
    fprintf(fid, "\n");
end

fprintf("wrote the file :)\n");

fclose(fid);