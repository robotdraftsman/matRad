function dij = matRad_calcPhotonDoseEgs(ct,stf,pln,cst,planName,nCasePerBixel,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,numOfParallelMCSimulations)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   planName:                   matRad plan name
%   nCasePerBixel:              number of photons simulated per bixel
%   calcDoseDirect:             boolian switch to bypass dose influence matrix
%                               computation and directly calculate dose; only makes
%                               sense in combination with matRad_calcDoseDirect.m%
% output
%   dij:                        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

%if need to check/move files around on/to/from cluster, enter your login information:
clusterUserName = '';
password = '';
address = 'tyr.physics.carleton.ca';
clusterPath = '/data/data060/shussain/egsnrc/dosxyznrc/';

% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 1;

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'vmc++'))
end

%egsinpbase = 'inputs';
egsinpbase = planName;
egsinpPath = 'EGSnrc/egsinpFiles';
%dosbase = 'inputs';
dosbase = planName;
dosPath = 'EGSnrc/EGSdosFiles';

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
   numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
   numOfBixelsContainer = 1;
end





% set up arrays for book keeping
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(1,1);
doseTmpContainerError = cell(1,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
    dij.physicalDoseError{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end


% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel

%the below is: TOHCC calibration (10x10 cm2 field, depth=5cm, SAD = 100cm) 6.495E-17
absCalibrationFactorEgs = 1/9.4457e-017;
absCalibrationFactorEgs_err = 1.022147627023901E-002.*absCalibrationFactorEgs;

% get default vmc options
EgsOptions = matRad_vmcOptions(pln,ct);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

writeCounter                  = 0;
readCounter                   = 0;


% initialize waitbar
figureWait = waitbar(0,'EGSnrc photon dose influence matrix calculation.. ');

fprintf("matRad: EGSnrc photon dose calculation... ");

% use this to match input files to phase space files
xmin = -45; %these are in mm. These are from the stf file
xmax = 45;
ymin = -45;
ymax = 45;
n = 1;
allBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);

% allBeamlets tracks which bixels are included in each beam
for x = xmin:5:xmax
    for y = ymin:5:ymax
        allBeamlets(n,:) = [x y];
        n = n + 1;
    end
end

% whichBeamlets which beamlet phsp file belongs to each bixel
whichBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1)  ,length(stf));  %beamlets included (index in a given row) for each beam (columns)
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        indexOfPhsp = find( (allBeamlets(:,1) == stf(i).ray(j).rayPos_bev(1)) & (allBeamlets(:,2) == -stf(i).ray(j).rayPos_bev(3)) );    %find index in reference vector that corresponds to this beamlet
        whichBeamlets(j,i) = indexOfPhsp;
    end
end


phantomPath = 'EGSnrc/phantoms';
phantomName = strcat(planName,'_matRad_CT','.egsphant');

runOnCluster = 0;   %if don't need to do anything on cluster (i.e. all the
%dos or at least 3ddose files exist) then it's 0.
%But if you need to run DOSXYZnrc on even one thing it's flipped to 1

for i = 1:dij.numOfBeams % loop over all beams

    for j = 1:stf(i).numOfRays % loop over all rays
        caught3ddose = 0;
        caughtEgsinp = 0;
        %check if dos exists here. If not, check cluster...
        dosfile = strcat(dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos');
        if (exist(fullfile(dosPath,dosfile),'file') == 0)
            fprintf("%s ain't here so I'm gonna get it\n", dosfile);
            %try to get the dos file from the cluster
            try
                scp_simple_get(address,clusterUserName,password,dosfile,'EGSnrc/EGSdosFiles',clusterPath);
            catch
                fprintf("couldn't get it from the cluster...\n");
                %if the try gives an error, the file's not there and we
                %have to go about making it. So we check if the ingredients
                %to make it are there (the 3ddose or at least the egsinp file).
                %If 3ddose file is there, tell user to convert it to .dos
                %format. If it's not there check here for the egsinp file.
                %If it's not here, make the egsinp file and transfer it.
                onCluster3ddt{1} = [];
                file3ddt = strcat(egsinpbase,'Beam',num2str(i),'Beamlet',num2str(j),'.3ddose');
                try
                    onCluster3ddt = ssh2_simple_command(address,clusterUserName,password,['ls ',clusterPath,file3ddt]);
                catch
                    fprintf("Could not access cluster to check for 3ddose file. Proceeding as if required resources are not on cluster...\n");
                    caught3ddose = 1;
                end
                if(isempty(onCluster3ddt{1}) || caught3ddose == 1) % if 3ddose file is not on cluster
                    runOnCluster = 1;
                    %then we wanna do the egsinp checking/creating:
                    fromCluster{1} = [];
                    egsinpFile = strcat(egsinpbase,'Beam',num2str(i),'Beamlet',num2str(j),'.egsinp');
                    try
                        fromCluster = ssh2_simple_command(address,clusterUserName,password,['ls ',clusterPath,egsinpFile]);
                    catch
                        fprintf("Could not access cluster to check for egsinp file. Proceeding as if required resources are not on cluster...\n");
                    end
                    
                    if (isempty(fromCluster{1}) || caughtEgsinp == 0)    %if it's empty -> no file on cluster or can't access cluster
                         if (exist(fullfile(egsinpPath,egsinpFile),'file') == 0) %if file doesn't exist locally, then make it
                            %make the egsinp file for this beamlet:
                            matRad_createEgsinp(stf,pln,phantomName,egsinpFile,whichBeamlets(j,i),i,j);
                         end
                        %scp it to the cluster:
                         try
                             fprintf("putting egsip file for Beam %d beamlet %d onto cluster...\n",i,j);
                             scp_simple_put(address,clusterUserName,password,egsinpFile,clusterPath,'EGSnrc/egsinpFiles',egsinpFile); 
                         catch
                             fprintf("Didn't transfer the file %s to the cluster: encountered a problem :| \n",egsinpFile);
                             fprintf("Please move the file to the cluster manually.\n");
                         end
                    end
                    fprintf("Please run dosxyznrc on %s on the cluster, convert resulting .3ddose file to a .dos file, and bring it back to me (in /matRad/EGSdosFiles/).\n",egsinpFile);
                    
                else %if it IS on the cluster:
                    fprintf("please convert %s to .dos format using 3ddose_to_dos.c or automate_3ddose_to_dos.c.\n",file3ddt);
                end                
            end
        end
    end
end


%this makes the CT phantom for DOSXYZnrc, if it's necessary to run
%DOSXYZnrc on the cluster and if the phantom doesn't already exist:
%NOTE: this creates a phantom with relative electron densities, not
%physical densities. The code to do physical densities can be uncommented
%in matRad_exportCtEgs.m to convert RED to physical density

caughtEgsphant = 0;
%if the user has to run DOSXYZnrc on the cluster, they're going to need the phantom:
if(runOnCluster == 1)
    %check on cluster for phantom.
    phantomOnCluster{1} = [];
    try
        phantomOnCluster = ssh2_simple_command(address,clusterUserName,password,['ls ',clusterPath,phantomName]);
    catch
        fprintf("Could not access cluster to check for egsphant file. Proceeding as if required resources are not on cluster...\n");
        caughtEgsphant = 1;
    end
    if (isempty(phantomOnCluster{1}) || caughtEgsphant == 0)    %if it's empty -> no file on cluster or can't access cluster
        %If not there check here and transfer to cluster if it is here. 
        if (exist(fullfile(phantomPath,phantomName),'file') == 0) %if file doesn't exist locally, then make it:
            % export CT cube as ASCII file for EGSnrc
            matRad_exportCtEgs(ct, fullfile(phantomPath, phantomName));
        end
        %scp it to the cluster:
        try
            fprintf("Putting the phantom onto cluster...\n");
            scp_simple_put(address,clusterUserName,password,egsinpFile,clusterPath,'EGSnrc/egsinpFiles',phantomName); 
        catch
            fprintf("Didn't transfer the file %s to the cluster: encountered a problem :| \n",phantomName);
            fprintf("Please move the file to the cluster manually.\n");
        end
    end
    fprintf("Please run dosxyznrc on the cluster using this phantom and the given .egsinp files.\n");
end



%%%%% reading in the .dos file from dosxyz into dij matrix:
%AKA the REAL MEAT of this code, everything else is mere prep

for i = 1:dij.numOfBeams % loop over all beams
   % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
    end
        
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        writeCounter = writeCounter + 1;

        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(writeCounter)  = i;
           dij.rayNum(writeCounter)   = j;
        end
        beamNum(writeCounter)  = i;
        rayNum(writeCounter)   = j;
        
        %%%%%%%%%%%%%%%%%%%% this is important info from the VMC things:
        % import calculated dose
            readCounter = readCounter + 1;

            % Display progress
            if verbose == 0
                matRad_progress(readCounter,dij.totalNumOfBixels);
            end

            % update waitbar
            waitbar(writeCounter/dij.totalNumOfBixels);
            
            %read in the dose and its error from the .dos file
            [bixelDose,bixelDoseError] = matRad_readDoseVmc(fullfile(dosPath, [dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos']),EgsOptions);

            if ~calcDoseDirect
                % if not calculating dose directly, sample dose
                % so we implement importance sampling and compute associated errors

                % determine cutoff
                doseCutoff          = relDoseCutoff*max(bixelDose);

                % determine which voxels to sample
                indSample = bixelDose < doseCutoff & bixelDose ~= 0;
                r = rand(nnz(indSample),1);

                % sample them
                thresRand = bixelDose(indSample)./doseCutoff;
                indKeepSampled = r < thresRand;
                indKeep = indSample;
                indKeep(indKeep) = indKeepSampled;

                clear indKeepSampled;
                bixelDose(indKeep) = doseCutoff;
                bixelDose(indSample & ~indKeep) = 0;
                clear indKeep;
                clear indSample;

                % coefficients and terms for the error
                erfArg = (doseCutoff-bixelDose)./(bixelDoseError.*sqrt(2));
                %can uncomment the following, then comment out the
                %bixelDoseError computation. Just rewrote it so wouldn't
                %have all these big matrices hanging around
                
%                     erfp_coeff  = (bixelDose.*doseCutoff - bixelDose.^2);
%                     erfp_term   = (1+erf(erfArg))./2;

%                     erfm_coeff  = bixelDoseError.^2;
%                     erfm_term   = (1-erf(erfArg))./2;

%                     gauss_coeff = bixelDose.*bixelDoseError.^2;
%                     gauss_term  = normpdf(doseCutoff,bixelDose,bixelDoseError);
%                     
                % clear NaNs from bixelDoseError = 0
%                     gauss_term(isnan(gauss_term)) = 0;

                %bixelDoseError = sqrt( erfp_coeff.*erfp_term + erfm_coeff.*erfm_term + gauss_coeff.*gauss_term );
                bixelDoseError = sqrt( (bixelDose.*doseCutoff - bixelDose.^2).*((1+erf(erfArg))./2) + (bixelDoseError.^2).*((1-erf(erfArg))./2) + (bixelDose.*bixelDoseError.^2).*(normpdf(doseCutoff,bixelDose,bixelDoseError)) );

            end

            % apply absolute calibration factor
            bixelDoseError  = sqrt((absCalibrationFactorEgs.*bixelDoseError).^2+(bixelDose.*absCalibrationFactorEgs_err).^2);
            bixelDose       = bixelDose*absCalibrationFactorEgs;

            % Save dose for every bixel in cell array
            doseTmpContainer{1,1}       = sparse(V,1,bixelDose(V),dij.numOfVoxels,1);
            doseTmpContainerError{1,1}  = sparse(V,1,bixelDoseError(V),dij.numOfVoxels,1);

            % save computation time and memory by sequentially filling the
            % sparse matrix dose.dij from the cell array
            if calcDoseDirect
                if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)),'weight')
                    % score physical dose
                    dij.physicalDose{1}(:,i)        = dij.physicalDose{1}(:,i) + stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight{1} * doseTmpContainer{1,1};
                    dij.physicalDoseError{1}(:,i)   = sqrt(dij.physicalDoseError{1}(:,i).^2 + (stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight{1} * doseTmpContainerError{1,1}).^2);
                else
                    error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                end
            else
                % fill entire dose influence matrix
                dij.physicalDose{1}(:,readCounter) = ...
                    [doseTmpContainer{1:mod(readCounter-1,numOfBixelsContainer)+1,1}];

                dij.physicalDoseError{1}(:,readCounter) = ...
                    [doseTmpContainerError{1:mod(readCounter-1,numOfBixelsContainer)+1,1}];
            end

    end
    fprintf('beam %d done',i);
    try
        text = strcat('got beam ' , num2str(i) , ' finished and it took ' , num2str(mins) , ' seconds');
        example_txtmsg('dos file read-in:',text);
    catch
        fprintf('%f seconds have passed',mins);
    end
end


% % delete temporary files if you're DARING enough and know that you won't 
% have to rerun anything :o
% delete(fullfile(phantomPath, phantomName));     % egs phantom file
% 
% for i = 1:length(stf)
%     for j = 1:stf(i).numOfRays
%         % deletes egsinp files:
%         delete(fullfile(egsinpPath, [egsinpbase,'Beam',num2str(i),'Beamlet',num2str(j),'.egsinp'])); % egs inputfile
%         
%         % deletes the .dos files we made from 3ddose files:
%         delete(fullfile(dosPath, [dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos']));    % dosxyz output file in vmc output file format
%     end
% end


try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
  delete(allWaitBarFigures);
  pause(0.1);
catch
end