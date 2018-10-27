function dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
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

%if need to check/move files around on/to/from cluster, put the password:
%password = 


% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 1;

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'vmc++'))
end

checkForFiles = 0;

egsinpbase = 'inputs';
egsinpPath = '/Users/sakinahussain/Documents/GitHub/matRad/EGSnrc/egsinpFiles';
dosbase = 'inputs';
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
   numOfBixelsContainer = 1;%ceil(dij.totalNumOfBixels/10);
end





% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(1,1);
doseTmpContainerError = cell(1,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
    dij.physicalDoseError{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end

%keeping these - check that the absolute calibration factor is good tho

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
% absCalibrationFactorVmc = 99.818252282632300;

%the below is: TOHCC calibration (10x10 cm2 field, depth=5cm, SAD = 100cm) 6.495E-17
absCalibrationFactorEgs = 1/9.4457e-017;
absCalibrationFactorEgs_err = 1.022147627023901E-002.*absCalibrationFactorEgs;

phantomPath = 'EGSnrc/phantoms';
phantomName = 'matRad_CT.egsphant';

% export CT cube as ASCII file for EGSnrc
%also get the x,y,z coordinates of the centre of the ct cube (don't think I
%need to keep this part though...can't see any later application for it?)

%##############################################
%ctCubeCentre = matRad_exportCtEgs(ct, fullfile(phantomPath, [phantomName,'.egsphant']));

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


%check if we have the files:
%first check for .dos files. If they're present, then great, go straight to
%dose read-in. Don't even need to generate an egsphant file
%If not, check for egsphant, egsinp and phsp files on the cluster. If they're there,
%tell user to run dosxyz on them.
%If the egsphant file is not there, check here. If here, transfer it to the
%cluster. If not here, then generate the CT egsphant and transfer it.
%If the egsip files are not there, check here. If they're here, transfer
%them to the cluster. If they're not here, then generate them and transfer
%them to the cluster.
%If the phsp files aren't there, check here and if here, transfer them. If
%they're not here either, then be very sad and tell the user to run
%readANDwriteBinaryPHSP on their desired phsp source.


xmin = -45; %these are in mm
xmax = 45;
ymin = -45;
ymax = 45;
n = 1;
allBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);

for x = xmin:5:xmax
    for y = ymin:5:ymax
        allBeamlets(n,:) = [x y];
        n = n + 1;
    end
end

whichBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1)  ,length(stf));  %beamlets included (index in a given row) for each beam (columns)
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        indexOfPhsp = find( (allBeamlets(:,1) == stf(i).ray(j).rayPos_bev(1)) & (allBeamlets(:,2) == -stf(i).ray(j).rayPos_bev(3)) );    %find index in reference vector that corresponds to this beamlet
        whichBeamlets(j,i) = indexOfPhsp;
    end
end



    clusterPath = '/data/data060/shussain/egsnrc/dosxyznrc/';
    for i = 1:dij.numOfBeams % loop over all beams
        
        for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
            %check if dos exists here. If not, check cluster...
            dosfile = strcat(dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos');
            if (exist(fullfile(dosPath,dosfile)) == 0)
                fprintf("%s ain't here so I'm gonna get it\n", dosfile);
                %try to get the dos file from the cluster
                try
                    scp_simple_get('tyr.physics.carleton.ca','shussain',password,dosfile,'/Users/sakinahussain/Documents/GitHub/matRad/EGSnrc/EGSdosFiles','/data/data060/shussain/egsnrc/dosxyznrc/');
                catch
                    fprintf("couldn't get it from the cluster...\n");
                    %if the try gives an error, the file's not there and we
                    %have to go about making it. So we check if the ingredients
                    %to make it are there (the egsinp file). If it's not here,
                    %make the file and transfer it.

                    egsinpFile = strcat(egsinpbase,'Beam',num2str(i),'Beamlet',num2str(j),'.egsinp');
                    fromCluster = ssh2_simple_command('tyr.physics.carleton.ca','shussain',password,['ls ',clusterPath,egsinpFile]);
                    if (length(fromCluster{1}) == 0)    %if it's empty -> no file on cluster
                         if (exist(fullfile(egsinpPath,egsinpFile)) == 0) %if file doesn't exist locally, then make it
                            %make the egsinp file for this beamlet:

                            matRad_createEgsinp(stf,pln,phantomName,egsinpFile,whichBeamlets(j,i),i,j);

                         end
                        %scp it to the cluster:
                         try
                             fprintf("putting egsip file for Beam %d beamlet %d onto cluster...\n",i,j);
                             scp_simple_put('tyr.physics.carleton.ca','shussain',password,egsinpFile,'/data/data060/shussain/egsnrc/dosxyznrc/','/Users/sakinahussain/Documents/GitHub/matRad/EGSnrc/egsinpFiles',thisegsinpfile); 
                         catch
                             fprintf("Didn't transfer the file %s to the cluster: encountered a problem :| \n",egsinpFile);
                         end
                    end
                    fprintf("Please run dosxyznrc on %s on the cluster, convert it to a .dos file, and bring it back to me (in /matRad/EGSdosFiles/).\n",egsinpFile);
                end
            end
        end
    end

%%%%% reading in the .dos file from dosxyz into dij matrix:

for i = 1:dij.numOfBeams % loop over all beams
   tic;
   % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end
        
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        %readCounter = readCounter + 1;
        
        %can do the checking in here instead? But rn can at least test the
        %checker up there...
        
        writeCounter = writeCounter + 1;
        
        %flag: donotneed?
        % create different seeds for every bixel
        %EgsOptions.McControl.rngSeeds = [randi(30000),randi(30000)];

        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(writeCounter)  = i;
           dij.rayNum(writeCounter)   = j;
           dij.bixelNum(writeCounter) = j;
        end
        beamNum(writeCounter)  = i;
        rayNum(writeCounter)   = j;
        bixelNum(writeCounter) = j;
        
        % set ray specific vmc++ parameters
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        rayCorner1 = (stf(i).ray(j).rayCorners_SCD(1,:) + stf(i).isoCenter)/10;              
        rayCorner2 = (stf(i).ray(j).rayCorners_SCD(2,:) + stf(i).isoCenter)/10;
        rayCorner3 = (stf(i).ray(j).rayCorners_SCD(3,:) + stf(i).isoCenter)/10; %vmc needs only three corners (counter-clockwise)
                
        beamSource = (stf(i).sourcePoint + stf(i).isoCenter)/10;
        
        % b) swap x and y (CT-standard = [y,x,z])
        rayCorner1 = rayCorner1([2,1,3]);              
        rayCorner2 = rayCorner2([2,1,3]);
        rayCorner3 = rayCorner3([2,1,3]);
        beamSource  = beamSource([2,1,3]);
        
        % c) set vmc++ parameters
        % EgsOptions.beamletSource.monoEnergy                = stf(i).ray(j).energy;                 % photon energy
        % %EgsOptions.beamletSource.monoEnergy                 = []                  ;                  % use photon spectrum
        % EgsOptions.beamletSource.beamletEdges               = [rayCorner1,rayCorner2,rayCorner3];    % counter-clockwise beamlet edges
        % EgsOptions.beamletSource.virtualPointSourcePosition = beamSource;                            % virtual beam source position
                                
            %%%%%%%%%%%%%%%%%%%% this is important stuff from the VMC things:
            % import calculated dose
                readCounter = readCounter + 1;
                
                % Display progress
                if verbose == 0
                   % matRad_progress(readCounter,dij.totalNumOfBixels);
                end
                
                % update waitbar
                waitbar(writeCounter/dij.totalNumOfBixels);
                
                [bixelDose,bixelDoseError] = matRad_readDoseVmc(fullfile(dosPath, [dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos']),EgsOptions);
                
                
                if ~calcDoseDirect
                    % if not calculating dose directly, sample dose
                    
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
                
                %{
                % apply absolute calibration factor
                %bixelDose = bixelDose*absCalibrationFactorVmc;
                bixelDose = bixelDose*absCalibrationFactorEgs;

                % Save dose for every bixel in cell array
                doseTmpContainer{1,1} = sparse(V,1,bixelDose(V),dij.numOfVoxels,1);
                
                % save computation time and memory by sequentially filling the 
                % sparse matrix dose.dij from the cell array
                    if calcDoseDirect
                        if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)),'weight')
                            % score physical dose
                            dij.physicalDose{1}(:,i) = dij.physicalDose{1}(:,i) + stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight * doseTmpContainer{1,1};
                        else
                            error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                        end
                        fprintf("lookie therey...\n");
                    else
                        % fill entire dose influence matrix
                        
                        dij.physicalDose{1}(:,readCounter) = ...
                            [doseTmpContainer{1,1}];
                    end
                %%%%%%%%%%%%%%%%%%%%%%
                %fprintf("ok that was beam %d beamlet %d\n",i,j);
                %}
            
    end
    mins = toc;
    fprintf('beam %d done',i);
    try
        text = strcat('got beam ' , num2str(i) , ' finished and it took ' , num2str(mins) , ' seconds');
        example_txtmsg('dos file read-in:',text);
    catch
        fprintf('%f seconds have passed',mins);
    end
end




% % delete temporary files
% delete(fullfile(phantomPath, 'matRad_CT.egsphant'));     % egs phantom file
% 
% for i = 1:length(stf)
%     for j = 1:stf(i).numOfRays
%         % deletes egsinp files:
%         delete(fullfile(egsinpPath, [egsinpbase,'Beam',num2str(i),'Beamlet',num2str(j),'.egsinp'])); % vmc inputfile
%         
%         % deletes the .dos files we made from 3ddose files:
%         delete(fullfile(dosPath, [dosbase,'Beam',num2str(i),'Beamlet',num2str(j),'.dos']));    % dosxyz output file in vmc output file format
%     end
% end
%     


try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end