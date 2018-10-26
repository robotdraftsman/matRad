function EgsOptions = matRad_egsOptions(pln,ct)

%% run options

% number of paralle MC simulations
if isfield(pln.propDoseCalc.EgsOptions,'numOfParMCSim')
    EgsOptions.run.numOfParMCSim    = pln.propDoseCalc.EgsOptions.numOfParMCSim;
else
    EgsOptions.run.numOfParMCSim    = 4;
end
if isunix && EgsOptions.run.numOfParMCSim > 1
    EgsOptions.run.numOfParMCSim = 1;
end

% number of histories per bixel
if isfield(pln.propDoseCalc.EgsOptions,'nCasePerBixel')
    EgsOptions.run.nCasePerBixel    = pln.propDoseCalc.EgsOptions.nCasePerBixel;
else
    EgsOptions.run.nCasePerBixel    = 5000;
end

% relative dose cutoff
EgsOptions.run.relDoseCutoff    = 10^(-3);

% version (Carleton, dkfz, etc.)
EgsOptions.run.version = pln.propDoseCalc.EgsOptions.version;

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
switch pln.propDoseCalc.EgsOptions.version
    case 'Carleton'
%         switch pln.propDoseCalc.EgsOptions.source
%             case 'phsp'
%                 load('CALIBRATION_PHANTOM_TOH_VMC.mat','d_50mm','d_50mm_error')
%                 EgsOptions.run.absCalibrationFactorVmc      = 1./d_50mm;
%                 EgsOptions.run.absCalibrationFactorVmc_err  = d_50mm_error./(d_50mm.^2);
%         end
    case 'dkfz'
        EgsOptions.run.absCalibrationFactorVmc  = 99.818252282632300;
end

%% source

EgsOptions.source.myName       = 'some_source';                                         % name of source
EgsOptions.source.monitorUnits = 1;
switch pln.propDoseCalc.EgsOptions.source
    case 'beamlet'
    EgsOptions.source.spectrum     = fullfile(runsPath,'spectra','var_6MV.spectrum');   % energy spectrum source (only used if no mono-Energy given)
    EgsOptions.source.charge       = 0;                                                 % charge (-1,0,1)
    EgsOptions.source.type         = 'beamlet';
    
    case 'phsp'
    EgsOptions.source.particleType  = 2;
    EgsOptions.source.type          = 'phsp';
end

%% transport parameters

EgsOptions.McParameter.automatic_parameter  = 'yes';                       % if yes, automatic transport parameters are used
EgsOptions.McParameter.spin                 = 0;                           % 0: spin effects ignored; 1: simplistic; 2: full treatment

%% MC control

EgsOptions.McControl.ncase  = EgsOptions.run.nCasePerBixel;                % number of histories
EgsOptions.McControl.nbatch = 10;                                          % number of batches

%% variance reduction

EgsOptions.varianceReduction.repeatHistory      = 0.041;
EgsOptions.varianceReduction.splitPhotons       = 1;   
EgsOptions.varianceReduction.photonSplitFactor  = -80;

%% quasi random numbers

EgsOptions.quasi.base      = 2;                                                 
EgsOptions.quasi.dimension = 60;                                             
EgsOptions.quasi.skip      = 1;

%% geometry
switch pln.propDoseCalc.EgsOptions.version
    case 'Carleton'
        EgsOptions.geometry.XyzGeometry.methodOfInput = 'MMC-PHANTOM';  % input method ('CT-PHANTOM', 'individual', 'groups')
    case 'dkfz'
        EgsOptions.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';   % input method ('CT-PHANTOM', 'individual', 'groups')
end
EgsOptions.geometry.dimensions          = ct.cubeDim;
EgsOptions.geometry.XyzGeometry.Ct      = 'CT';                         % name of geometry

%% scoring manager
EgsOptions.scoringOptions.startInGeometry               = 'CT';            % geometry in which partciles start their transport
EgsOptions.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
EgsOptions.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
EgsOptions.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
EgsOptions.scoringOptions.outputOptions.dumpDose        = pln.propDoseCalc.EgsOptions.dumpDose;               % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)



end