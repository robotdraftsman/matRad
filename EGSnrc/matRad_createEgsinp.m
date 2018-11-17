function matRad_createEgsinp(stf, pln, phantomName, filebase, whichPhspBeamlet,n,m)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% like createVmcInput except for making .egsinp files
% these files are necessary to run DOSXYZnrc. One file per beamlet, so they
% have info on beamlet location, associated phase space file, the phantom
% the dose is being calculated in, and some run parameters.
% 
% stf:                  matRad stf struct
% pln:                  matRad pln struct
% phantomName:          name of the CT phantom .egsphant file used
% filebase:             the name of the files according to the convention:
%                       filebaseBeamIBeamletJ.egsinp. This filebase is used for
%                       all the files specific to a given run of matRad with
%                       EGSnrc.
% whichPhspBeamlet:     which phase space file (specified by a number) is
%                       connected to this beamlet (phsp files are numbered
%                       sequentially; these numbers don't correspond to
%                       beamlet numbers in the stf file)
% n:                    the beam number of the egsinp file being generated
% m:                    the beamlet number of the egsinp file being generated
% 
% NOTE: most settings copied straight from a DOSXYZ input file; only a few
% paramters are changed - those that aren't just direct text
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phspfilebase = "dividedPhsp36M";
egsinpPath = 'EGSnrc/egsinpFiles';

rec2 = 0;

thisegsinpfile = strcat(filebase, 'Beam',num2str(n),'Beamlet',num2str(m), '.egsinp');
sourcephspfile = strcat(phspfilebase,num2str(whichPhspBeamlet), '.egsphsp1');

%The angles in the matRad (DICOM) coordinate system:
thetaC = 0; %madRad collimator angle
thetaG = stf(n).gantryAngle;
thetaT = stf(n).couchAngle;

dsource = 50;   %distance in cm

%angles in dosxyznrc coordinate system:
%from: Lixin Zhan, Runqing Jiang and Ernest K Osei, "Beam coordinate transformations from DICOM to DOSXYZnrc"
dosxyzTheta = acosd(-sind(thetaT)*sind(thetaG));            
dosxyzPhi = atan2d( -cosd(thetaG),cosd(thetaT)*sind(thetaG) );
dosxyzPhiCol = thetaC - 90 + atan2d( (-sind(thetaT)*cosd(thetaG)),(cosd(thetaT)));

%isocentre in matRad CT phantom coordinate system-> don't need to transform
%divide by 10 to convert from mm to cm
isocentre = [pln.propStf.isoCenter(n,1) pln.propStf.isoCenter(n,2) pln.propStf.isoCenter(n,3) ]/10;

%valuess from 4th line/record in egsinp file:
%ECUTIN,PCUTIN: Electron (total) and photon global cutoff energies in MeV.
%SMAX: Dummy input
rec4 = [0.7 0.01 0];   

%zeroairdose,doseprint,MAX20: (binary 1/0):
rec5 = [0 0 0];

nhistories = 10000000;

%open the file and get writing:
file = fopen(fullfile(egsinpPath, thisegsinpfile),'w');

fprintf(file,"                                                                                 #!GUI1.0\n");
fprintf(file, "%d\n",rec2);
fprintf(file,"/data/data060/shussain/egsnrc/dosxyznrc/%s\n",phantomName);
fprintf(file,"%f, %f, %f\n", rec4);
fprintf(file,"%d, %d, %d\n", rec5);
fprintf(file,"2, 2, %f, %f, %f, %f, %f, %f, %f, 0, 0, 0, 0, 0\n",isocentre(1),isocentre(2),isocentre(3),dosxyzTheta,dosxyzPhi,dsource,dosxyzPhiCol);
fprintf(file,"2, 0, 0, 0, 0, 0, 0, 0\n");
fprintf(file,"/data/data060/shussain/egsnrc/dosxyznrc/%s\n",sourcephspfile);

%print record 13:
fprintf(file,"%d, 0, 99, 33, 97, 100.0, 0, 0, 0, 0, , -1, 0, 0, 1, 0, 0\n",nhistories);

mctparameter = strcat('#########################\n',...
    [':Start MC Transport Parameter:\n\n'],...
    ['Global ECUT= 0.7\n'],...
    ['Global PCUT= 0.01\n'],...
    ['Global SMAX= 5\n'],...
    ['ESTEPE= 0.25\n'],...
    ['XIMAX= 0.5\n'],...
    ['Boundary crossing algorithm= PRESTA-I\n'],...
    ['Skin depth for BCA= 0\n'],...
    ['Electron-step algorithm= PRESTA-II\n'],...
    ['Spin effects= On\n'],...
    ['Brems angular sampling= Simple\n'],...
    ['Brems cross sections= BH\n'],...
    ['Bound Compton scattering= Off\n'],...
    ['Compton cross sections= default\n'],...
    ['Pair angular sampling= Simple\n'],...
    ['Pair cross sections= BH\n'],...
    ['Photoelectron angular sampling= Off\n'],...
    ['Rayleigh scattering= Off\n'],...
    ['Atomic relaxations= Off\n'],...
    ['Electron impact ionization= Off\n'],...
    ['Photon cross sections= xcom\n'],...
    ['Photon cross-sections output= Off\n\n'],...
    [':Stop MC Transport Parameter:\n'],...
    ['#########################\n']);

fprintf(file,mctparameter);

fclose(file);

end