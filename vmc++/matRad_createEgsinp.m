%like createVmcInput except for making .egsinp files
%function matRad_createEgsinp(n,i,filebase)
%where n is beam number and i is beamlet number
%make it a function later when integrated into calcPhotonDoseEgs

%transform xyz coords using this matrix:
% T = [(cos(thetaT)*cos(thetaG)*sin(thetaC)-sin(thetaT)*cos(thetaC)) (cos(thetaT)*cos(thetaG)*sin(thetaC)+sin(thetaT)*cos(thetaC)) (cos(thetaT)*sin(thetaG));
%     sin(thetaG)*sin(thetaC) sin(thetaG)*cos(thetaC) -cos(thetaG);...
%     (-sin(thetaT)*cos(thetaG)*sin(thetaC) - cos(thetaT)*cos(thetaC)) (-sin(thetaT)*cos(thetaG)*sin(thetaC) + cos(thetaT)*cos(thetaC)) (-sin(thetaT)*sin(thetaG))];


%first make matrix that tells which beamlets (numbered according to the
%phsp files) are used in a given beam (is binary; 1 = used, 0 = not used):

xmin = -40; %these are in mm
xmax = 40;
ymin = -45;
ymax = 45;
n = 1;
allBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1),2);

for x = xmin:5:xmax
    for y = ymin:5:ymax
        allBeamlets(n,:) = [x y];
        n = n + 1
    end
end

whichBeamlets = zeros(((xmax-xmin)/5 + 1 )*( (ymax-ymin)/5 + 1)  ,length(stf));  %beamlets included (index in a given row) for each beam (columns)
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        indexOfPhsp = find( (allBeamlets(:,1) == stf(i).ray(j).rayPos_bev(1)) & (allBeamlets(:,2) == stf(i).ray(j).rayPos_bev(3)) );    %find index in reference vector that corresponds to this beamlet
        whichBeamlets(indexOfPhsp,i) = 1;
    end
end


%okay now make the phasespace files by looping over the beams and beamlets:

filebase = 'inputs';
phspfilebase = "dividedPhsp";
egsinpPath = 'egsinpFiles';

%just for first beam rn (hence n = 1) otherwise loop over n beams (here, 5)
% n=1;
% i = 1;
for n = 1:length(stf)
    for i = 1:length(whichBeamlets) %so it loops over ALL beamlets, not just those in the beam
        %so bc loops over all, I'll just make the egsinp file if the
        %relevant entry in the whichBeamlets matrix is 1. Else I keep going
        rec2 = 0;
        
        %to only include the beamlets in that beam, check if index of
        %whichBeamlets(row=i,column=n) is 1 or 0. If it's 1 we include it,
        %else we exclude it (so do nothing).
        %This means our beamlet numbers in the file names have gaps
        
        %rn the c script assumes file names are like:
        %filebaseBeamnBeamleti with n = beam number and i = beamlet number
        
        if(whichBeamlets(i,n) == 1)
            phantomName = "matRad_CT.egsphant"; % = strcat(ct.dicomInfo.PatientName.GivenName,ct.dicomInfo.PatientName.MiddleName,ct.dicomInfo.PatientName.FamilyName,"_CT",".egsphant");
            thisegsinpfile = strcat(filebase, 'Beam',num2str(n),'Beamlet',num2str(i), '.egsinp')
            sourcephspfile = strcat(phspfilebase,num2str(i), '.egsphsp1');

            %The angles in the matRad (DICOM) coordinate system:
            thetaC = 0; %madRad collimator angle
            thetaG = stf(n).gantryAngle;
            thetaT = stf(n).couchAngle;

            dsource = 50;   %distance in cm

            %angles in dosxyznrc coordinate system:
            dosxyzTheta = radtodeg(acos(-sin(thetaT)*sin(thetaG)));
            dosxyzPhi = radtodeg(atan( (-cos(thetaG))/(cos(thetaT)*sin(thetaG))));
            dosxyzPhiCol = 0; %thetaC- 90 + radtodeg(atan( (-sin(thetaT)*cos(thetaG))/(cos(thetaT))));


            %isocentre in matRad CT phantom coordinate system-> don't need to transform
            %divide by 10 to convert from mm to cm

            isocentre = [pln.propStf.isoCenter(n,1) pln.propStf.isoCenter(n,2) pln.propStf.isoCenter(n,3) ]/10;

            %valuess from 4th line/record in egsinp file:
            %ECUTIN,PCUTIN: Electron (total) and photon global cutoff energies in MeV.
            % SMAX: Dummy input
            rec4 = [0.7 0.01 0];   

            %zeroairdose,doseprint,MAX20: (binary 1/0):
            rec5 = [0 0 0];

            file = fopen(fullfile(egsinpPath, thisegsinpfile),'w');

            fprintf(file,"                                                                                 #!GUI1.0\n");
            fprintf(file, "%d\n",rec2);
            fprintf(file,"/data/data060/shussain/egsnrc/dosxyznrc/%s\n",phantomName);
            fprintf(file,"%f, %f, %f\n", rec4);
            fprintf(file,"%d, %d, %d\n", rec5);
            fprintf(file,"2, 2, %f, %f, %f, %f, %f, %f, %f, 0, 0, 0, 0, 0\n",isocentre(1),isocentre(2),isocentre(3),dosxyzTheta,dosxyzPhi,dsource,dosxyzPhiCol);
            fprintf(file,"2, 0, 0, 0, 0, 0, 0, 0\n");
            fprintf(file,"/data/data060/shussain/egsnrc/dosxyznrc/%s\n",sourcephspfile);

           
            nhistories = 100000;
            
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
    end
end
