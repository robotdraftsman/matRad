 

%get number of particles in each beamlet

stuffinBeamlet = zeros(361,4);
filebase = 'dividedPhsp';
filepath = '/Users/sakinahussain/Documents/GitHub/beamletPHSPfiles-OLD';
%for i = 1:361
i = 156;
    phspFile = fullfile(filepath,strcat(filebase,num2str(i),'.egsphsp1'));
    m = memmapfile(phspFile,...
         'Format',{'int8',[1,5],'mode';...
            'int32',[1,1],'NUM_PHSP_TOT';...
            'int32',[1,1],'PHOT_PHSP_TOT';...
            'single',[1,1],'EKMAX_PHSP_SHORT';...
            'single',[1,1],'EKMINE_PHSP_SHORT';...
            'single',[1,1],'NINC_PHSP_SHORT'},'Repeat',1);
    
        stuffinBeamlet(i,1) = m.Data.EKMAX_PHSP_SHORT;
        stuffinBeamlet(i,2) = m.Data.EKMINE_PHSP_SHORT;
        stuffinBeamlet(i,3) = m.Data.NUM_PHSP_TOT;
        stuffinBeamlet(i,4) = m.Data.PHOT_PHSP_TOT;
  % clear m;
  
%end

% phspFile = '1x1.egsphsp1';
% [header2 phspData2 charges2 lastParticle2 numParticlesLeft2] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip);
% "ok read in the old phsp file"

% BeamHeaders2 = writeBinaryPHSP_optimized(phspData2,charges2,header2);
