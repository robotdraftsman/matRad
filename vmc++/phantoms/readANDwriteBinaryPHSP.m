%just a control file for the optimized read/write binary things so things
%are a bit easier to handle...

%so idea is call readBinary on a portion of the file. Give it the number of
%particles to skip, each of which consitutes 28 bytes. So the offset for
%reading things is 28 + numParticles2skip*28.
%And specify when to stop. Unfortinately I have it set to read in an entire
%thing in one go, so the "repeat" number is useless; it's basically just 1.
%Then we have the values of this chunk in the workspace.


clear;

statusFile = "matRad_running_status.txt";

%phspFile = '40x40.egsphsp1';
%phspFile = 'beamletTestphsp.egsphsp1';
%phspFile = '1x1_old.egsphsp1';
phspFile = '5x5_at_50cm.egsphsp1';
readThisMuch = inf;
numParticlesToSkip = 0;

% statusfid = fopen(statusFile,'w');
% fwrite(statusfid,"starting the phase space dividing process...");
% fclose(statusfid);
% scp_simple_put('tyr.physics.carleton.ca','shussain','Sands5rocks',statusFile,'/data/data060/shussain/egsnrc/dosxyznrc/','/Users/sakinahussain/Documents/GitHub/matRad/egsinpFiles/',statusFile);
fprintf("here we go with reading in the big 'un!\n");

[header phspData charges lastParticle numParticlesLeft] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip);

fprintf("Wow just finished reading in the large phsp file!\nNow to split it up and write the little guis.\n");


% statusfid = fopen(statusFile,'w');
% fwrite(statusfid,"finished reading in the phase space");
% fclose(statusfid);
scp_simple_put('tyr.physics.carleton.ca','shussain','Sands5rocks',statusFile,'/home/shussain/','/Users/sakinahussain/Documents/GitHub/matRad/',statusFile);


writeBinaryPHSP_optimized();

% statusfid = fopen(statusFile,'w');
% fwrite(statusfid,"finished writing the beamlet phase space files");
% fclose(statusfid);
% scp_simple_put('tyr.physics.carleton.ca','shussain','Sands5rocks',statusFile,'/home/shussain/','/Users/sakinahussain/Documents/GitHub/matRad/',statusFile);
fprintf("OK I just finished making the smaller phase space files. Next up, put them on the cluster!\n");

filebase = 'dividedPhsp36M';
for n = 1:361
    phspFileToTransfer = fullfile(strcat(filebase,num2str(n),'.egsphsp1'));
    scp_simple_put('tyr.physics.carleton.ca','shussain','Sands5rocks',phspFileToTransfer,'/data/data060/shussain/egsnrc/dosxyznrc/','/Users/sakinahussain/Documents/GitHub/matRad/beamletPHSPfiles',phspFileToTransfer); 
end

% fprintf("so that's a wrap!\n");
% statusfid = fopen(statusFile,'w');
% fwrite(statusfid,"All done! :D");
% fclose(statusfid);
scp_simple_put('tyr.physics.carleton.ca','shussain','Sands5rocks',statusFile,'/home/shussain/','/Users/sakinahussain/Documents/GitHub/matRad/',statusFile);