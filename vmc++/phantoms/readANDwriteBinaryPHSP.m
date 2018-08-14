%just a control file for the optimized read/write binary things so things
%are a bit easier to handle...

%so idea is call readBinary on a portion of the file. Give it the number of
%particles to skip, each of which consitutes 28 bytes. So the offset for
%reading things is 28 + numParticles2skip*28.
%And specify when to stop. Unfortinately I have it set to read in an entire
%thing in one go, so the "repeat" number is useless; it's basically just 1.
%Then we have the values of this chunk in the workspace.

password = input("what's your cluster password?\n");

% clear;
% 
% 
% %phspFile = '40x40.egsphsp1';
% %phspFile = 'dividedPhsp36M201.egsphsp1';
% %phspFile = '1x1_old.egsphsp1';
% phspFile = '5x5_at_50cm.egsphsp1';
% readThisMuch = inf;
% numParticlesToSkip = 0;
% 
% fprintf("here we go with reading in the big 'un!\n");
% example_txtmsg('phsp file stuff','beginning');
% try
%   [header phspData charges lastParticle numParticlesLeft] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip);
% catch
%   fprintf("ran into an error when reading in the large file :(\n");
%   example_txtmsg('phsp file stuff','ran into an error when reading in the large file :(');
%   error('something went wrong in readBinaryPHSP_optimized...');
% end
% 
% fprintf("Wow just finished reading in the large phsp file!\nNow to split it up and write the little guis.\n");

example_txtmsg('phsp file stuff','got the big file. Now split it up...');

fclose all;
try
    BeamHeaders = writeBinaryPHSP_optimized(phspData,charges,header);
catch
    fprintf("ran into an error when reading in the large file :(\n");
    example_txtmsg('phsp file stuff','ran into an error when writing the smaller files :(');
    error('something went wrong in writeBinaryPHSP_optimized...');
end

example_txtmsg('phsp file stuff','finished dividing it up. Now transferring to cluster...');

fprintf("OK I just finished making the smaller phase space files. Next up, put them on the cluster!\n");

try
    filebase = 'dividedPhsp36M';
    for n = 1:361
        phspFileToTransfer = fullfile(strcat(filebase,num2str(n),'.egsphsp1'));
        scp_simple_put('tyr.physics.carleton.ca','shussain',password,phspFileToTransfer,'/data/data060/shussain/egsnrc/dosxyznrc/','/Users/sakinahussain/Documents/GitHub/matRad/beamletPHSPfiles',phspFileToTransfer); 
    end
catch
    fprintf("ran into an error when transferring the beamlet phsp files :(\n");
    example_txtmsg('phsp file stuff','ran into an error when writing the smaller files :(');
    error('something went wrong when transferring ot the cluster...');
end

example_txtmsg('phsp file stuff','done! Now you can run dosxyz :)');

fprintf("so that's a wrap!\n");


