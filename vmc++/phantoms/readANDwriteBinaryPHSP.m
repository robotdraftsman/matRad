%just a control file for the optimized read/write binary things so things
%are a bit easier to handle...

%so idea is call readBinary on a portion of the file. Give it the number of
%particles to skip, each of which consitutes 28 bytes. So the offset for
%reading things is 28 + numParticles2skip*28.
%And specify when to stop. Unfortinately I have it set to read in an entire
%thing in one go, so the "repeat" number is useless; it's basically just 1.
%Then we have the values of this chunk in the workspace.

clear;
%phspFile = '40x40.egsphsp1';
%phspFile = 'beamletTestphsp.egsphsp1';
%phspFile = '1x1_old.egsphsp1';
phspFile = '5x5_at_50cm.egsphsp1';
readThisMuch = inf;
numParticlesToSkip = 0;

[header phspData charges lastParticle numParticlesLeft] = readBinaryPHSP_optimized(phspFile,readThisMuch,numParticlesToSkip);