% just some scp stuff



%move egsinp file from here to the cluster:
filebase = 'inputs';
for n = 1:length(stf)
   for i = 1:stf(n).numOfRays
       egsinp_file = strcat(filebase, 'Beam',num2str(n),'Beamlet',num2str(i), '.egsinp');
       scp_simple_put('tyr.physics.carleton.ca','[username]','[password]',egsinp_file,'/data/data060/shussain/egsnrc/dosxyznrc/','/Users/sakinahussain/Documents/GitHub/matRad/egsinpFiles/',egsinp_file);
   end
end