%read in a binary file (the phase space file)
%EXCEPT DON'T USE THIS, use readBinaryPHSP_optimized.m instead! This one is
%less efficient and somewhat wrong. So get out!

%clear;

%phspFile = '40x40.egsphsp1';
%phspFile = 'beamletTestphsp.egsphsp1';
phspFile = '1x1_old.egsphsp1';
%phspFile = '5x5_at_50cm.egsphsp1';

% read the header:
 m = memmapfile(phspFile,...
     'Format',{'int8',[1,5],'mode';...
        'int32',[1,1],'NUM_PHSP_TOT';...
        'int32',[1,1],'PHOT_PHSP_TOT';...
        'single',[1,1],'EKMAX_PHSP_SHORT';...
        'single',[1,1],'EKMINE_PHSP_SHORT';...
        'single',[1,1],'NINC_PHSP_SHORT'},'Repeat',1);

    
%read everything after the header:
%but whether we read in zlast or not depends on the mode:
if (m.Data.mode == 'MODE0')
    m2 = memmapfile(phspFile,'Offset',28,...
        'Format',{'uint8',[1,4],'LATCH';...
            'single',[1,1],'ESHORT';...
            'single',[1,1],'X_PHSP_SHORT';...
            'single',[1,1],'Y_PHSP_SHORT';...
            'single',[1,1],'U_PHSP_SHORT';...
            'single',[1,1],'V_PHSP_SHORT';...
            'single',[1,1],'WT_PHSP_SHORT'},'Repeat',inf);
    "Running in MODE0: no ZLAST"
elseif (m.Data.mode == 'MODE2')
    m2 = memmapfile(phspFile,'Offset',28,...
        'Format',{'uint8',[1,4],'LATCH';...
            'single',[1,1],'ESHORT';...
            'single',[1,1],'X_PHSP_SHORT';...
            'single',[1,1],'Y_PHSP_SHORT';...
            'single',[1,1],'U_PHSP_SHORT';...
            'single',[1,1],'V_PHSP_SHORT';...
            'single',[1,1],'WT_PHSP_SHORT';...
            'single',[1,1],'ZLAST_PHSP_SHORT'},'Repeat',inf);
    "Running in MODE2: ZLAST included"
else
    "Uhh it's not mode0 or mode2 dude??\n"
end

strMode = char(m.Data.mode);

fprintf("memmapped everything; getting charges now...\n");

%now to get the required information out of latch
charges = zeros(length(m2.Data),1);
%for complement representation of signed integers:
for i = 1:length(m2.Data)
    getBits = dec2bin(m2.Data(i).LATCH(4));
    abyte = '00000000';
    abyte(9-length(getBits):8) = getBits(:);
    if(abyte(2) == '1' && abyte(3) == '0')
        charge = -1;
    elseif (abyte(2) == '0' && abyte(3) == '1')
        charge = 1;
    elseif (abyte(2) == '0' && abyte(3) == '0')
        charge = 0;
    else
        charge = 42;    %just so it's easily found
        fprintf("You have a 11 in binary for charge, dude!");
    end
    charges(i) = charge;
end

testArray = zeros(267,1);
testArray = m2.Data.ESHORT

%the index of charges matches the index of m2.Data from whence it came    

%A quick note: can calculate z angle cosine from the x and y ones, and its
%sign is the sign of the corresponding WT


% will have structure like:
%     '00000001'
%     '00000110'
%     '10000000'
%     '00010111'

% largest number is on left, i.e. 128 is 10000000

% bits 29 and 30 are charge -> so look at this last block
% in above block, the 11 near the end is the charge (goes from 0 to 31)
% now must determine what this is and what other charges are

% last block is mostly 23, but also 74, 55, 42, and 10

% 10: 00001010 -> 01 if doing later bits = smaller ones. else: 00
% 23: 00010111 -> 11 if doing later bits = smaller ones. else: 00
% 42: 00101010 -> 01 if doing later bits = smaller ones. else: 10
% 55: 00110111 -> 11 if doing later bits = smaller ones. else: 10
% 74: 01001010 -> 01 if doing later bits = smaller ones. else: 01\

% out of 365 particles, have 362 photons -> 23 should have charge 0
%using one's or two's complement gives right magnitude and sign

% orr is it the two that start second from left? What order does it count?

% !!!!
% assume that just take the binary number and convert. First thing is sign,
% second is whether it's 1 or 0. Just look up that stuff to see which, and
% of course make sure it's not backwards

%^^^^^^^^^^^^^^^^^^^^^


%mortran stuff to read in the header (~line 750):
% "Input parameters:"
% "{P1}=unit number"
% "{P2}=MODE0 or MODE2"
% "{P3}=total number of particles"
% "{P4}=no. of photons"
% "{P5}=max k.e. of particles"
% "{P6}=min. k.e. of electrons"
% "{P7}=no. of particles incident from original source"
% 
% READ({P1},REC=1,IOSTAT=IERR_PHSP){P2},NUM_PHSP_TOT,PHOT_PHSP_TOT,
% EKMAX_PHSP_SHORT,EKMINE_PHSP_SHORT,NINC_PHSP_SHORT;
% {P3}=NUM_PHSP_TOT;    --> an integer I think it defaults to integer*4
% {P4}=PHOT_PHSP_TOT;   --> so these guys are 4 bytes each = 32 bits
% {P5}=EKMAX_PHSP_SHORT; -> this is a real*4 so again 32-bit, or 4 bytes
% {P6}=EKMINE_PHSP_SHORT;   ->applies to this and the next as well
% {P7}=NINC_PHSP_SHORT;

% P2 is character*5 so 5 characters = 5 bytes
% so I think read it in as 5 chars, one at a time
% and then do the others

%so each variable is nonetheless 32 bits


%the actual data for each particle:

%interpretation of data has 2 cases: in one you record Z_last and in one
%you don't
%copy-pasted from mortran code: (this is what you do if first case):

%   READ({P2},REC={P3},IOSTAT=IERR_PHSP) {P6},ESHORT,X_PHSP_SHORT,Y_PHSP_SHORT,
% U_PHSP_SHORT,V_PHSP_SHORT,WT_PHSP_SHORT,ZLAST_PHSP_SHORT;

%just delete that last variable for it to be the second case

%now, sizes of these variables:
%    ESHORT,    "single precision E read from/written to phsp"
%    X_PHSP_SHORT, "single precision x read from phsp"
%    Y_PHSP_SHORT, "single precision y read from phsp"
%    U_PHSP_SHORT, "single precision u read from phsp"
%    V_PHSP_SHORT, "single precision v read from phsp"
%    WT_PHSP_SHORT, "single precision wt read from phsp"
%    ZLAST_PHSP_SHORT, "single precision zlast read from phsp"

%LATCH, E, X, Y, U, V, WT, (ZLAST)

%all single precision. that's 32 bits each


% oh my goodness. So these guys are just singles. This makes everything so
% much simpler! I can just read in the header using the repeat = 1 with the
% relevant format, then read the whole file and just delete the first few
% values, which correspond to the header.
% Then I can go through the remaining values and make a matrix with say
% rows being the different particles, and the columns being the particle
% data

% but don't forget to disassemble the LATCH and get the relevant stuff from
% it

% oh but make sure use the correct one between single and real*4








