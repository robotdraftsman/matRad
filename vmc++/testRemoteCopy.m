%practise with remotecopy

%remotecopy -local copyMePls.txt -to /home/shussain -remotehost thor -remoteplatform { unix }

[status, cmdout] = system("cd");
cmdout
%"ssh -Y shussain@thor")