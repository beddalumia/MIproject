%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOW TO USE?
%
%  - Put this script, together with an input-file for your driver in a path
%    > This path will contain directories for all the U values you set-up.
%  - Set-up the name of your driver program (without .f90 extension)
%    > e.g. driver = 'cdn_hm_2dsquare'; 
%  - Set Nbath to your desire: --> you will get a fixed-Nbath linear span
%  - Adjust Umin and Umax to your desire --> U \in [Umin, Umax]
%  - Adjust Uold to catch a 'restart-folder' in the path [!applies -> -1]
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Run everything with $ matlab -batch MITline_autostep
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_cdmft.txt which is just a mirror of the CDMFT output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole CDMFT
%  - Also a additional output files in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%    > possibly some error-flag files in the format 'ERROR_U=%f'
%      if you see them you *may* have to discard the corresponding folder
%      > check it for convergence! (look at LOG_cdmft.txt)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'cdn_hm_2dsquare';	doMPI = true;

% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ulist = fopen('U_list.txt','a');

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbath  = 1;                      % Input Bath size 
Umin = 0; Umax = 10;             % Input Hubbard 
Wmix = 0.3;		 	               % Input Self-Mixing

Ustep = [0.5, 0.25, 0.1, 0.05, 0.01]; % To be auto-determined...hopefully!
NUstep = length(Ustep);

notConvFlag = false;		         % Convergence-fail *flag*
notConvCount = 0;		            % Convergence-fail *counter*
notConvThreshold = NUstep-1;     % Maximum #{times} we accept CDMFT to fail

U = Umin; Uold = -1;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);      % ------------------------------------
if isfolder(oldDIR)                  % If it exist a "previous" folder: 
restartpack = [oldDIR,'/*.restart']; % Copy all the restart files from the
copyfile(restartpack);               % last cdmft evaluation...
end                                  % ------------------------------------

copyfile ../input*             % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';				% Control of MPI
else						      % boolean flag...
mpi = [];
end
HUBBARD =sprintf(' uloc=%f',U);		      % OVERRIDE
NBATH =sprintf(' nbath=%d',Nbath);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
outLOG = ' > LOG_cdmft.txt';
cdmft_setup = [mpi,driver,HUBBARD,NBATH,MIXING];
cdmft_ed_call = [cdmft_setup,' Nloop=100 DM_FLAG=F ',outLOG];
tic
system(cdmft_ed_call);				         % Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');
fprintf(file_id,'%f\n', chrono);		      % Write on time-log
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) CDMFT LOOP
if isfile('ERROR.README')
    notConvFlag = true;
    notConvCount = notConvCount + 1;
    movefile('ERROR.README',sprintf('../ERROR_U=%f',U));
else
    fprintf(Ulist,'%f\n', U);	               % Write on U-log
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAST LOOP [Evaluates the RDMs]
outLOG_last = ' > LOG_lastloop.txt';
last_loop_call = [cdmft_setup,' Nloop=1 DM_FLAG=T ',outLOG_last];
% .restart -> .restart_original
restartpack = dir([pwd, '/*.restart']);
for i = 1:numel(restartpack)
   file = fullfile(pwd, restartpack(i).name);
   [tempDir, tempFile] = fileparts(file); 
   copyfile(file, fullfile(tempDir, [tempFile, '.restart_original']));
   delete(file);
end
% .used -> .restart
usedpack = dir([pwd, '/*.used']);
for i = 1:numel(usedpack)
   file = fullfile(pwd, usedpack(i).name);
   [tempDir, tempFile] = fileparts(file); 
   copyfile(file, fullfile(tempDir, [tempFile, '.restart']));
end
% Last LOOP call
system(last_loop_call);
delete('ERROR.README');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd ..                                        % Exit the U-folder


if notConvCount > notConvThreshold
   error('CDMFT not converged: phase-span stops now!');         
end 

if notConvFlag == true
   U = Uold; 			% if nonconverged we don't want to update 
   notConvFlag = false;	% > but we want to reset the flag(!)
else
   Uold = U; 			% if converged we update Uold and proceed to...         
end

U = U + Ustep(notConvCount+1); % ...Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);



