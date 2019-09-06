%% Spin Splitter
%% Converts m.out file in a set of oomf files organized according to a given structure configuration file
%% 
function spin_splitter(spin_file, mydirectory, conf_file, ti, tf, step, layer, output_dir)
%% Revision 1.4.1 2017-04-26 Introduced another parameter: output folder (output_dir)
cRel = '1.4.1 - 2017-04-26';
display(['Spin_Splitter ',cRel]);
current_dir = cd;

%% Configuration Block
disp('START Configuration');
structure = load_ascii(conf_file);
x = structure(1);
y = structure(2);
z = structure(3);

%% First Validation Block
if nargin<5
    error('Insufficient Params, only the Layer param is optional!');
elseif nargin < 6
    step = 1;
    % By default, compute all the layers (FLAG_COMPLETE = 1)
    layer = 0;
    % If not included use the root folder as the basis to create the
    % .\Spin_XX subfolder
    output_dir = mydirectory;
elseif nargin < 7
    % By default, compute all the layers (FLAG_COMPLETE = 1)
    layer = 0;
    % If not included use the root folder as the basis to create the
    % .\Spin_XX subfolder
    output_dir = mydirectory;
elseif nargin < 8
    % If not included use the root folder as the basis to create the
    % .\Spin_XX subfolder
    output_dir = mydirectory;
end

%% Second Validation Block
if tf <= ti
     error('Tf <= Ti');
elseif step <= 0
     error('Step must be positive!');
elseif ti == 0 || tf == 0
     error('Ti or Tf cannot be 0!');
elseif exist([mydirectory,filesep,spin_file])~=2 || exist(mydirectory)~=7|| exist(conf_file)~=2
     error('Files or Folder do not exist');
elseif layer>z || layer <0
     error(['Insert a layer with a non-negative value smaller than ',num2str(z)]);
end
%% End of Validation Block

%% Validation Block
display(['SpinFile[',spin_file,'], Folder[',mydirectory,'], OutputFolder[',output_dir,'], Conf_File[',conf_file,...
    '], Layer=',num2str(layer),', Ti=',num2str(ti),', Tf=',num2str(tf),', Step=',num2str(step)]);

% Create an identifier to main spin file (usually it is m.out)
fid_spin=fopen([mydirectory,filesep,spin_file],'r');

% Evaluate single of all-layers computation
if layer == 0
    % Whole structure Z
    FLAG_COMPLETE = 1;
    display(['Processing of all layers has been selected...']);
else
    % Single Layer Only SEE 1
    FLAG_COMPLETE = 0;
    display(['Processing of single layer [',num2str(layer),'] has been selected...']);
end

% Can support up to 999999 frames 
str_const='000000';
% Declare the output folders
Spin_single_layer=['Spin_of_Layer_',int2str(layer)];
Spin_all_layers = 'Spin_Complete';
% Check the flag
if FLAG_COMPLETE==0
    new_dir=[output_dir,filesep,Spin_single_layer];
    % Just a single layer
    temp_spin=zeros(x*y*1,3);
elseif FLAG_COMPLETE==1
    new_dir=[output_dir,filesep,Spin_all_layers];
    % All the layers
    temp_spin=zeros(x*y*z,3);
end
tspin = temp_spin;
mkdir(new_dir);
cd(new_dir); 
disp('END of Configuration');

%%
disp('Starting computation...');

% Start Preprocessing
time=tf-ti+1;

% Move to the beginning of the Ti selected (which is not necessarily Ti=1)
for k=1:x*y*z*ti
    fscanf(fid_spin, '%f %f %f', 3);
end

% Define the coarse time using the above defined step
cTime = 1:step:time;
% Process
for jj=1:time
    switch (FLAG_COMPLETE)
        %% Specific scenario
        case 0,            
            for kk=1:z
                nEL = x*y;
                tmp_ = fscanf(fid_spin, '%f %f %f', 3*nEL);
                temp_spin = (reshape(tmp_,[3,nEL]))';
                %% If we are analysing the right layer
                if kk==layer 
                    %% Save the oomf file only if the corresponding time index is 
                    %% equals to the coarse time as defined by the step parameter
                    if ~isempty(find(cTime == jj,1,'first'))
                        % Set the filename in a way that it is ALWAYS
                        % related to the whole time sequence (1,2,...N)
                        ind=num2str(jj+ti-1);
                        str_temp=str_const;
                        str_temp(end-length(ind)+1:end)=ind;
                        nome_file=['SL',str_temp,'.omf'];
                        write_oomf_file(nome_file,temp_spin,x,y,kk,FLAG_COMPLETE);
                        display(['Generated File ',nome_file,' iteration ',num2str(jj),' of ',num2str(time),' completed ',num2str(100*jj/time),'%']);
                    end
                end
            end
        %% Common case    
        case 1,
            nEL = x*y*z;
            tmp_ = fscanf(fid_spin, '%f %f %f', 3*nEL);
            temp_spin = (reshape(tmp_,[3,nEL]))';
            %% Save the oomf file only if the corresponding time index is 
            %% equals to the coarse time as defined by the step parameter
            if ~isempty(find(cTime == jj,1,'first'))
                % Imposta il nome del file
                ind=num2str(jj+ti-1);
                str_temp=str_const;
                str_temp(end-length(ind)+1:end)=ind;
                nome_file=['SL',str_temp,'.omf'];
                write_oomf_file(nome_file,temp_spin,x,y,z,FLAG_COMPLETE);
                display(['Generated File ',nome_file,' iteration ',num2str(jj),' of ',num2str(time),' completed ',num2str(100*jj/time),'%']);
            end
    end
end
fclose(fid_spin);

disp('... End of Computation');
% back to the original folder
cd(current_dir);
end

%% Many improvements
function write_oomf_file(nome_file,temp_spin,x,y,z,FLAG_COMPLETE)
if FLAG_COMPLETE==0
    z1=1;
elseif FLAG_COMPLETE==1
    z1=z;
end

fid=fopen(nome_file,'W');
fprintf(fid,'# OOMMF: rectangular mesh v0.99 \n');
fprintf(fid,'# Segment count: 1 \n');
fprintf(fid,'# Begin: Segment \n');
fprintf(fid,'# Begin: Header \n');
fprintf(fid,'# Title: model.omf \n');
fprintf(fid,'# Desc: \n');
fprintf(fid,'# meshunit: unknown \n');
fprintf(fid,'# xbase: 0 \n');
fprintf(fid,'# ybase: 0 \n');
fprintf(fid,'# zbase: 0 \n');
fprintf(fid,'# xstepsize: 1 \n');
fprintf(fid,'# ystepsize: 1 \n');
fprintf(fid,'# zstepsize: 1 \n');
fprintf(fid,'# xnodes: %d \n', x);
fprintf(fid,'# ynodes: %d \n', y);
fprintf(fid,'# znodes: %d \n', z1);
fprintf(fid,'# xmin: -0.5 \n');
fprintf(fid,'# ymin: -0.5 \n');
fprintf(fid,'# zmin: -0.5 \n');
fprintf(fid,'# xmax:   %f \n',x-0.5);
fprintf(fid,'# ymax:   %f \n',y-0.5);
fprintf(fid,'# zmax:  %f \n', z*FLAG_COMPLETE-0.5); % in order to have the correct configuration
fprintf(fid,'# valueunit: Ms \n');
fprintf(fid,'# valuemultiplier: 1 \n');
fprintf(fid,'# ValueRangeMaxMag:  1. \n');
fprintf(fid,'# ValueRangeMinMag:  1e-8 \n');
fprintf(fid,'# End: Header \n');
fprintf(fid,'# Begin: data text \n');

% Write data
if FLAG_COMPLETE==0
    k = x*y*1;
elseif FLAG_COMPLETE==1
    k = x*y*z;
end
A = temp_spin(1:k,:);
% Write the content ONCE forever!
fprintf(fid,'%f %f %f \n',A');
% Finalize file
fprintf(fid,'# End: data text \n');
fprintf(fid,'# End: segment \n');
% Close file
fclose(fid);
end
