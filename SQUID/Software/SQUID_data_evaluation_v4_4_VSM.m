%% Initialize and data input
clear all
close all
set(0,'defaultAxesFontSize',10)
set(0, 'DefaultLineLineWidth', 1);

% constants
uB = 9.274*10^-21;  % 1 Bohr magneton in emu
aLAO = 3.793*10^-10;    % lattice constant for LAO
aNGO = 3.858*10^-10;    % lattice constant for NGO 
aLSAT = 3.868*10^-10;    % lattice constant for LSAT
aLGO = 3.89*10^-10;    % lattice constant for LGO
aSTO = 3.905*10^-10;    % lattice constant for STO

aDSO = 3.948*10^-10;    % lattice constant for DSO

% select files containing MH and MT data
[filenames,pathname] = uigetfile('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/*.dat','MultiSelect','on');
path = textscan(pathname,'%s','delimiter','/') ;    % extract path
sample_name = char(path{1}(end));                   % extract sample name from folder name
sample_title = strrep(sample_name,'_','\_');        % replace '_' with '\_' for correct printing

% input for sample area, film thickness, substrate, mag. field for fit
prompt = {'Sample area in mm^2', 'Film thickness in u.c.', 'Substrate (LAO, NGO, LSAT, LGO, STO, DSO)', 'Min mag. field for fit in T', 'limit for point removal in %'};
input_title = 'Input';
input_size = [1 35];
default_input = {'','','','4', '0.4'};
input = inputdlg(prompt, input_title, input_size, default_input);

% write user input into corresponding variables
A = str2double(char(input(1)))*10^-6;   % m^2
thickness_uc = str2double(char(input(2)));
substrate = char(input(3));
Hmin = str2double(char(input(4)));
correction_limit = str2double(char(input(5)))/100;

% select substrate lattice constant corresponding to substrate input [m]
if strcmp(substrate,'LAO')
    a = aLAO;
elseif strcmp(substrate,'NGO')
    a = aNGO;
elseif strcmp(substrate,'LSAT')
    a = aLSAT;
elseif strcmp(substrate,'LGO')
    a = aLGO;
elseif strcmp(substrate,'STO')
    a = aSTO;
elseif strcmp(substrate,'DSO')
    a = aDSO;
else
    error('invalid substrate')
end

% calculate number of u.c. in film
Nuc = A/a^2*thickness_uc;

% sort files in MH and MT
% 'if ischar()' is necessary to distinguish the case where only one file is loaded. In this case 'filename' is a char and not a cell array.
if ischar(filenames)
    if contains(filenames,'_MH_')
        MH_filenames = filenames;
    end
    if contains(filenames,'_MT_')
        MT_filename = filenames;
    end
else
    if sum(contains(filenames,'_MH_')) ~= 0
        MH_filenames = filenames(contains(filenames,'_MH_'));
    end
    if sum(contains(filenames,'_MT_')) ~= 0
        MT_filename = filenames(contains(filenames,'_MT_'));
    end
end

% create a new folder to export the processed data to
mkdir([pathname '/data_export']);


%% Load substrate data

answer = questdlg('Use specific substrate data for correction?');
switch answer
    case 'Yes'
        sub = true;
        [sub_filenames,sub_pathname] = uigetfile('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/*.dat','MultiSelect','on');
        sub_path = textscan(sub_pathname,'%s','delimiter','/') ;    % extract path

        % sort file in MH and MT
        if ischar(sub_filenames)
            if contains(sub_filenames,'_MH_')
                MH_sub_filenames = sub_filenames;
            end
            if contains(sub_filenames,'_MT_')
                MT_sub_filename = sub_filenames;
            end
        else
            if sum(contains(sub_filenames,'_MH_')) ~= 0
                MH_sub_filenames = sub_filenames(contains(sub_filenames,'_MH_'));
            end
            if sum(contains(sub_filenames,'_MT_')) ~= 0
                MT_sub_filename = sub_filenames(contains(sub_filenames,'_MT_'));
            end
        end
    case 'No'
        sub = false;
        average_substrate = false;
        answer2 = questdlg('Remove average substrate contribution for MT? If No, constant value gets removed.');
        switch answer2
            case 'Yes'
                % load average data for correct substrate material
                average_substrate = true;
                if strcmp(substrate,'STO')
                    MT_average = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_average_MT.txt',2);
                elseif strcmp(substrate,'LSAT')
                    MT_average = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/LSAT001_average_MT_low.txt',2);
                elseif strcmp(substrate,'NGO')
                    MT_average = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/NGO110_average_MT.txt',2);
                elseif strcmp(substrate,'LGO')
                    MT_average = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/LGO110_average_MT.txt',2);
                end 
                temp_average = MT_average(:,1);     % temperature in K
                mag_average = MT_average(:,2);      % mag in mu_B / mm^2
                mag_average = mag_average*A*10^6/Nuc*2;   % mag in mu_B / f.u.
                [temp_average, sortIndex_temp_average] = sort(temp_average);
                mag_average = mag_average(sortIndex_temp_average);
                % shift curve to zero at high T
                mag345 = mag_average(temp_average>345);
                avg_mag_average_high_T = mean(mag345(~isnan(mag345)));
                mag_average = mag_average - avg_mag_average_high_T;
        end
end


answer3 = questdlg('Use average substrate data for MH? Only available for STO');
MH_subst_avg = false;
switch answer3
    case 'Yes'
        MH_subst_up_high = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_highT_average_MH_up.txt',2);
        H_subst_up_high = MH_subst_up_high(:,1);     % field in T
        mag_subst_up_high = MH_subst_up_high(:,2);      % mag in mu_B / mm^2
        mag_subst_up_high = mag_subst_up_high*A*10^6/Nuc*2;   % mag in mu_B / f.u.
        MH_subst_down_high = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_highT_average_MH_down.txt',2);
        H_subst_down_high = MH_subst_down_high(:,1);     % field in T
        mag_subst_down_high = MH_subst_down_high(:,2);      % mag in mu_B / mm^2
        mag_subst_down_high = mag_subst_down_high*A*10^6/Nuc*2;   % mag in mu_B / f.u.
        
        MH_subst_up_20 = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_20K_average_MH_up.txt',2);
        H_subst_up_20 = MH_subst_up_20(:,1);     % field in T
        mag_subst_up_20 = MH_subst_up_20(:,2);      % mag in mu_B / mm^2
        mag_subst_up_20 = mag_subst_up_20*A*10^6/Nuc*2;   % mag in mu_B / f.u.
        MH_subst_down_20 = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_20K_average_MH_down.txt',2);
        H_subst_down_20 = MH_subst_down_20(:,1);     % field in T
        mag_subst_down_20 = MH_subst_down_20(:,2);      % mag in mu_B / mm^2
        mag_subst_down_20 = mag_subst_down_20*A*10^6/Nuc*2;   % mag in mu_B / f.u.
        
        MH_subst_up_5 = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_5K_average_MH_up.txt',2);
        H_subst_up_5 = MH_subst_up_5(:,1);     % field in T
        mag_subst_up_5 = MH_subst_up_5(:,2);      % mag in mu_B / mm^2
        mag_subst_up_5 = mag_subst_up_5*A*10^6/Nuc*2;   % mag in mu_B / f.u.
        MH_subst_down_5 = csvread('/Users/annabelladrewanowski/Documents/UZH/Masterarbeit/Analysis/SQUID/substrate_reference_data/STO001_5K_average_MH_down.txt',2);
        H_subst_down_5 = MH_subst_down_5(:,1);     % field in T
        mag_subst_down_5 = MH_subst_down_5(:,2);      % mag in mu_B / mm^2
        mag_subst_down_5 = mag_subst_down_5*A*10^6/Nuc*2;   % mag in mu_B / f.u.
      
        MH_subst_avg = true;
        
        answer4 = questdlg('Correct the average substrate curve by subtracting high T sample data?');
        MH_subst_avg_highT_correction = false;
        switch answer4
            case 'Yes'
                MH_subst_avg_highT_correction = true;
        end
end




%% MT Evaluation

if exist('MT_filename','var')   % only executed if an MT file is loaded
    % read data
    MT_filename = char(MT_filename);
    MT = csvread([pathname MT_filename],28);
    temp = MT(:,3);
    mag = MT(:,5);
    [temp, sortIndex_temp] = sort(temp);
    mag = mag(sortIndex_temp);
    mag = medfilt1(mag,3);
    
    mag = mag/(uB*Nuc)*2;       % mag in mu_B / f.u.
    % shift curve to zero at high T
    % mag345 = mag(temp>345);
    mag345 = mag(temp>290);
    avg_mag_high_T = mean(mag345(~isnan(mag345)));
    % subtract substrate data
    mag = mag - avg_mag_high_T;
    
    % subtract data measured on the substrate if applicable
    substrate_check = false;
    if sub && exist('MT_sub_filename','var')    % check if specific MT substrate data is loaded

        MT_sub_filename = char(MT_sub_filename);    % load MT substrate data
        MT_sub = csvread([sub_pathname MT_sub_filename],28);
        temp_sub = MT_sub(:,3);
        mag_sub = MT_sub(:,5);
        [temp_sub, sortIndex_temp_sub] = sort(temp_sub);
        mag_sub = mag_sub(sortIndex_temp_sub);
        mag_sub = medfilt1(mag_sub,3);
        
        mag_sub = mag_sub/(uB*Nuc)*2;
        mag345 = mag_sub(temp_sub>345);
        avg_mag_sub_high_T = mean(mag345(~isnan(mag345)));
        mag_sub = mag_sub - avg_mag_sub_high_T;
        
        % interpolate substrate data so it has the same x scale as the
        % sample data (necessary for subtraction)
        mag_sub_interpol = interp1(temp_sub,mag_sub,temp);   
        
        % subtract substrate data from sample data
        mag_corrected = mag - mag_sub_interpol;
        
        substrate_check = true;
    
    % if no substrate data is available for this specific substrate
    else    
        % subtract average substrate data
        if average_substrate
            [temp_average, index] = unique(temp_average); 
            mag_average = mag_average(index);
            mag_average = medfilt1(mag_average,3);
            mag_average_interpol = interp1(temp_average,mag_average,temp);

            % subtract substrate data from sample data
            mag_corrected = mag - mag_average_interpol;
            
        % don't subtract substrate data
        else
            mag_corrected = mag;
        end
    end
        
    % Calculate Tc (still in development)
    
    temp_Tc = temp((20<temp)&(temp<330));
    temp_interp = min(temp_Tc):0.01:max(temp_Tc);
    mag_Tc = mag_corrected((20<temp)&(temp<330));              % cut high and low T (to get nice differentiation)
    mag_lp = lowpass(mag_Tc,0.01);   
    mag_interp = interp1(temp_Tc,mag_Tc,temp_interp);
    
    mag_diff = diff(mag_lp) ./ diff(temp_Tc);                    % first derivative of data
    mag_diff_lp = lowpass(mag_diff,0.01);                       % low pass filter first derivative
    mag_diff_interp = interp1(temp_Tc(2:end),mag_diff_lp,temp_interp);
    
    mag_diff2 = diff(mag_diff_lp) ./ diff(temp_Tc(2:end));       % second derivative of data
    mag_diff2_lp = lowpass(mag_diff2,0.01);                     % low pass filter second derivative
    mag_diff2_interp = interp1(temp_Tc(3:end),mag_diff2_lp,temp_interp);
    
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    zx = zci(mag_diff2_interp);
    
    for xx = length(zx):-1:1
        xpos = zx(xx);
        temp_infl = temp_interp(xpos);
        if temp_infl < 300 && mag_diff_interp(xpos) < -10^(-3)
            break
        end
    end
    mag_at_infl = mag_interp(xpos);
    slope = mag_diff_interp(xpos);
    
    Tc = (slope * temp_infl - mag_at_infl) / slope;
    
    
    % plot data
    f1 = figure;
    plot(temp, mag_corrected);
    xlabel('Temperature (K)')
    ylabel('Magnetization (\mu_B / f.u.)')
    if substrate_check
        title([sample_title ' M vs. T, substrate corrected'])
    else
        title([sample_title ' M vs. T'])
    end
    %annotation('textbox', [0.6, 0.8, 0.1, 0.1], 'String', ['T_C = ' num2str(Tc,4) ' K'])
    
    % save plot
    saveas(gcf,[pathname sample_name '_MT.png'])
    savefig([pathname sample_name '_MT.fig'])
    
    % export data
    output = [temp mag_corrected];
    fileID = fopen([pathname '/data_export/' MT_filename(1:end-4) '.txt'],'w');
    fprintf(fileID,'%-20s %-s\n','Temperature (K)','Magnetization (mu_B / f.u.)');
    fprintf(fileID,'%-20.4f %-.4f\n',output');
    fclose(fileID);
end


%% MH Evaluation

if exist('MH_filenames','var')      % only executed if an MH file is loaded
    
    % define size of plot window
    nplots = length(MH_filenames);
    if ischar(MH_filenames)         % if only one MH file is loaded this will be true because then the filename will be stored as char and not cell
        plotwidth = 1;
        plotheigth = 1;
        nplots = 1;
        f2 = figure('Position', [10 10 1200 800]);
        f3 = figure('Position', [10 10 1200 800]);
    elseif nplots < 7
        plotwidth = ceil(nplots/2);
        plotheigth = 2;
        f2 = figure('Position', [10 10 1440 960]);
        sgtitle([sample_title ' M vs. H'])
        f3 = figure('Position', [10 10 1440 960]);
    else
        plotwidth = ceil(nplots/3);
        plotheigth = 3;
        f2 = figure('Position', [10 10 1440 960]);
        sgtitle([sample_title ' M vs. H'])
        f3 = figure('Position', [10 10 1440 960]);
    end
    
    
    plotnumber = 1;
    
    % generate temperature list
    temp_list = [];
     for file = MH_filenames
         if ischar(MH_filenames)
            file = MH_filenames;  
         end
         MH = csvread([pathname char(file)],28);
         temp = round(MH(1,3));
         temp_list = [temp_list; temp];
     end
    
    % loop through MH files
    MH_filenames_reverse = flip(MH_filenames);
    for file = MH_filenames_reverse
        
        % the following if statement is necessary because in the case that only one
        % MH file is loaded 'MH_filenames' will be a char variable. in this case, looping
        % through 'MH_filenames' will loop through each of the filenames
        % characters. the filename has to be rewritten. in this case the
        % for loop is stoped after one iteration.
        if ischar(MH_filenames)
            file = MH_filenames;  
        end
        
        % read data
        MH = csvread([pathname char(file)],28);
        temp = round(MH(1,3));
        H = MH(:,4)/10000;  % field in Tesla
        M = MH(:,5);        % magnetization in emu
        
        % if two successive points are measured at the same field average them
        i = 1;
        while i < length(H)
            if abs(H(i+1) - H(i)) < 1e-3
                H(i) = mean([H(i) H(i+1)]);
                M(i) = mean([M(i) M(i+1)]);
                H(i+1) = [];
                M(i+1) = [];
            end
            i = i+1;
        end
        
        %%% At high fields (>3T) there can appear artifacts from the
        %%% measurement where the points will slightly deviate from the
        %%% straight line. These points have to be removed before the linear
        %%% fit is subtracted from the data. The data is split in two parts,
        %%% for pos and neg H field. For fields >3T the deviation of each point
        %%% from a fit done for fields >Hmin is calculated. If it's above a
        %%% threshold the point gets removed. Then a new fit for fields >Hmin
        %%% is made to then be substracted from the whole data range.

        % extract data used for first fit
        M_pos_fit = M(H>2);
        H_pos_fit = H(H>2);
        M_neg_fit = M(H<-2);
        H_neg_fit = H(H<-2);

        % index is needed to keep track of the order of the points. So the pos
        % and neg data can be recombined in the end
        index = 1:length(H);

        % separate data in pos and neg H field
        M_pos = M(H>0);
        H_pos = H(H>0);
        index_pos = index(H>0);
        M_neg = M(H<0);
        H_neg = H(H<0);
        index_neg = index(H<0);    

        % make fit for data correction
        fitpos = polyfit(H_pos_fit, M_pos_fit, 1);
        fitneg = polyfit(H_neg_fit, M_neg_fit, 1);

        % remove points that deviate more than 0.4 % from fit
        i = 1;
        while i <= length(H_pos)
            if abs(1-polyval(fitpos,H_pos(i))/M_pos(i)) > correction_limit && H_pos(i)>3 || abs(1-polyval(fitpos,H_pos(i))/M_pos(i)) > 100 && 1<H_pos(i) && H_pos(i)<3
                H_pos(i) = [];
                M_pos(i) = [];
                index_pos(i) = [];
            else
            i = i+1;
            end
        end

        i = 1;
        while i <= length(H_neg)
            if abs(1-polyval(fitneg,H_neg(i))/M_neg(i)) > correction_limit && H_neg(i)<-3 || abs(1-polyval(fitneg,H_neg(i))/M_neg(i)) > 100 && -1>H_neg(i) && H_neg(i)>-3
                H_neg(i) = [];
                M_neg(i) = [];
                index_neg(i) = [];
            else
            i = i+1;
            end
        end

        % recombine cleaned data
        H = [H_neg; H_pos];
        M = [M_neg; M_pos];
        index = [index_neg index_pos];
        [~, index] = sort(index);
        H = H(index);
        M = M(index);
        
        % medfilt for data between +- 1T
        M_up = [];
        H_up = [];
        M_down = [];
        H_down = [];
        for i = 1:length(H)-1
            if H(i+1) > H(i)
                H_up = [H_up H(i)];
                M_up = [M_up M(i)];
                if i == length(H)-1
                    last_point = 1;
                end
            elseif H(i+1) < H(i)
                H_down = [H_down H(i)];
                M_down = [M_down M(i)];
                if i == length(H)-1
                    last_point = 2;
                end
            end
        end
        if last_point == 1
            H_up = [H_up H(end)];
            M_up = [M_up M(end)];
        elseif last_point == 2
            H_down = [H_down H(end)];
            M_down = [M_down M(end)];
        end
        
        M_up_neg = M_up(H_up<-1);
        M_up_center = M_up(abs(H_up)<1);
        M_up_pos = M_up(H_up>1);
        M_up_center_start = M_up_center(1);
        M_up_center_end = M_up_center(end);
        M_up_center_filtered = medfilt1(M_up_center,3,'includenan','truncate');
        M_up_center_filtered(1) = M_up_center_start;
        M_up_center_filtered(end) = M_up_center_end;
        M_up = [M_up_neg M_up_center_filtered M_up_pos];
        
        M_down_neg = M_down(H_down<-1);
        M_down_center = M_down(abs(H_down)<1);
        M_down_pos = M_down(H_down>1);
        M_down_center_start = M_down_center(1);
        M_down_center_end = M_down_center(end);
        M_down_center_filtered = medfilt1(M_down_center,3,'includenan','truncate');
        M_down_center_filtered(1) = M_down_center_start;
        M_down_center_filtered(end) = M_down_center_end;
        M_down = [M_down_pos M_down_center_filtered M_down_neg];
        
        M = [M_down M_up]';
        
        
        
        % subtract data measured on the specific substrate if applicable
        substrate_check = false;
        if exist('MH_sub_filenames','var') % check if MH substrate data is loaded
            
            % Check if for the temperature of the current iteration there
            % is substrate data loaded. If yes, select corresponding file
            % name
            if ischar(MH_sub_filenames)
                if contains(MH_sub_filenames,[num2str(temp) 'K'])
                    MH_sub_filename_it = MH_sub_filenames;
                    substrate_check = true;
                end
            else
                if sum(contains(MH_sub_filenames,[num2str(temp) 'K'])) ~= 0
                    MH_sub_filename_it = MH_sub_filenames(contains(MH_sub_filenames,[num2str(temp) 'K']));
                    substrate_check = true;
                end
            end
            
            % if substrate data for the current iteration is loaded, start
            % removal process
            if substrate_check
                % load substrate data
                MH_sub_filename_it = char(MH_sub_filename_it);
                MH_sub = csvread([sub_pathname MH_sub_filename_it],28);
                H_sub = MH_sub(:,4)/10000;  % field in Tesla
                M_sub = MH_sub(:,5);        % mag in emu
                %M_sub = medfilt1(M_sub,3);
                
                % if two successive points are measured at the same field average them
                i = 1;
                while i < length(H_sub)
                    if abs(H_sub(i+1) - H_sub(i)) < 1e-3
                        H_sub(i) = mean([H_sub(i) H_sub(i+1)]);
                        M_sub(i) = mean([M_sub(i) M_sub(i+1)]);
                        H_sub(i+1) = [];
                        M_sub(i+1) = [];
                    end
                    i = i+1;
                end
                
                
                % initalize vectors
                M_up = [];
                H_up = [];
                M_down = [];
                H_down = [];
                
                M_sub_up = [];
                H_sub_up = [];
                M_sub_down = [];
                H_sub_down = [];
                              
                
                % devide sample data into parts where field is ramped up or down
                for i = 1:length(H)-1
                    if H(i+1) > H(i)
                        M_up = [M_up M(i)];
                        H_up = [H_up H(i)];
                    elseif H(i+1) < H(i)
                        M_down = [M_down M(i)];
                        H_down = [H_down H(i)];
                    end
                end
                % put last value in same vector as second but last
                if H(i+1) > H(i)
                    M_up = [M_up M(i+1)];
                    H_up = [H_up H(i+1)];
                elseif H(i+1) < H(i)
                    M_down = [M_down M(i+1)];
                    H_down = [H_down H(i+1)];
                end
                
                % sort sample data for ascending field
                [sorted_H_up, sortIndex_H_up] = sort(H_up);
                M_up = M_up(sortIndex_H_up);
                H_up = sorted_H_up;
                
                [sorted_H_down, sortIndex_H_down] = sort(H_down);
                M_down = M_down(sortIndex_H_down);
                H_down = sorted_H_down;
                
                % devide substrate data into parts where field is ramped up or down
                for i = 1:length(H_sub)-1
                    if H_sub(i+1) > H_sub(i)
                        M_sub_up = [M_sub_up M_sub(i)];
                        H_sub_up = [H_sub_up H_sub(i)];
                    elseif H_sub(i+1) < H_sub(i)
                        M_sub_down = [M_sub_down M_sub(i)];
                        H_sub_down = [H_sub_down H_sub(i)];
                    end
                end
                if H_sub(i+1) > H_sub(i)
                    M_sub_up = [M_sub_up M_sub(i+1)];
                    H_sub_up = [H_sub_up H_sub(i+1)];
                elseif H_sub(i+1) < H_sub(i)
                    M_sub_down = [M_sub_down M_sub(i+1)];
                    H_sub_down = [H_sub_down H_sub(i+1)];
                end
                % put last value in same vector as second but last
                
                
                
%                 % remove duplicates
%                 [H_up, index] = unique(H_up); 
%                 M_up = M_up(index);
%                 
%                 [H_down, index] = unique(H_down); 
%                 M_down = M_down(index);
%                 
%                 [H_sub_up, index] = unique(H_sub_up); 
%                 M_sub_up = M_sub_up(index);
%                 
%                 [H_sub_down, index] = unique(H_sub_down); 
%                 M_sub_down = M_sub_down(index);
                
                
                % sort substrate data for ascending field
                [sorted_H_sub_up, sortIndex_H_sub_up] = sort(H_sub_up);
                M_sub_up = M_sub_up(sortIndex_H_sub_up);
                H_sub_up = sorted_H_sub_up;
                
                [sorted_H_sub_down, sortIndex_H_sub_down] = sort(H_sub_down);
                M_sub_down = M_sub_down(sortIndex_H_sub_down);
                H_sub_down = sorted_H_sub_down;
                
                % interpolate substrate data so that the x scaling matches
                % that of the sample data
                M_sub_up_interpol = interp1(H_sub_up,M_sub_up,H_up,'linear','extrap');
                M_sub_down_interpol = interp1(H_sub_down,M_sub_down,H_down,'linear','extrap');
                
                % subtract substrate data from sample data
                M_up_corrected = M_up - M_sub_up_interpol;
                M_down_corrected = M_down - M_sub_down_interpol;
                
                % combine up and down data into one vector
                M_corrected = [flip(M_down_corrected) M_up_corrected];
                H = [flip(H_down) H_up];
                
                % convert emu to uB / f.u.
                M_corrected_specific = M_corrected/(uB*Nuc)*2;
                MM = M_corrected_specific;  % used for determination of sat. mag, coercive field, sat. mag.
            end 
        end
        
        
        % if no MH substrate data is loaded for the temperature of the
        % current iteration, perform a linear fit for |H| > |Hmin|
        if ~substrate_check
            % interpolate linear part of the MH curve starting at Hmin
            M = M/(uB*Nuc)*2;  %convert M to uB / f.u.
            %M = M/uB/A/10^6;       %convert M to uB / mm^2
            Mpos = M(H>Hmin);
            Hpos = H(H>Hmin);
            Mneg = M(H<-Hmin);
            Hneg = H(H<-Hmin);

            fitpos = polyfit(Hpos, Mpos, 1);
            fitneg = polyfit(Hneg, Mneg, 1);
            slope = mean([fitpos(1), fitneg(1)]);
            M_corrected_linear = M - H*slope;    % subtract linear fit from the data
            MM = M_corrected_linear;  % used for determination of sat. mag, coercive field, sat. mag.
            
            % subtract average subst data if available
            if  strcmp(substrate,'STO') && MH_subst_avg
                if temp > 20
                    H_subst_up = H_subst_up_high;
                    mag_subst_up = mag_subst_up_high;
                    H_subst_down = H_subst_down_high;
                    mag_subst_down = mag_subst_down_high;
                elseif temp == 20
                    H_subst_up = H_subst_up_20;
                    mag_subst_up = mag_subst_up_20;
                    H_subst_down = H_subst_down_20;
                    mag_subst_down = mag_subst_down_20;
                elseif temp == 5
                    H_subst_up = H_subst_up_5;
                    mag_subst_up = mag_subst_up_5;
                    H_subst_down = H_subst_down_5;
                    mag_subst_down = mag_subst_down_5;
                end
                
                H_up = [];
                M_up = [];
                H_down = [];
                M_down = [];
                for i = 1:length(H)-1
                    if H(i+1) > H(i)
                        H_up = [H_up H(i)];
                        M_up = [M_up M_corrected_linear(i)];
                        if i == length(H)-1
                            last_point = 1;
                        end
                    elseif H(i+1) < H(i)
                        H_down = [H_down H(i)];
                        M_down = [M_down M_corrected_linear(i)];
                        if i == length(H)-1
                            last_point = 2;
                        end
                    end
                end
                if last_point == 1
                    H_up = [H_up H(end)];
                    M_up = [M_up M_corrected_linear(end)];
                elseif last_point == 2
                    H_down = [H_down H(end)];
                    M_down = [M_down M_corrected_linear(end)];
                end
                
                % check if high T data can be used to improve average substrate data
                if plotnumber == 1 && MH_subst_avg_highT_correction
                    max_temp = max(temp_list);
                    if temp ~= max_temp
                        error('last MH file is not the one with the highest temperature')
                    end
                    if max_temp >= 100
                        msgbox(strcat('average MH loops are corrected for the ', num2str(max_temp), ' K value of the sample'))
                        M_highT_pos_field_avg = mean(M_corrected_linear(H>5.4));
                        M_highT_neg_field_avg = abs(mean(M_corrected_linear(H<-5.4)));
                        M_highT_avg = mean([M_highT_pos_field_avg M_highT_neg_field_avg]);
                        
                        M_subst_up_highT_pos_field_avg = mean(mag_subst_up(H_subst_up>5.4));
                        M_subst_down_highT_pos_field_avg = mean(mag_subst_down(H_subst_down>5.4));
                        M_subst_up_highT_neg_field_avg = abs(mean(mag_subst_up(H_subst_up<-5.4)));
                        M_subst_down_highT_neg_field_avg = abs(mean(mag_subst_down(H_subst_down<-5.4)));
                        M_subst_highT_avg = mean([M_subst_up_highT_pos_field_avg; M_subst_down_highT_pos_field_avg; M_subst_up_highT_neg_field_avg; M_subst_down_highT_neg_field_avg]);
                        
                        subst_correction_factor = M_highT_avg/M_subst_highT_avg;
                    else
                        msgbox('No MH data >= 100 K. No correction of avg. substrate data performed.')
                    end
                end
                
                if max_temp >= 100  && MH_subst_avg_highT_correction
                    mag_subst_up = mag_subst_up * subst_correction_factor;
                    mag_subst_down = mag_subst_down * subst_correction_factor;
                end                    
                    
                 mag_subst_up_interpol = interp1(H_subst_up,mag_subst_up,H_up);
                 M_up_corrected = M_up -  mag_subst_up_interpol;
                 mag_subst_down_interpol = interp1(H_subst_down,mag_subst_down,H_down);
                 M_down_corrected = M_down -  mag_subst_down_interpol;
                 
                 
                 
                 M_corrected_avg = [M_down_corrected M_up_corrected]';
                 MM = M_corrected_avg;  % used for determination of sat. mag, coercive field, sat. mag.
            end
        
        end
        
        % find rem. mag., corecive field and sat. mag
        if plotnumber==1
            coercive_field = [];
            rem_mag = [];
            sat_mag = [];
        end
        H_rounded = round(H,2);
        % rem mag
        rem_mag = [rem_mag; mean(abs(MM(H_rounded==0)))];
        % sat mag
        sat_mag = [sat_mag; mean(abs(MM(H_rounded==-5 | H_rounded==5)))];
        % coercive field
        H_up = [];
        M_up = [];
        H_down = [];
        M_down = [];
        for i = 1:length(H)-1
            if H(i+1) > H(i)
                H_up = [H_up H(i)];
                M_up = [M_up MM(i)];
                if i == length(H)-1
                    last_point = 1;
                end
            elseif H(i+1) < H(i)
                H_down = [H_down H(i)];
                M_down = [M_down MM(i)];
                if i == length(H)-1
                    last_point = 2;
                end
            end
        end
        if last_point == 1
            H_up = [H_up H(end)];
            M_up = [M_up MM(end)];
        elseif last_point == 2
            H_down = [H_down H(end)];
            M_down = [M_down MM(end)];
        end
        
        H_highres_up = min(H_up):0.001:max(H_up);
        H_highres_down = max(H_down):-0.001:min(H_down);
        MM_up_interpol = interp1(H_up,M_up,H_highres_up);
        MM_down_interpol = interp1(H_down,M_down,H_highres_down);
        
        [MM_min_up, MM_min_up_position] = min(abs(MM_up_interpol));
        [MM_min_down, MM_min_down_position] = min(abs(MM_down_interpol));
        coercive_field = [coercive_field; mean([abs(H_highres_up(MM_min_up_position)), abs(H_highres_down(MM_min_down_position))])];

        % plot data
        figure(f2)
        hold on
        subplot(plotheigth,plotwidth,plotnumber)
        if substrate_check
            plot(H,M_corrected_specific)
        elseif strcmp(substrate,'STO') && MH_subst_avg
            hold on
            plot(H,M_corrected_avg)
            %plot(H_up,mag_subst_up_interpol,'r')
            %plot(H_down,mag_subst_down_interpol,'r')
        else
            plot(H,M_corrected_linear)
        end
        axis([-8 8 -6 6])
        xlabel('Magnetic Field (T)')
        ylabel('Magnetization (/mu_B / f.u.)')
        %ylabel('Magnetization (/mu_B / mm^2)')
        if nplots > 1
            if substrate_check
                title([num2str(temp) ' K, specific substrate corrected'])
            elseif strcmp(substrate,'STO') && MH_subst_avg
                title([num2str(temp) ' K, average substrate corrected'])
            else
                title([num2str(temp) ' K, linear correction'])
            end
        else
            if substrate_check
                sgtitle([sample_title ' M vs. H ' num2str(temp) ' K, specific substrate corrected'])
            elseif strcmp(substrate,'STO') && MH_subst_avg
                title([num2str(temp) ' K, average substrate corrected'])
            else
                sgtitle([sample_title ' M vs. H ' num2str(temp) ' K, linear correction'])
            end
        end
        
        % plot uncorrected data together with the linear fits (ommited if substrate data is used for correction in the current iteration)
        if ~substrate_check
            figure(f3)
            subplot(plotheigth,plotwidth,plotnumber)
            hold on
            plot(H,M)
            plot(-8:8,polyval(fitpos,-8:8),'g--')
            plot(-8:8,polyval(fitneg,-8:8),'g--')
            plot(-8:8,[-8:8]*slope,'r')
            xlabel('Magnetic Field (T)')
            ylabel('Magnetization (/mu_B / f.u.)')
            hold off
        end
        
        % export data
       
        if substrate_check
            export_name = char(file);
            output = [H; M_corrected_specific];
            fileID = fopen([pathname '/data_export/' export_name(1:end-4) '_specific_subst_correction.txt'],'w');
            fprintf(fileID,'%-20s %-s\n','Magnetic field (T)','Magnetization (mu_B / f.u.)');
            fprintf(fileID,'%-20.4f %-.4f\n',output');
            fclose(fileID);
        elseif strcmp(substrate,'STO') && MH_subst_avg
            export_name = char(file);
            output = [H M_corrected_avg];
            fileID = fopen([pathname '/data_export/' export_name(1:end-4) '_avg_subst_correction.txt'],'w');
            fprintf(fileID,'%-20s %-s\n','Magnetic field (T)','Magnetization (mu_B / f.u.)');
            fprintf(fileID,'%-20.4f %-.4f\n',output');
            fclose(fileID);
        else
            export_name = char(file);
            output = [H M_corrected_linear];
            fileID = fopen([pathname '/data_export/' export_name(1:end-4) '_lin_subst_correction.txt'],'w');
            fprintf(fileID,'%-20s %-s\n','Magnetic field (T)','Magnetization (mu_B / f.u.)');
            fprintf(fileID,'%-20.4f %-.4f\n',output');
            fclose(fileID);
        end


        plotnumber = plotnumber + 1;

        % break the for loop if only one MH file is loaded
        if ischar(MH_filenames)
            break
        end

    end
    
    % save plot
    figure(f2)
    saveas(gcf,[pathname sample_name '_MH.png'])
    savefig([pathname sample_name '_MH.fig'])
    
    % export sat. mag, coercive field, rem. mag
    output = [temp_list sat_mag rem_mag coercive_field];
    fileID = fopen([pathname '/data_export/satmag_remmag_coercivefield.txt'],'w');
    fprintf(fileID,'%-20s %-40s %-40s %-s \n','Temperature (K)','Saturation Magnetization (mu_B / f.u.)', 'Remanent Magentization (mu_B / f.u.)', 'Coercive Field (T)');
    fprintf(fileID,'%-20.4f %-40.4f %-40.4f %-.4f \n',output');
    fclose(fileID);
    
end


