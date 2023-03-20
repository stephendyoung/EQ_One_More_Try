%% Generate .sac files for python analysis 
%custom functions used: IsWdwAvailable, InfraData, OKutc2ct,DateFolderCk,
%chronilogical

year1='2018';
month1='7';
day1='5';
hour1='2';
min1='0';
wdw='2';

micnumber=datagetter(year1,month1,day1, hour1,min1, wdw, 0) ;



function micnum = datagetter(year1,month1,day1, hour1,min1, wdw, iswdwavailoverride) 
    %clearvars -except windcat

    tstart=tic;

    %% controls
    cmpt = 1;                       %1=EFPL server; 2=BRE laptop



                           %window increment in hours
    %dates1=datetime(2018,04,09,22,0,0);
    filename=strcat(year1,month1,day1,hour1,wdw);
    sacmic=ones(3,1);
    disp(strcat('filename= ',filename))



    wdw1=str2num(wdw)
    year1=str2num(year1)
    month1=str2num(month1)
    day1=str2num(day1)
    hour1=str2num(hour1)
    min1=str2num(min1)
    if cmpt == 1
        Drive = '';
        MainFolder = '/run/user/1000/gvfs/smb-share:server=efpl-fs.ceat.okstate.edu,share=projects/Infrasound/Raw Data/Matlab/IA1/';
    end

    location = strcat(MainFolder);
    YEARS = dir(strcat(MainFolder));

    event = [year1,month1,day1,hour1,min1,0];
    %% identify start and end times for data files
    %START,END,EVENT = [year,month,day,hour,minute,second]
    %T_start,T_end,T_event = matlab date number format

    if event(1,1) > 2017
        START = event;
        END = event;
        END(4) = END(4)+wdw1;
        EVENT = START;
    else

        START = OKutc2ct(event(1,1:3),event(1,4:6));
        event_end = event;
        event_end(4) = event_end(4) + wdw1;
        END = OKutc2ct(event_end(1,1:3),event_end(1,4:6));
        EVENT = START;
    end
     %disp(location)
     Dck = IsWdwAvailable(START,wdw1,location);        %1=data is available on that day; 0 = not available
     
    if Dck == 0
        disp('data no available')
        Dck = iswdwavailoverride
    elseif Dck == 1
        
            T_start = datenum(START);           %start time (used to save files; CT or UTC)
            T_end = datenum(END);               %end time (used to save files; CT or UTC)
            T_event = datenum(START);           %event time (used to save files; CT or UTC)

            %% make time trace of data range of interest
            clear DATA RR
            [DATA,Nmic,fs] = InfraData(START,END,EVENT,location);

            %% identify mic data of interest
            
            for imic = 1:int32(Nmic)

                    MIC = strcat('mic',num2str(imic));
                    
                    pA = DATA.(MIC)(:,2);           %sound pressure level (Pa)
                    ts = DATA.(MIC)(:,1)*60;        %relative time in seconds
                    if imic == 1
                        S(:,4)=ts(:);
                    end
                    S(:,imic)=pA(:);
                    
                    maxmic=max(S(:,imic));
                    if maxmic==0
                        sacmic(imic)=0
                    end
                    clear MIC pA ts

            end
            S(isnan(S))=0;
    %         end
        
            

        
    end

    %po=S(1,4);
    %hi=time(S(1,4))
    po= datenum([year1,month1,day1,hour1,min1,S(1,4)]);
    %disp(po)
    if Nmic == 1
        
        mksac(strcat(filename,'.OSU1.STW.TW.sac'),S(:,1),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','W', 'KHOLE', 'STW', 'KCMPNM', 'TW', 'STLA',  36.1344,  'STLO', -97.0819 )
    end
    if Nmic ==2
        mksac(strcat(filename,'.OSU1.W.sac'),S(:,1),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','W', 'KHOLE', 'STW', 'KCMPNM', 'TW', 'STLA',  36.1344,  'STLO', -97.0819 )
        mksac(strcat(filename,'.OSU1.E.sac'),S(:,2),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','E', 'KHOLE', 'STW', 'KCMPNM', 'TW', 'STLA',  36.1342,  'STLO', -97.0813 )
    end
    if Nmic ==3
        
        mksac(strcat(filename,'.OSU1.W.sac'),S(:,1),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','W', 'KHOLE', 'STW', 'KCMPNM', 'OK', 'STLA',  36.1344,  'STLO', -97.0819 )
        mksac(strcat(filename,'.OSU1.E.sac'),S(:,2),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','E', 'KHOLE', 'STW', 'KCMPNM', 'OK', 'STLA',  36.1342,  'STLO', -97.0813 )
        mksac(strcat(filename,'.OSU1.N.sac'),S(:,3),po,'DELTA',1/1000, 'KNETWK', 'OSU1','KSTNM','N', 'KHOLE', 'STW', 'KCMPNM', 'OK', 'STLA',  36.1347,  'STLO', -97.0814 )
    end

    %% clean up workspace
    % hold off

    tend=toc(tstart);
    
    disp(strcat('datagetter.m runtime = ',' ', num2str(tend), ' seconds...'))
    micnum=Nmic;
end






function [ordered] = Chronological(files)
ordered = files;
times = NaN(length(files),1);
for i = 1:length(files)
    filename = files(i).name;
    if strcmp(filename(1:4),'logD') ~= 1
        time = FileNameTime(filename);
        times(i) = datenum(time);
    end
end

times = times - times(1);
%times(isnan(times))=1;
Tsort = sortrows(times);        

for i = 1:length(Tsort)
    for j = 1:length(times)
        if Tsort(i) == times(j)
            ordered(i) = files(j);
        end
    end
end


end


function [] = DateFolderCk(year,month,day,output_location)
%check if year folder exists
YEAR = num2str(year); 
dirX = dir(strcat(output_location,YEAR));
if isempty(dirX) == 1
    mkdir(strcat(output_location,YEAR));
end

%check if month folder exists
if month < 10
    MON = strcat('0',num2str(month));
else
    MON = num2str(month);
end
dirX = dir(strcat(output_location,YEAR,'\',MON));
if isempty(dirX) == 1
    mkdir(strcat(output_location,YEAR,'\',MON));
end

%check if day folder exists
if day < 10
    DAY = strcat('0',num2str(day));
else
    DAY = num2str(day);
end
dirX = dir(strcat(output_location,YEAR,'\',MON,'\',DAY));
if isempty(dirX) == 1
    mkdir(strcat(output_location,YEAR,'\',MON,'\',DAY));
end  
end


function [DATA,N,fs] = InfraData(START,END,EVENT,location)

ct = 0; ck = 0;
YEARS = dir(location);
for yr = 3:length(YEARS)
    Year = str2double(YEARS(yr).name(1:4));
    if YEARS(yr).isdir == 1 && Year >= START(1) && Year <= END(1)
        folder = strcat(YEARS(yr).name,'/');
        MONTHS = dir(strcat(location,folder));
        for mon = 3:length(MONTHS)
            Month = str2double(MONTHS(mon).name(1:2));
            % identify max and min date range for given folder
            Dmin = datenum(Year,Month,1,0,0,0);
            if mon == 12
                Dmax = datenum(Year+1,1,1,0,0,0);
            else
                Dmax = datenum(Year,Month+1,1,0,0,0);
            end
            if MONTHS(mon).isdir == 1 && Dmax>datenum(START) && Dmin <= datenum(END)
                folder = strcat(YEARS(yr).name,'/',MONTHS(mon).name,'/');
                DAYS = dir(strcat(location,folder));
                for dy = 3:length(DAYS)
                    Day = str2double(DAYS(dy).name(1:2));
                    Dmin = datenum(Year,Month,Day,0,0,0);
                    Dmax = datenum(Year,Month,Day,23,59,59);
                    if DAYS(dy).isdir == 1 && Dmax>=datenum(START) && Dmin <= datenum(END)
                        folder = strcat(YEARS(yr).name,'/',MONTHS(mon).name,'/',DAYS(dy).name,'/');
                        files = dir(strcat(location,folder,'*.mat'));
                        files = Chronological(files);
                        logD = load(strcat(location,folder,'logD.mat')); logD = logD.logD;
                        
                        for i = 1:length(files)
                            filename = files(i).name;
                            if strcmp(filename(1:4),'logD') ~= 1
                                varT = FileNameTime(filename);
                                Y = varT(1); M = varT(2); D = varT(3);
                                h = varT(4); m = varT(5); s = varT(6);
                                
                                %% extract information from log file for give run
                                for j = 1:size(logD,1)
                                    if isequal(varT,logD(j,1:6)) == 1
%                                         DAYLIGHT = logD(j,7);           %1=CST; 2=CDT
%                                         DAQ = logD(j,8:10);              %1=yes,0=no; [SoundPressure,Accelerometer,Geophone]
                                        Period = logD(j,11);            %sample period in seconds
                                        fs = logD(j,12);                %sampling frequency (Hz)
%                                         Tseg = logD(j,13);              %length of each segment that gets saved to the computer
                                    end
                                end
                                
                                %% Establish times relative to event
                                T_int = datenum(varT);
                                T_final = datenum(Y,M,D,h,m,s+Period);
                                
                                if T_int <= datenum(END) && T_final >= datenum(START)
                                    Ts = (T_int-datenum(EVENT))*24;            %start time in hours relative to event (start of the day)
                                    
                                    %load data from the given file
                                    varS = load(strcat(location,folder,filename));
                                    S = varS.SP.rData(:,:);
                                    [numrows,numcols]=size(S);
%                                     if numcols > 1
%                                         disp('help')
%                                     end
                                    time = (0:1/fs:size(S,1)/fs-1/fs)';     %time in seconds
                                    time = time/(60*60) + Ts;               %time in hours relative to start of file
                                    
                                    %% determine whether to use entire segment or part of it [aa:bb]
                                    Range = [0,1];          %percent of total sample to include
                                    if T_int < datenum(START)
                                        Range(1) = (datenum(START)-T_int)/(T_final-T_int);
                                    end
                                    if T_final > datenum(END)
                                        Range(2) = 1 - (T_final-datenum(END))/(T_final-T_int);
                                    end
                                    
                                    aa=floor(Range(1)*(size(S,1)-1)/(1-0)+1);
                                    bb=ceil(Range(2)*(size(S,1)-1)/(1-0)+1);
                                    aa = max(1,aa); bb = min(size(S,1),bb);
                                    
                                    %% combine segments into a single DATA file
                                    ct = ct + 1; ck = 1;
                                    tm = time(aa:bb)*60;     %relative time (min)
                                    p = S(aa:bb,:);      %pressure (Pa)
                                    
                                    if ct == 1
                                        data = [tm,p];
                                    else
                                        Ng = round((tm(1)-data(end,1))*60*fs)-1;
                                        gap(:,1) = data(end,1)+1/(fs*60):1/(fs*60):Ng/(fs*60)+data(end,1);
                                        
                                        for j = 1:size(S,2)
                                            gap(:,j+1) = NaN(length(gap(:,1)),1);
                                        end
                                        if round(1/((gap(1,1)-data(end,1))*60))-fs > 0.2*fs || round(1/((tm(1)-gap(end,1))*60))-fs > 0.1*fs
                                            display(['deviation of ' num2str(max(10^-2*round(100*fs^-1/((gap(1,1)-data(end,1))*60))-1,10^-2*round(100*fs^-1/((tm(1)-gap(end,1))*60))-1))])
                                        end
                                        
                                        Dnew = [tm,p];
                                        data = [data;gap;Dnew];
                                        clear gap Dnew varS S time tm p
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

N = size(data,2)-1; 

if ck == 0
    DATA = NaN;
    disp(['No data available from ',datestr(datenum(START)),' to ',datestr(datenum(END))])  
else
    for i = 1:N
        MIC = strcat('mic',num2str(i));
        DATA.(MIC) = [data(:,1) data(:,i+1)];
    end
end
end



function [out1] = IsWdwAvailable(Time,wdw,location)
out1 = 0;
YEAR = dir(location);
for yr = 3:length(YEAR)
    Year = str2double(YEAR(yr).name(1:4));
    if YEAR(yr).isdir == 1 && Time(1) == Year
        folder = strcat(YEAR(yr).name,'/');
        MONTHS = dir(strcat(location,folder));
        
        for mon = 3:length(MONTHS)
            Month = str2double(MONTHS(mon).name(1:2));
            if MONTHS(mon).isdir == 1 && Time(2) == Month
                folder = strcat(YEAR(yr).name,'/',MONTHS(mon).name,'/');
                DAYS = dir(strcat(location,folder));
                
                for dy = 3:length(DAYS)
                    Day = str2double(DAYS(dy).name(1:2));
                    if DAYS(dy).isdir == 1 && Time(3) == Day
                        folder = strcat(YEAR(yr).name,'/',MONTHS(mon).name,'/',DAYS(dy).name,'/');
                        FILES = dir(strcat(location,folder,'*.mat'));
                        %disp(FILES)
                        FILES = Chronological(FILES);
                        logD = load(strcat(location,folder,'logD.mat')); logD = logD.logD;

                        for file = 1:length(FILES)
                            filename = FILES(file).name;
                            if strcmp(filename(1:4),'logD') ~= 1
                                varT = FileNameTime(filename);
                                Y = varT(1); M = varT(2); D = varT(3);
                                h = varT(4); m = varT(5); s = varT(6);
                                
                                % extract information from log file for give run
                                for j = 1:size(logD,1)
                                    if isequal(varT,logD(j,1:6)) == 1
                                        DAYLIGHT = logD(j,7);           %1=CST; 2=CDT
                                        DAQ = logD(j,8:10);              %1=yes,0=no; [SoundPressure,Accelerometer,Geophone]
                                        Period = logD(j,11);            %sample period in seconds
                                        fs = logD(j,12);                %sampling frequency (Hz)
                                        Tseg = logD(j,13);              %length of each segment that gets saved to the computer
                                    end
                                end
                                
                                %interval recorded in the given file
                                T_int = datenum(varT);
                                T_final = datenum(Y,M,D,h,m,s+Period);
                                
                                event = Time;
                                for seg = 0:10
                                    hr = Time(4) + seg*wdw/10;
                                    event(4) = hr;
                                    T_event = datenum(event);
                                    if T_event >= T_int && T_event <= T_final
                                        out1 = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

end



function mksac(f,d,varargin)
%MKSAC Write data in SAC format file.
%	MKSAC(FILENAME,D) writes a file FILENAME from data vector D in the 
%	Seismic Analysis Code (SAC) format. But additional options are highly
%	recommended to produce a real consistent seismic file:
%
%	MKSAC(FILENAME,D,T0,H) will use origin time T0 (DATENUM format), and
%	header fields in structure H (as defined by the IRIS/SAC format).
%
%	MKSAC(FILENAME,D,T0,'HEADER1',header1,'HEADER2',header2, ...) is an 
%	alternative to define the header fields.
%
%	MKSAC(FILENAME,D,T,...) where T is a time vector of same size as D, in
%	datenum format, will define DELTA sampling value from the time interval
%	median but will be overwritten by header field DELTA if specified.
%
%	MKSAC will produce a SAC file in any case, using default values for 
%	any missing header field. But we strongly suggest to define at least 
%	the following header fields:
%	- DELTA (increment between evenly spaced samples, in seconds), default
%	  will be 1 second;
%	- KSTNM (station code), KHOLE (location code), KCMPNM (channel code),
%	  and KNETWK (network code) to define the station component;
%
%	Header field names that are not recognized are ignored, while some will
%	be overwritten in order to keep the file consistency:
%	- NPTS, DEPMIN, DEPMAX, DEPMEN are inferred from data D
%	- NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC are defined from time T0
%	  rounded to the closest millisecond.
%	- B is forced to 0, E to the last sample relative time (in seconds)
%	- NVHDR is forced to 6 (header version number)
%
%	Example:
%	   mksac('test.sac',sin(linspace(0,10*pi)),now,'DELTA',1/100,'KSTNM','TEST')
%
%
%	Reference: http://www.iris.edu/files/sac-manual/
%
%	Author: F. Beauducel <beauducel@ipgp.fr>
%	Created: 2015-11-12
%	Updated: 2020-11-11

%	Release history:
%	[2020-11-11] v1.3
%		- fix an issue in origin time NZMSEC (thanks to Randall Plate)
%	[2020-01-15] v1.2
%		- fix a missing DELTA setting from T vector (thanks to Randall
%		  Plate comment)
%	[2016-03-05] v1.1
%		- fix a problem with date of origin time
%	[2015-11-12] v1.0
%
%	Copyright (c) 2020, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the 
%	     distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
%	IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
%	TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%	PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
%	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
%	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
%	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
%	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
%	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin < 2
	error('Not enough input argument.')
end

if nargin < 3
	t0 = now;
else
	t0 = varargin{1};
end

if nargin > 3 && isstruct(varargin{2})
	H = varargin{2};
end

if nargin > 4
	for n = 2:2:(nargin-3)
		if ischar(varargin{n})
			try
				H.(varargin{n}) = varargin{n+1};
			catch
				fprintf('Warning: invalid header name "%s" or value.\n',varargin{n});
			end
		end
	end
end

% header sampling rate
if numel(t0) > 1 && (~exist('H','var') || ~isfield(H,'DELTA'))
	dt = 86400*median(diff(t0));
	if dt > 0
		H.DELTA = dt;
	end
end


% header default values
H0 = struct('DELTA',1,'NVHDR',6,'IFTYPE',1,'LEVEN',1);
for h = fieldnames(H0)'
	if ~exist('H','var') || ~isfield(H,h{:})
		H.(h{:}) = H0.(h{:});
	end
end

% header origin time values
msinday = 24*60*60*1000; % milliseconds in a day
tv = datevec(round(t0(1)*msinday)/msinday); % rounds to the closest ms
H.NZYEAR = tv(1);
H.NZJDAY = datenum(tv(1:3)) - datenum(tv(1),1,1) + 1;
H.NZHOUR = tv(4);
H.NZMIN = tv(5);
H.NZSEC = floor(tv(6));
H.NZMSEC = (tv(6) - H.NZSEC)*1e3;

% other header overwritten fields
H.NPTS = length(d);
H.DEPMIN = min(d);
H.DEPMAX = max(d);
H.DEPMEN = mean(d);
H.B = 0;
H.NVHDR = 6;

% splits KEVNM in 2 bytes
if isfield(H,'KEVNM') && ischar(H.KEVNM) && length(H.KEVNM) == 16
	H.KEVNM0 = H.KEVNM(1:8);
	H.KEVNM1 = H.KEVNM(9:16);
	H = rmfield(H,'KEVNM');
end

fid = fopen(f, 'wb', 'ieee-le');
if fid == -1
	error('Cannot write the data file %s',f);
end

writeheader(fid,H);
fwrite(fid,single(d),'float32');	% writes data as single class

fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeheader(fid,H)

novalue = -12345;

% --- classifies header fields
% floating variables (single)
v = { ...
'DELTA',   'DEPMIN',  'DEPMAX',  'SCALE',  'ODELTA';
'B',       'E',       'O',       'A',      'INTERNAL';
'T0',      'T1',      'T2',      'T3',     'T4';
'T5',      'T6',      'T7',      'T8',     'T9';
'F',       'RESP0',   'RESP1',   'RESP2',  'RESP3';
'RESP4',   'RESP5',   'RESP6',   'RESP7',  'RESP8';
'RESP9',   'STLA',    'STLO',    'STEL',   'STDP';
'EVLA',    'EVLO',    'EVEL',    'EVDP',   'MAG';
'USER0',   'USER1',   'USER2',   'USER3',  'USER4';
'USER5',   'USER6',   'USER7',   'USER8',  'USER9';
'DIST',    'AZ',      'BAZ',     'GCARC',  'INTERNAL';
'INTERNAL','DEPMEN',  'CMPAZ',   'CMPINC', 'XMINIMUM';
'XMAXIMUM','YMINIMUM','YMAXIMUM','UNUSED', 'UNUSED';
'UNUSED',  'UNUSED',  'UNUSED',  'UNUSED', 'UNUSED';
}';

hf = single(repmat(novalue,size(v)));
for n = 1:numel(v)
	if isfield(H,v{n})
		hf(n) = single(H.(v{n}));
	end
end

% integer variables (int32)
v = { ...
'NZYEAR',  'NZJDAY',  'NZHOUR',  'NZMIN',  'NZSEC';
'NZMSEC',  'NVHDR',   'NORID',   'NEVID',  'NPTS';
'INTERNAL','NWFID',   'NXSIZE',  'NYSIZE', 'UNUSED';
'IFTYPE',  'IDEP',    'IZTYPE',  'UNUSED', 'IINST';
'ISTREG',  'IEVREG',  'IEVTYP',  'IQUAL',  'ISYNTH';
'IMAGTYP', 'IMAGSRC', 'UNUSED',  'UNUSED', 'UNUSED';
'UNUSED',  'UNUSED',  'UNUSED',  'UNUSED', 'UNUSED';
}';

hi = int32(repmat(novalue,size(v)));
for n = 1:numel(v)
	if isfield(H,v{n})
		% case of enumerated fields I* that may contain 'Description {N}'
		if strncmp(v{n},'I',1) && ischar(H.(v{n}))
			d = regexp(H.(v{n}),' {(\d*)}','tokens');
			if ~isempty(d{1}{1})
				hi(n) = str2double(d{1}{1});
			end
		end
		if isnumeric(H.(v{n}))
			hi(n) = int32(H.(v{n}));
		end
	end
end

% logical variables (int32)
v = { ...
'LEVEN',   'LPSPOL',  'LOVROK',  'LCALDA', 'UNUSED';
}';

hl = int32(zeros(size(v)));
for n = 1:numel(v)
	if isfield(H,v{n})
		hl(n) = int32(H.(v{n}));
	end
end


% alphanumerical variables (char)
v = { ...
'KSTNM',  'KEVNM0', 'KEVNM1';
'KHOLE',  'KO',     'KA';
'KT0',    'KT1',    'KT2';
'KT3',    'KT4',    'KT5';
'KT6',    'KT7',    'KT8';
'KT9',    'KF',     'KUSER0';
'KUSER1', 'KUSER2', 'KCMPNM';
'KNETWK', 'KDATRD', 'KINST';
}';

hc = char(repmat(sprintf('%-8d',novalue),1,numel(v)));
for n = 1:numel(v)
	if isfield(H,v{n})
		s = char(H.(v{n}));
		if length(s) > 8
			s = s(1:8);
		end
		hc((n-1)*8 + (1:8)) = sprintf('%-8s',s);
	end
end


% writes the header
fwrite(fid,hf,'float32');
fwrite(fid,[hi,hl],'int32');
fwrite(fid,hc','char');


end



function [CT] = OKutc2ct(Date,Time)

DST = [2010,3,14,11,7; 2011,3,13,11,6;2012,3,11,11,4;2013,3,10,11,3;2014,3,9,11,2;2015,3,8,11,1; 2016,3,13,11,6;2017,3,12,11,5; 2018,3,11,11,4;2019,3,10,11,3; 2020,3,8,11,1;2021,3,14,11,7;2022,3,13,11,6];

yr = Date(1); mon = Date(2); dy = Date(3);
hr = Time(1); mn = Time(2); sec = Time(3);

if yr >= 2010 && yr <= 2022
    %locate the dates and times for time changes for the given year
    for i = 1:length(DST)
        if DST(i,1) == yr
            dst = [DST(i,1:3),2,0,0; DST(i,1),DST(i,4:5),2,0,0];
        end
    end
    
    if datenum([Date,hr-6,mn,sec]) >= datenum(dst(1,:)) && datenum([Date,hr-5,mn,sec]) < datenum(dst(2,:))
        ct = [yr,mon,dy,hr-5,mn,sec];
    else
        ct = [yr,mon,dy,hr-6,mn,sec];
    end
    
    
    CT = [ct(1), ct(2), ct(3), ct(4), 0, 0];
    if mon == dst(2,2) && dy == dst(2,3) 
        if hr == dst(2,4)-1+5 || hr == dst(2,4)-1+6 
            disp('The CT time is NOT unique')
        end
    end
else
    disp(['Outside of Range (',num2str(DST(1,1)),' to ',num2str(DST(end,1)),')'])
end
end

function [out1] = FileNameTime(filename)
Y = str2double(filename(5:8));
M = str2double(filename(1:2));
D = str2double(filename(3:4));
h = str2double(filename(10:11));
m = str2double(filename(12:13));
s = str2double(filename(14:15));

if strcmp('AM',filename(17:18)) == 1 || strcmp('PM',filename(17:18)) == 1
    if h == 12
        if strcmp('AM',filename(17:18)) == 1
            h = 0;
        else
            h = 12;
        end
    else
        if strcmp('PM',filename(17:18)) == 1
            h = h + 12;
        end
    end
end
out1 = [Y,M,D,h,m,s];
end