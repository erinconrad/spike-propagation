function portalTimings


%% Parameters to change each time
day0 = '01/02/00';
time0 = '22:06:08';
daysz = '1/05/00';
timesz = '1:41 AM';

%% Figure out how many seconds into portal viewer seizure is
sz_date = datetime([daysz,' ',timesz]);
start_date = datetime([day0, ' ',time0]);

between_time = seconds(sz_date - start_date);

fprintf('There are %d seconds between\n%s and %s\n',...
    between_time,datestr(start_date),datestr(sz_date));


end