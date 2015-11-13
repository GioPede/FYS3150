x = [1/20, 1/40, 1/60, 1/80, 1/100, 1/150, 1/200]

filename = '20.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T20 = dataArray{:, 1};
E20 = dataArray{:, 2};
C20 = dataArray{:, 3};
M20 = dataArray{:, 4};
S20 = dataArray{:, 5};

filename = '40.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T40 = dataArray{:, 1};
E40 = dataArray{:, 2};
C40 = dataArray{:, 3};
M40 = dataArray{:, 4};
S40 = dataArray{:, 5};


filename = '60.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T60 = dataArray{:, 1};
E60 = dataArray{:, 2};
C60 = dataArray{:, 3};
M60 = dataArray{:, 4};
S60 = dataArray{:, 5};


filename = '80.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T80 = dataArray{:, 1};
E80 = dataArray{:, 2};
C80 = dataArray{:, 3};
M80 = dataArray{:, 4};
S80 = dataArray{:, 5};


filename = '100.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T100 = dataArray{:, 1};
E100 = dataArray{:, 2};
C100 = dataArray{:, 3};
M100 = dataArray{:, 4};
S100 = dataArray{:, 5};

filename = '150.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T150 = dataArray{:, 1};
E150 = dataArray{:, 2};
C150 = dataArray{:, 3};
M150 = dataArray{:, 4};
S150 = dataArray{:, 5};


filename = '200.txt';
formatSpec = '%5f%12f%10f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
T200 = dataArray{:, 1};
E200 = dataArray{:, 2};
C200 = dataArray{:, 3};
M200 = dataArray{:, 4};
S200 = dataArray{:, 5};



T_max20C = find(C20 == max(C20));
T_max20S = find(S20 == max(S20));

T_max40C = find(C40 == max(C40));
T_max40S = find(S40 == max(S40));

T_max60C = find(C60 == max(C60));
T_max60S = find(S60 == max(S60));

T_max80C = find(C80 == max(C80));
T_max80S = find(S80 == max(S80));

T_max100C = find(C100 == max(C100));
T_max100S = find(S100 == max(S100));

T_max150C = find(C150 == max(C150));
T_max150S = find(S150 == max(S150));

T_max200C = find(C200 == max(C200));
T_max200S = find(S200 == max(S200));


y = [(T20(T_max20C)   + T20(T_max20C))  /2, 
     (T40(T_max40C)   + T40(T_max40C))  /2, 
     (T60(T_max60C)   + T60(T_max60C))  /2, 
     (T80(T_max80C)   + T80(T_max80C))  /2, 
     (T100(T_max100C) + T100(T_max100C))/2,
     (T150(T_max150C) + T150(T_max150C))/2,
     (T200(T_max200C) + T200(T_max200C))/2]




