clc;
clear;
close all;

%% Select folder
folderPath = uigetdir(pwd, 'Select Zemax folder');

if isequal(folderPath, 0)
    return;
end

% Find ZMX files whose names start with "MC_T"
filePrefix = 'MC_T';
files = dir(fullfile(folderPath, [filePrefix '*.zmx']));
filesUpper = dir(fullfile(folderPath, [filePrefix '*.ZMX']));
files = [files; filesUpper];

% Windows dir() is case-insensitive; keep only names starting exactly with the prefix
files = files(startsWith(string({files.name}), filePrefix));

if isempty(files)
    error('No ZMX files starting with "%s" found.', filePrefix);
end

% Remove duplicate filenames, if any
fileNamesAll = string({files.name});
fileNamesLower = lower(fileNamesAll);
[~, uniqueIdx] = unique(fileNamesLower, 'stable');
files = files(uniqueIdx);

%% Initialize ZOS-API connection
% InitConnection and CleanupConnection are from the
% official Zemax MATLAB standalone sample files.

TheApplication = InitConnection();

if isempty(TheApplication)
    error('Could not connect to OpticStudio.');
end

%% Read EFL values
fileNames = strings(numel(files), 1);
EFL_mm = nan(numel(files), 1);
status = strings(numel(files), 1);

fprintf('Processing %d files...\n', numel(files));

try
    TheSystem = TheApplication.PrimarySystem;

    for k = 1:numel(files)

        filePath = fullfile(folderPath, files(k).name);
        fileNames(k) = string(files(k).name);

        try
            % Load optical system
            loaded = TheSystem.LoadFile(filePath, false);

            if ~loaded
                error('OpticStudio could not load this file.');
            end

            % Calculate effective focal length
            EFL_mm(k) = TheSystem.MFE.GetOperandValue( ...
                ZOSAPI.Editors.MFE.MeritOperandType.EFFL, ...
                0, 0, 0, 0, 0, 0, 0, 0);

            if isfinite(EFL_mm(k))
                status(k) = "OK";
            else
                status(k) = "Invalid EFL";
            end

        catch ME
            status(k) = string(ME.message);
        end
    end

catch ME
    CleanupConnection(TheApplication);
    rethrow(ME);
end

CleanupConnection(TheApplication);

%% Results table (file -> EFL)
Results = table(fileNames, EFL_mm, status, ...
    'VariableNames', {'FileName', 'EFL_mm', 'Status'});

%% Keep valid results only
valid = isfinite(EFL_mm) & status == "OK";
eflValid = EFL_mm(valid);

nFailed = sum(~valid);
if nFailed > 0
    fprintf('%d file(s) could not be evaluated (see Status column in Excel).\n', nFailed);
end

if isempty(eflValid)
    warning('No valid EFL values were obtained.');
    return;
end

%% Statistical summary
statNames  = ["Count"; "Mean"; "Median"; "Std"; "Min"; "Max"];
statValues = [numel(eflValid); mean(eflValid); median(eflValid); ...
              std(eflValid); min(eflValid); max(eflValid)];

fprintf('\nNumber of valid files: %d\n', numel(eflValid));
fprintf('Mean EFL:   %.4f mm\n', statValues(2));
fprintf('Median EFL: %.4f mm\n', statValues(3));
fprintf('Std EFL:    %.4f mm\n', statValues(4));
fprintf('Min EFL:    %.4f mm\n', statValues(5));
fprintf('Max EFL:    %.4f mm\n', statValues(6));

%% Histogram: EFL distribution
fig = figure('Name', 'EFL Distribution', 'Color', 'w');

histogram(eflValid, 'BinMethod', 'auto');

xlabel('Effective Focal Length (mm)');
ylabel('Number of optical systems');
title(sprintf('EFL Distribution - %d files', numel(eflValid)));
grid on;

% Save histogram as image to the selected folder
imageFile = fullfile(folderPath, 'EFL_Histogram.png');
try
    exportgraphics(fig, imageFile, 'Resolution', 300);   % R2020a+
catch
    saveas(fig, imageFile);
end

%% Export Excel
outputFile = fullfile(folderPath, 'EFL_Results.xlsx');

% Remove old file so stale rows/sheets are not kept
if isfile(outputFile)
    delete(outputFile);
end

writetable(Results, outputFile, 'Sheet', 'EFL by File');

Stats = table(statNames, statValues, ...
    'VariableNames', {'Statistic', 'EFL_mm'});
writetable(Stats, outputFile, 'Sheet', 'Statistics');

fprintf('\nHistogram saved to:\n%s\n', imageFile);
fprintf('Results saved to:\n%s\n', outputFile);


function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = 'C:\Users\gayeisik\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll';

% This is the path to OpticStudio
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\Ansys Zemax OpticStudio 2025 R2.04');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    LogMessage('NOT-Found OpticStudio');
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
TheConnection.ConnectionTimeoutSeconds = 20;
app = TheConnection.CreateNewApplication();

if isempty(app)
    HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MException('zosapi:HandleError', error);
throw(ME);
end

function CleanupConnection(TheApplication)
% Note - this will close down the connection.
TheApplication.CloseApplication();
end
