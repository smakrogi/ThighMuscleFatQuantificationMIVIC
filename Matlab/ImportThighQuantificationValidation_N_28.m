%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\smakrogiannis\Documents\MATLAB\MRIThighAnalysis\MetS_Experiments\ThighQuantificationValidation_NS_only_FCM_20100614_left_Matching_MetS_BMI.xlsx
%    Worksheet: Geanie&3T
%
% Auto-generated by MATLAB on 13-Mar-2020 11:31:15

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 53);

% Specify sheet and range
opts.Sheet = "Geanie&3T";
opts.DataRange = "B3:BB18";

% Specify column names and types
opts.VariableNames = ["VarName2", "SubjectID", "Visit", "Age", "Gender", "AcqDate", "Comments", "Leg", "CTTHCSAM", "CTTHDenM", "CTTHCSAF", "CTTHCSAO", "CTTHCSAA", "CTTHCSA", "VarName16", "VarName17", "Leg1", "FirstUsedSlice", "LastUsedSlice", "SATvolumemm3", "AverageSATareamm", "Musclevolumemm3", "AverageMuscleareamm", "IMFATvolumemm3", "AverageIMFATareamm", "TotalFatareamm", "TotalVolume", "MuscleFat", "MuscleIMAT", "MuscleTotalVolume", "M_Volheight2", "M_Volweight", "M_VolBMI", "VarName35", "idno", "visit", "LATEST_DOV", "obesity", "Hyperglycemia", "Dyslipidemia_tgl", "Dyslipidemia_hdl", "Hypertension", "metS_score", "B_M_I", "weight", "height", "bmi", "WtKg", "HtCm", "waist_circ", "hip_circ", "WHR", "waistdep"];
opts.VariableTypes = ["categorical", "double", "double", "double", "categorical", "datetime", "string", "categorical", "double", "double", "double", "double", "double", "double", "string", "string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "double", "double", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Comments", "VarName16", "VarName17", "VarName35"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName2", "Gender", "Comments", "Leg", "VarName16", "VarName17", "Leg1", "VarName35"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "AcqDate", "InputFormat", "");
opts = setvaropts(opts, "LATEST_DOV", "InputFormat", "");

% Import the data
DataMatrix = readtable("C:\Users\smakrogiannis\Documents\MATLAB\MRIThighAnalysis\MetS_Experiments\ThighQuantificationValidation_NS_only_FCM_20100614_left_Matching_MetS_BMI.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts