classdef ReliaCompGUI < matlab.apps.AppBase
%***********************************************************************
% Copyright (C) 2020  Mayank Chetan
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% ***********************************************************************
% This runs a Monte Carlo Simulations (MCS) approximation to find the
% reliability related to the MECH6338 project. The function uses the
% Cholesky decomposition method[1]. The functional approximation MAT file
% obtained by running the "FitTANA.m" script is used here.
%
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************
    % Properties that correspond to app components
    properties (Access = public)
        ReliaComp                       matlab.ui.Figure
        Menu                            matlab.ui.container.Menu
        ResetFieldsMenu                 matlab.ui.container.Menu
        LoadMenu                        matlab.ui.container.Menu
        SaveMenu                        matlab.ui.container.Menu
        ExitMenu                        matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        DocumentationMenu               matlab.ui.container.Menu
        AboutMenu                       matlab.ui.container.Menu
        ControlPanel                    matlab.ui.container.Panel
        NumberofVariablesSpinnerLabel   matlab.ui.control.Label
        NoOfVarsSpinner                 matlab.ui.control.Spinner
        FuncInputField                  matlab.ui.control.EditField
        SelectFunctionButton            matlab.ui.control.Button
        ComputationMethodDropDownLabel  matlab.ui.control.Label
        ComputationMethodDropDown       matlab.ui.control.DropDown
        RUNButton                       matlab.ui.control.Button
        MethodSpecificInputsPanel       matlab.ui.container.Panel
        TabGroup2                       matlab.ui.container.TabGroup
        HLRFTab                         matlab.ui.container.Tab
        GradientMethodDropDownLabel     matlab.ui.control.Label
        GradDropDown                    matlab.ui.control.DropDown
        ConvergenceAlphaEditFieldLabel  matlab.ui.control.Label
        AlphaEditField                  matlab.ui.control.NumericEditField
        ConvergenceBetaEditFieldLabel   matlab.ui.control.Label
        BetaEditField                   matlab.ui.control.NumericEditField
        GradientStepLabel               matlab.ui.control.Label
        GradPertubEditField             matlab.ui.control.NumericEditField
        MonteCarloTab                   matlab.ui.container.Tab
        NumberofsamplesforMCEditFieldLabel  matlab.ui.control.Label
        samplesforMCEditField           matlab.ui.control.NumericEditField
        MVFOSMTab                       matlab.ui.container.Tab
        WarningThismethodassumesalltheLabel  matlab.ui.control.Label
        InteractionPanel                matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        InputsTab                       matlab.ui.container.Tab
        InputTable                      matlab.ui.control.Table
        CorrTable                       matlab.ui.control.Table
        CorrelationMatrixLabel          matlab.ui.control.Label
        OutputTab                       matlab.ui.container.Tab
        OutputTextArea                  matlab.ui.control.TextArea
        ResultsTab                      matlab.ui.container.Tab
        ResultsTable                    matlab.ui.control.Table
    end

    
    properties (Access = public)
        DataStruct = struct('input',[],'model',[]); % Description
    end
    
    methods (Access = public)
        
        function UpdateOutputWindow(app,nextLine,type)
            if type == 1
                DispString = sprintf('[%s] \t %s',datestr(now,'ddmmmyyyy HH:MM:SS'),nextLine);
                app.OutputTextArea.Value = [DispString;app.OutputTextArea.Value];
                pause(0.01)
            else
                app.OutputTextArea.Value{1} = sprintf('[%s] \t %s',datestr(now,'ddmmmyyyy HH:MM:SS'),nextLine);
            end
        end
        
        function Switch2InputTab(app)
            app.TabGroup.SelectedTab = app.InputsTab;
        end
        
    end
    
    methods (Access = private)
        
        function results = VerifyInputData(app)
        
            results = 0;
            
            % Getting the information from all the inputs and storing it in
            % the public property structure!
            
            UpdateOutputWindow(app,'Checking input function....',1)
            
            if isempty(app.DataStruct.input)
                UpdateOutputWindow(app,'No input function Found! please select input function',1)
                % Alerting user and switching over to the Input tab!
                uialert(app.ReliaComp,'Input function not found, please select input function','Error')
                return
            else
                UpdateOutputWindow(app,'Checking input function.... Done',0)
            end
            
            UpdateOutputWindow(app,'Adding path......',1)
            addpath(app.DataStruct.input.Path)
            UpdateOutputWindow(app,'Adding path...... Done',0)
            model.funcName = app.DataStruct.input.Func(1:end-2);
            
            
            UpdateOutputWindow(app,'Checking input data....',1)
            if app.NoOfVarsSpinner.Value == 0
                UpdateOutputWindow(app,'No input data provided!',1)
                uialert(app.ReliaComp,'Input variables not defined, please define input variables','Error')
                return
            else
                UpdateOutputWindow(app,sprintf('Checking input data.... %d variables defined!',app.DataStruct.input.NoOfVars),0)
            end
            
            
            % Reading in the specs of the iterations!
            model.gradFlag = strcmp(app.GradDropDown,'Central Diff');
            model.gradDelta = app.GradPertubEditField.Value; % numerical gradient perturbation in %
            model.BetaDiff = app.BetaEditField.Value; % Convergence criteria for Beta
            model.AlphaDiff = app.AlphaEditField.Value; % Convergence criteria for Alpha
            model.MCSamples = app.samplesforMCEditField.Value; % Number of samples for MC method.
            
            
            
            % Read the data and populate the model structure
            UpdateOutputWindow(app,'Reading in tabular data....',1)
            
            if any(isnan(table2array(app.InputTable.Data(:,[2,3]))),'all') ||...
                    any(isnan(app.CorrTable.Data),'all')
                UpdateOutputWindow(app,'NaN found in the tabular input data!',1)
                % Alerting user and switching over to the Input tab!
                uialert(app.ReliaComp,'NaN found in tabular input data please correct and re-run','Error')
                return
            elseif ~ValidateCorrMatrix(app)
                return
            else
                
                model.varsMeans = app.InputTable.Data.Mean;
                model.varsCOV = app.InputTable.Data.COV;
                
                model.varSD = model.varsMeans.*model.varsCOV; % Calc of SD's
                
                DistMap = {'Normal','Log Normal','Uniform'};
                for idx = 1:app.NoOfVarsSpinner.Value
                    model.Types(idx) = find(strcmp(string(app.InputTable.Data.Dist(idx)),DistMap)) - 1;
                end
                
                model.Cprime = app.CorrTable.Data;
                
                model.LNSD = sqrt(log(1 + model.varsCOV.^2));
                model.LNmean = log(model.varsMeans) - 0.5 .* model.LNSD .^2;
                
                UpdateOutputWindow(app,'Reading in tabular data.... Done!',0)
            end
            
            % Test the input function for the right number of input
            % variables!
            UpdateOutputWindow(app,'Testing input function for validity....',1)
            
            try
                feval(model.funcName,model.varsMeans);
                [A,B] = FindGrads(model.funcName,model.varsMeans,model.gradDelta);
                UpdateOutputWindow(app,'Testing input function for validity.... Done!',0)
            catch ME
                UpdateOutputWindow(app,'Testing input function for validity.... Failed!',0)
                switch ME.identifier
                    
                    case 'MATLAB:UndefinedFunction'
                        uialert(app.ReliaComp,'The input function is undefined!','Error')
                        return
                        
                    case 'MATLAB:badsubscript'
                        uialert(app.ReliaComp,'Missmatch in number of variable provided and number of inputs taken in by function.','Error')
                        return
                        
                    otherwise
                        uialert(app.ReliaComp,ME.message,'MATLAB Error')
                        return
                end
            end
            
            app.DataStruct.model = model;
            app.DataStruct.input.VarInTable = app.InputTable.Data;
            app.DataStruct.input.VarCorrTable = app.CorrTable.Data;
            results = 1;
            
        end
        
        function results = ValidateCorrMatrix(app)
            
            results = 0;
            
            if ~issymmetric(app.CorrTable.Data)
                
                UpdateOutputWindow(app,'Reading in tabular data.... Correlation Matrix is not Symmertric',0)
                uialert(app.ReliaComp,'Correlation Matrix is not Symmertric, Please correct it','Error')
                return
                
            elseif any(app.CorrTable.Data > 1)
                
                UpdateOutputWindow(app,'Reading in tabular data.... Correlation Matrix has element with value greater than 1',0)
                uialert(app.ReliaComp,'Correlation Matrix has element with value greater than 1, Please correct it','Error')
                return
                
            else
                
                [~,D] = eig(app.CorrTable.Data);
                
                D = real(diag(D));
                
                if ~all(D >= 0)
                    
                    UpdateOutputWindow(app,'Reading in tabular data.... Correlation Matrix is not positive semi-definite',0)
                    uialert(app.ReliaComp,'Correlation Matrix is not positive semi-definite, Please correct it','Error')
                    return
                else
                    
                    results = 1;
                    
                end
                
            end
            
            
        end
        
        function results = VerifyHLRFInput(app)
            results = 0;
            
            if app.DataStruct.model.gradDelta > 1e-3
                UpdateOutputWindow(app,'Please set grad delta to below 1e-3 to compute accurate gradients',1)  
                return
            end
            
            if app.DataStruct.model.BetaDiff > 1e-5
                UpdateOutputWindow(app,'Please set beta to below 1e-5 to have better convergence',1)    
                return
            end
            
            if app.DataStruct.model.AlphaDiff > 1e-4
                UpdateOutputWindow(app,'Please set alpha to below 1e-5 to have better convergence',1)       
                return
            end
            
            results = 1;
            
        end
        
        function results = VerifyMCInput(app)
            results = 0;
            if app.DataStruct.model.MCSamples > 1e9
                UpdateOutputWindow(app,'Please set lower number of samples for MC',1)
                return
            end
            results = 1;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.InputTable.RowName = 'numbered';
            app.CorrTable.RowName = 'numbered';
            clc
            addpath('utils') % Need to address this for packaging!
            
            %% Read in the standard values that are populated!
            
            app.DataStruct.input.CompMethod         = app.ComputationMethodDropDown.Value;
            app.DataStruct.input.GradMethod         = app.GradDropDown.Value;
            app.DataStruct.input.GradPeturb         = app.GradPertubEditField.Value ;
            app.DataStruct.input.Alpha              = app.AlphaEditField.Value;
            app.DataStruct.input.Beta               = app.BetaEditField.Value;
            app.DataStruct.input.MCsamples          = app.samplesforMCEditField.Value;
            
        end

        % Button pushed function: SelectFunctionButton
        function SelectFunctionButtonPushed(app, event)
            [file,path] = uigetfile('*.m');
            if isequal(file,0)
                disp('User selected Cancel');
            else
                app.FuncInputField.Value = file;
                app.DataStruct.input.Func = file;
                app.DataStruct.input.Path = path;
            end
            
            app.Switch2InputTab
        end

        % Value changed function: NoOfVarsSpinner
        function NoOfVarsSpinnerValueChanged(app, event)
            value = app.NoOfVarsSpinner.Value;
            
            cnames = categorical({'Normal'},{'Normal','Log Normal','Uniform'},'Ordinal',true);
            TableTemplate = table(" ",0,0,cnames,'VariableNames',{'VarName'; 'Mean'; 'COV'; 'Dist'});
            
            if value > 0
                
                app.InputTable.ColumnEditable = true;
                app.CorrTable.ColumnEditable = true;
                
                covTableVals = app.CorrTable.Data;
                
                if size(app.InputTable.Data,1) < value % The number is being increased
                    app.InputTable.Data = [app.InputTable.Data; TableTemplate]; % Setting up for Data input table
                    
                    if value == 1
                        app.CorrTable.Data = 1;
                        app.CorrTable.ColumnName = 1;
                    else
                        app.CorrTable.Data = [covTableVals zeros(size(covTableVals,1),1);zeros(1,value-1) 1];
                        app.CorrTable.ColumnName = 1:value;
                    end
                    
                    app.CorrTable.ColumnWidth = num2cell((app.CorrTable.InnerPosition(3) ./ (value+1)).* ones(1,value));
                    
                else % The number is being decreased or same
                    app.InputTable.Data = app.InputTable.Data(1:value,:);
                    app.CorrTable.Data = app.CorrTable.Data(1:value,1:value);
                end
                
            else % Has a value of Zero!
                app.InputTable.Data =[];
                app.CorrTable.Data =[];
                
            end
            
            app.DataStruct.input.NoOfVars = value;
        end

        % Button pushed function: RUNButton
        function RUNButtonPushed(app, event)
            % Here is where the fun starts!!!!!
            UpdateOutputWindow(app,'---------------------------------------',0)
            
            % Switching over to the output tab!
            app.TabGroup.SelectedTab = app.OutputTab;
            
            % Verify the input Data!
            if VerifyInputData(app)
                UpdateOutputWindow(app,'Input Checks Completed!',1)
            else
                UpdateOutputWindow(app,'Checks Failed! PLease follow error prompts to correct issues',1)
                return
            end
            
            switch(app.ComputationMethodDropDown.Value)
                
                case 'HL-RF'
                    
                    % Verifying input information related to HLRF!
                    UpdateOutputWindow(app,'Verifying input related to HLRF.......',1)
                    if ~VerifyHLRFInput(app)
                        UpdateOutputWindow(app,'Verifying input related to HLRF.......FAILED!',0)
                        uialert(app.ReliaComp,'HLRF input verification has failed, please check output','Error')
                        return
                    else
                    UpdateOutputWindow(app,'Verifying input related to HLRF.......PASSED!',0)
                    end
                    
                    
                    % Run the HL-RF Reliability analysis !!!!! need errror
                    % handling here !!!
                    UpdateOutputWindow(app,'Begining Computation......',1)
                    app.DataStruct.HLRF.rel = findReliabilityHLRF(app);
                    
                    % Changing the tab to the Output......!!!!
                    app.TabGroup.SelectedTab = app.ResultsTab;
                    tempTable = struct2table(app.DataStruct.HLRF.rel);
                    ResTable = table(tempTable.beta,tempTable.R .* 100,tempTable.Pf .* 100);
                    app.ResultsTable.Data = ResTable;
                    app.ResultsTable.ColumnName = {'Beta (-)'; 'Reliability (%)'; 'Probability of Failure (%)'};
                    
                    % Can be done in a better way, this is a temp fix!
                    app.ResultsTable.ColumnName(end+1:end+app.NoOfVarsSpinner.Value) =...
                        table2cell(app.InputTable.Data(:,1));
                    tabLen = length(app.ResultsTable.ColumnName);
                    for jdx = 1:length(app.DataStruct.HLRF.rel)
                        app.ResultsTable.Data(jdx,tabLen-app.NoOfVarsSpinner.Value+1:tabLen) = num2cell(tempTable.xStar{jdx})';
                    end
                    
                    app.ResultsTable.RowStriping = true;
                    
                    app.ResultsTable.BackgroundColor = repmat([1 1 1],[size(app.ResultsTable.Data,1),1]);
                    app.ResultsTable.BackgroundColor(end,:) = [0.5 1 0.5];
                    
                    assignin('base','Results',app.DataStruct);
                    
                case 'Monte Carlo'
                    
                    % Verifying input information related to Monte Carlo!
                    UpdateOutputWindow(app,'Verifying input related to Monte Carlo.......',1)
                    if ~VerifyMCInput(app)
                        UpdateOutputWindow(app,'Verifying input related to Monte Carlo.......FAILED!',0)
                        uialert(app.ReliaComp,'Monte Carlo input verification has failed, please check output','Error')
                        return
                    else
                        UpdateOutputWindow(app,'Verifying input related to Monte Carlo.......PASSED',0)
                    end
                   
                    
                    % Run the HL-RF Reliability analysis !!!!! need errror
                    % handling here !!!
                    UpdateOutputWindow(app,'Begining Computation......',1)
                    app.DataStruct.MCS.rel = findReliabilityMCS(app);
                    
                    disp = sprintf('The reliability for %2.2E MCS runs is %2.6f with 95%% CI (%2.6f,%2.6f)  \n',...
                        app.DataStruct.model.MCSamples,app.DataStruct.MCS.rel.R,app.DataStruct.MCS.rel.ConfDown,app.DataStruct.MCS.rel.ConfUp);
                    UpdateOutputWindow(app,disp,1)
                    
                    % Changing the tab to the Results......!!!!
                    app.TabGroup.SelectedTab = app.ResultsTab;
                    app.ResultsTable.ColumnName = {'Reliability (%)'; 'Upper 95% conf';'Lower 95% conf';'No Of Samples'};
                    app.ResultsTable.Data = [app.DataStruct.MCS.rel.R .*100, app.DataStruct.MCS.rel.ConfUp .*100,...
                                             app.DataStruct.MCS.rel.ConfDown .*100, app.DataStruct.model.MCSamples];
                                         
                     if any(isnan(app.DataStruct.MCS.rel.ConfUp)) || any(isnan(app.DataStruct.MCS.rel.ConfDown))
                        UpdateOutputWindow(app,'Number of Monte Carlo samples might be low to yield reliable results.',1)
                        uialert(app.ReliaComp,'Number of Monte Carlo samples might be low to yield reliable results!','MCS samples low','Icon','warning')
                     end
                    
                     assignin('base','Results',app.DataStruct);
                    
                     
                case 'MVFSOM'
                     
                    % Verifying input information related to HLRF!
                    UpdateOutputWindow(app,'Verifying input related to MVFSOM.......',1)
                    if 0
                        % Need to add verificaitons here!, warning about
                        % uisng non normal distributions.
                    else
                        UpdateOutputWindow(app,'Verifying input related to MVFSOM.......PASSED!',0)
                    end
                    
                    
                    % Run the HL-RF Reliability analysis !!!!! need errror
                    % handling here !!!
                    UpdateOutputWindow(app,'Begining Computation......',1)
                    app.DataStruct.MVFSOM.rel = findReliabilityMVFSOM(app);
                    
                    
                    % Changing the tab to the Output......!!!!
                    app.TabGroup.SelectedTab = app.ResultsTab;
                    tempTable = struct2table(app.DataStruct.MVFSOM.rel,'AsArray',true);
                    ResTable = table(tempTable.beta,tempTable.R .* 100,tempTable.Pf .* 100);
                    app.ResultsTable.Data = ResTable;
                    app.ResultsTable.ColumnName = {'Beta (-)'; 'Reliability (%)'; 'Probability of Failure (%)'};
                    
                    % Can be done in a better way, this is a temp fix!
                    app.ResultsTable.ColumnName(end+1:end+app.NoOfVarsSpinner.Value) =...
                        table2cell(app.InputTable.Data(:,1));
                    tabLen = length(app.ResultsTable.ColumnName);
                    for jdx = 1:length(app.DataStruct.MVFSOM.rel)
                        app.ResultsTable.Data(jdx,tabLen-app.NoOfVarsSpinner.Value+1:tabLen) = num2cell(tempTable.xStar{jdx})';
                    end
                    
                    app.ResultsTable.RowStriping = true;
                    
                    app.ResultsTable.BackgroundColor = repmat([1 1 1],[size(app.ResultsTable.Data,1),1]);
                    app.ResultsTable.BackgroundColor(end,:) = [0.5 1 0.5];
                    
                    assignin('base','Results',app.DataStruct);
                     
                otherwise
                    UpdateOutputWindow(app,'Invalid method! Please select a method that has been implimented.',1)
                    uialert(app.ReliaComp,'Selected method is not yet impimented!.','Method Error','Icon','warning')
                    return
            end
            
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, event)
            delete(app)
        end

        % Menu selected function: ResetFieldsMenu
        function ResetFieldsMenuSelected(app, event)
            % Clearing all fields and memory!
            app.FuncInputField.Value='';
            app.NoOfVarsSpinner.Value=0;
            app.ComputationMethodDropDown.Value='HL-RF';
            app.GradDropDown.Value ='Central Diff';
            app.GradPertubEditField.Value = 1e-6;
            app.AlphaEditField.Value = 1e-6;
            app.BetaEditField.Value = 1e-9;
            app.InputTable.Data = [];
            app.CorrTable.Data = [];
            app.OutputTextArea.Value={' '};
            app.ResultsTable.Data = [];
            app.DataStruct = struct();
            clc
        end

        % Menu selected function: SaveMenu
        function SaveMenuSelected(app, event)
            reliaCompData = app.DataStruct;
            DefaultFile = 'NewProject.mat'; ifile=1;
            while exist(DefaultFile,'file') == 2
                DefaultFile = sprintf('NewProject_%d.mat',ifile);
                ifile = ifile+1;
            end
            uisave({'reliaCompData'},DefaultFile)
        end

        % Menu selected function: DocumentationMenu
        function DocumentationMenuSelected(app, event)
            winopen(fullfile('docs', 'ReliaComp.pdf'));
        end

        % Menu selected function: LoadMenu
        function LoadMenuSelected(app, event)
            
            try
                uiopen('load')
                
                app.DataStruct = reliaCompData;
                
                % Common Data
                app.FuncInputField.Value            = app.DataStruct.input.Func;
                app.NoOfVarsSpinner.Value           = app.DataStruct.input.NoOfVars;
                app.ComputationMethodDropDown.Value = app.DataStruct.input.CompMethod;
                app.InputTable.Data                 = app.DataStruct.input.VarInTable;
                app.CorrTable.Data                  = app.DataStruct.input.VarCorrTable;
                %app.OutputTextArea.Value            = {' '};
                app.ResultsTable.Data               = [];
                
                % HLRF Related
                app.GradDropDown.Value              = app.DataStruct.input.GradMethod;
                app.GradPertubEditField.Value       = app.DataStruct.input.GradPeturb;
                app.AlphaEditField.Value            = app.DataStruct.input.Alpha;
                app.BetaEditField.Value             = app.DataStruct.input.Beta;
                
                % MCS Relate
                app.samplesforMCEditField.Value              = app.DataStruct.input.MCsamples;
                
                
                % HLRF??
                
                
                
                % Switching to input tab
                Switch2InputTab(app);
                
                % Adding the load path to MATLAB path
                addpath(app.DataStruct.input.Path);
                UpdateOutputWindow(app,'**************** Finished Loading File ****************',1)
                uialert(app.ReliaComp,'The file was sucessfully loaded!.','Success!!','Icon','success')
                
                
            catch ME
                
                uialert(app.ReliaComp,'There was an error loading the file, try loading a different file.','Error')
                return
                
            end
        end

        % Value changed function: ComputationMethodDropDown
        function ComputationMethodDropDownValueChanged(app, event)
            value = app.ComputationMethodDropDown.Value;
            app.DataStruct.input.CompMethod = value;
            
            switch value
                
                case 'HL-RF'
                    
                    app.TabGroup2.SelectedTab = app.HLRFTab;
                    
                case 'Monte Carlo'
                    app.TabGroup2.SelectedTab = app.MonteCarloTab;
                    
                case 'MVFOSM'
                    
                    app.TabGroup2.SelectedTab = app.MVFOSMTab;
                
            end
            
            
            
        end

        % Value changed function: GradDropDown
        function GradDropDownValueChanged(app, event)
            value = app.GradDropDown.Value;
            app.DataStruct.input.GradMethod = value;
        end

        % Value changed function: GradPertubEditField
        function GradPertubEditFieldValueChanged(app, event)
            value = app.GradPertubEditField.Value;
            app.DataStruct.input.GradPeturb = value;
        end

        % Value changed function: AlphaEditField
        function AlphaEditFieldValueChanged(app, event)
            value = app.AlphaEditField.Value;
            app.DataStruct.input.Alpha = value;
        end

        % Value changed function: BetaEditField
        function BetaEditFieldValueChanged(app, event)
            value = app.BetaEditField.Value;
            app.DataStruct.input.Beta = value;
        end

        % Menu selected function: AboutMenu
        function AboutMenuSelected(app, event)
            web("https://github.com/mayankchetan/ReliaComp")
        end

        % Cell edit callback: CorrTable
        function CorrTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
            % This kind of solves the need to validate the correlation
            % matrix!
            
            if indices(1) ~= indices(2) % Forcing the symmetry of the Correlation Matrix!
                app.CorrTable.Data(indices(2),indices(1)) = newData;
            else % Forcing the Diagonal elements to be one! LOL!
                app.CorrTable.Data(indices(1),indices(2)) = 1;
            end

        end

        % Value changed function: samplesforMCEditField
        function samplesforMCEditFieldValueChanged(app, event)
            value = app.samplesforMCEditField.Value;
            app.DataStruct.input.MCsamples = value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ReliaComp and hide until all components are created
            app.ReliaComp = uifigure('Visible', 'off');
            app.ReliaComp.Position = [100 100 884 653];
            app.ReliaComp.Name = 'ReliaComp v1.0';

            % Create Menu
            app.Menu = uimenu(app.ReliaComp);
            app.Menu.Text = 'Menu';

            % Create ResetFieldsMenu
            app.ResetFieldsMenu = uimenu(app.Menu);
            app.ResetFieldsMenu.MenuSelectedFcn = createCallbackFcn(app, @ResetFieldsMenuSelected, true);
            app.ResetFieldsMenu.Text = 'Reset Fields';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.Menu);
            app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadMenuSelected, true);
            app.LoadMenu.Text = 'Load';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.Menu);
            app.SaveMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveMenuSelected, true);
            app.SaveMenu.Text = 'Save';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.Menu);
            app.ExitMenu.MenuSelectedFcn = createCallbackFcn(app, @ExitMenuSelected, true);
            app.ExitMenu.Text = 'Exit';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.ReliaComp);
            app.HelpMenu.Text = 'Help';

            % Create DocumentationMenu
            app.DocumentationMenu = uimenu(app.HelpMenu);
            app.DocumentationMenu.MenuSelectedFcn = createCallbackFcn(app, @DocumentationMenuSelected, true);
            app.DocumentationMenu.Text = 'Documentation';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.MenuSelectedFcn = createCallbackFcn(app, @AboutMenuSelected, true);
            app.AboutMenu.Text = 'About';

            % Create ControlPanel
            app.ControlPanel = uipanel(app.ReliaComp);
            app.ControlPanel.Title = 'Control Panel';
            app.ControlPanel.Position = [11 399 344 245];

            % Create NumberofVariablesSpinnerLabel
            app.NumberofVariablesSpinnerLabel = uilabel(app.ControlPanel);
            app.NumberofVariablesSpinnerLabel.HorizontalAlignment = 'right';
            app.NumberofVariablesSpinnerLabel.Position = [21 119 114 22];
            app.NumberofVariablesSpinnerLabel.Text = 'Number of Variables';

            % Create NoOfVarsSpinner
            app.NoOfVarsSpinner = uispinner(app.ControlPanel);
            app.NoOfVarsSpinner.Limits = [0 Inf];
            app.NoOfVarsSpinner.ValueChangedFcn = createCallbackFcn(app, @NoOfVarsSpinnerValueChanged, true);
            app.NoOfVarsSpinner.Position = [150 119 100 22];

            % Create FuncInputField
            app.FuncInputField = uieditfield(app.ControlPanel, 'text');
            app.FuncInputField.Position = [150 171 140 22];

            % Create SelectFunctionButton
            app.SelectFunctionButton = uibutton(app.ControlPanel, 'push');
            app.SelectFunctionButton.ButtonPushedFcn = createCallbackFcn(app, @SelectFunctionButtonPushed, true);
            app.SelectFunctionButton.Position = [21 171 100 22];
            app.SelectFunctionButton.Text = 'Select Function';

            % Create ComputationMethodDropDownLabel
            app.ComputationMethodDropDownLabel = uilabel(app.ControlPanel);
            app.ComputationMethodDropDownLabel.HorizontalAlignment = 'right';
            app.ComputationMethodDropDownLabel.Position = [21 67 117 22];
            app.ComputationMethodDropDownLabel.Text = 'Computation Method';

            % Create ComputationMethodDropDown
            app.ComputationMethodDropDown = uidropdown(app.ControlPanel);
            app.ComputationMethodDropDown.Items = {'HL-RF', 'Monte Carlo', 'MVFSOM'};
            app.ComputationMethodDropDown.ValueChangedFcn = createCallbackFcn(app, @ComputationMethodDropDownValueChanged, true);
            app.ComputationMethodDropDown.Position = [150 67 158 22];
            app.ComputationMethodDropDown.Value = 'HL-RF';

            % Create RUNButton
            app.RUNButton = uibutton(app.ControlPanel, 'push');
            app.RUNButton.ButtonPushedFcn = createCallbackFcn(app, @RUNButtonPushed, true);
            app.RUNButton.Position = [122 16 100 22];
            app.RUNButton.Text = 'RUN!';

            % Create MethodSpecificInputsPanel
            app.MethodSpecificInputsPanel = uipanel(app.ReliaComp);
            app.MethodSpecificInputsPanel.Title = 'Method Specific Inputs';
            app.MethodSpecificInputsPanel.Position = [367 399 505 245];

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.MethodSpecificInputsPanel);
            app.TabGroup2.Position = [0 0 505 225];

            % Create HLRFTab
            app.HLRFTab = uitab(app.TabGroup2);
            app.HLRFTab.Title = 'HL-RF';

            % Create GradientMethodDropDownLabel
            app.GradientMethodDropDownLabel = uilabel(app.HLRFTab);
            app.GradientMethodDropDownLabel.HorizontalAlignment = 'right';
            app.GradientMethodDropDownLabel.Position = [35 149 95 22];
            app.GradientMethodDropDownLabel.Text = 'Gradient Method';

            % Create GradDropDown
            app.GradDropDown = uidropdown(app.HLRFTab);
            app.GradDropDown.Items = {'Central Diff', 'User Provided'};
            app.GradDropDown.ValueChangedFcn = createCallbackFcn(app, @GradDropDownValueChanged, true);
            app.GradDropDown.Position = [165 148 124 22];
            app.GradDropDown.Value = 'Central Diff';

            % Create ConvergenceAlphaEditFieldLabel
            app.ConvergenceAlphaEditFieldLabel = uilabel(app.HLRFTab);
            app.ConvergenceAlphaEditFieldLabel.HorizontalAlignment = 'right';
            app.ConvergenceAlphaEditFieldLabel.Position = [35 65 111 22];
            app.ConvergenceAlphaEditFieldLabel.Text = 'Convergence Alpha';

            % Create AlphaEditField
            app.AlphaEditField = uieditfield(app.HLRFTab, 'numeric');
            app.AlphaEditField.ValueChangedFcn = createCallbackFcn(app, @AlphaEditFieldValueChanged, true);
            app.AlphaEditField.Position = [165 65 100 22];
            app.AlphaEditField.Value = 1e-06;

            % Create ConvergenceBetaEditFieldLabel
            app.ConvergenceBetaEditFieldLabel = uilabel(app.HLRFTab);
            app.ConvergenceBetaEditFieldLabel.HorizontalAlignment = 'right';
            app.ConvergenceBetaEditFieldLabel.Position = [35 24 105 22];
            app.ConvergenceBetaEditFieldLabel.Text = 'Convergence Beta';

            % Create BetaEditField
            app.BetaEditField = uieditfield(app.HLRFTab, 'numeric');
            app.BetaEditField.ValueChangedFcn = createCallbackFcn(app, @BetaEditFieldValueChanged, true);
            app.BetaEditField.Position = [165 24 100 22];
            app.BetaEditField.Value = 1e-09;

            % Create GradientStepLabel
            app.GradientStepLabel = uilabel(app.HLRFTab);
            app.GradientStepLabel.HorizontalAlignment = 'right';
            app.GradientStepLabel.Position = [35 109 76 22];
            app.GradientStepLabel.Text = 'Gradient Step';

            % Create GradPertubEditField
            app.GradPertubEditField = uieditfield(app.HLRFTab, 'numeric');
            app.GradPertubEditField.ValueChangedFcn = createCallbackFcn(app, @GradPertubEditFieldValueChanged, true);
            app.GradPertubEditField.Position = [165 109 100 22];
            app.GradPertubEditField.Value = 1e-06;

            % Create MonteCarloTab
            app.MonteCarloTab = uitab(app.TabGroup2);
            app.MonteCarloTab.Title = 'Monte Carlo';

            % Create NumberofsamplesforMCEditFieldLabel
            app.NumberofsamplesforMCEditFieldLabel = uilabel(app.MonteCarloTab);
            app.NumberofsamplesforMCEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofsamplesforMCEditFieldLabel.Position = [98 101 149 22];
            app.NumberofsamplesforMCEditFieldLabel.Text = 'Number of samples for MC';

            % Create samplesforMCEditField
            app.samplesforMCEditField = uieditfield(app.MonteCarloTab, 'numeric');
            app.samplesforMCEditField.ValueChangedFcn = createCallbackFcn(app, @samplesforMCEditFieldValueChanged, true);
            app.samplesforMCEditField.Position = [262 101 100 22];
            app.samplesforMCEditField.Value = 1000000;

            % Create MVFOSMTab
            app.MVFOSMTab = uitab(app.TabGroup2);
            app.MVFOSMTab.Title = 'MVFOSM';

            % Create WarningThismethodassumesalltheLabel
            app.WarningThismethodassumesalltheLabel = uilabel(app.MVFOSMTab);
            app.WarningThismethodassumesalltheLabel.Position = [22 129 469 43];
            app.WarningThismethodassumesalltheLabel.Text = {'Warning: This method assumes all the input variables are independant and normally'; 'distributed'};

            % Create InteractionPanel
            app.InteractionPanel = uipanel(app.ReliaComp);
            app.InteractionPanel.Title = 'Interaction Panel';
            app.InteractionPanel.Position = [11 10 861 381];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.InteractionPanel);
            app.TabGroup.Position = [1 0 860 359];

            % Create InputsTab
            app.InputsTab = uitab(app.TabGroup);
            app.InputsTab.Title = 'Variable Inputs';

            % Create InputTable
            app.InputTable = uitable(app.InputsTab);
            app.InputTable.ColumnName = {'Variable Name'; 'Mean'; 'COV'; 'Distributions'};
            app.InputTable.RowName = {''};
            app.InputTable.ColumnSortable = false;
            app.InputTable.ColumnEditable = true;
            app.InputTable.Position = [11 9 573 315];

            % Create CorrTable
            app.CorrTable = uitable(app.InputsTab);
            app.CorrTable.ColumnName = {''};
            app.CorrTable.ColumnWidth = {'auto'};
            app.CorrTable.RowName = {''};
            app.CorrTable.ColumnSortable = false;
            app.CorrTable.ColumnEditable = true;
            app.CorrTable.CellEditCallback = createCallbackFcn(app, @CorrTableCellEdit, true);
            app.CorrTable.Position = [592 36 251 185];

            % Create CorrelationMatrixLabel
            app.CorrelationMatrixLabel = uilabel(app.InputsTab);
            app.CorrelationMatrixLabel.HorizontalAlignment = 'center';
            app.CorrelationMatrixLabel.Position = [668 229 100 22];
            app.CorrelationMatrixLabel.Text = 'Correlation Matrix';

            % Create OutputTab
            app.OutputTab = uitab(app.TabGroup);
            app.OutputTab.Title = 'Output';

            % Create OutputTextArea
            app.OutputTextArea = uitextarea(app.OutputTab);
            app.OutputTextArea.Position = [11 9 835 315];

            % Create ResultsTab
            app.ResultsTab = uitab(app.TabGroup);
            app.ResultsTab.Title = 'Results';

            % Create ResultsTable
            app.ResultsTable = uitable(app.ResultsTab);
            app.ResultsTable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.ResultsTable.RowName = {};
            app.ResultsTable.Position = [10 10 836 319];

            % Show the figure after all components are created
            app.ReliaComp.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ReliaCompGUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ReliaComp)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ReliaComp)
        end
    end
end