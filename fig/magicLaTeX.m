function varargout = magicLaTeX(options)
%MAGICLATEX Sets LaTeX everywhere
%   This function globally sets LaTeX everywhere, plus a few other settings
%   to produce paper-quality plots.
%
%   Usage:
%
%     cleanLaTeX = MAGICLATEX;    % Sets LaTeX everywhere
%         % --- do your plotting stuff here --- %
%     clear cleanLaTeX            % Restores previous settings
% 
%   If you don't specify an output variable, the settings are set globally
%   and are persistent until you clear MAGICLATEX's persistent variables
%   via (note that this is not the same variable as before):
% 
%     clear MAGICLATEX
% 
%   The syntax without an output variable is especially useful in
%   'startup.m'; also note that it resets the settings to the first-ever
%   settings detected by MAGICLATEX, not to the ones detected in the latest
%   function call.
%
%   Options:
% 
%     cleanLaTeX = MAGICLATEX(options);
% 
%   where 'options' are specified as MAGICLATEX(...,'key',value). Accepted
%   options are:
%     'FontSize'    Font size for axes and colorbar.
%                   Default: 16.
% 
%   See also TEX.

% (c) 2023 Matteo Seclì. All Rights Reserved.


%% Input handling
% Preliminary input validation
% Format:
%     argName (dimensions) dataType {validators} = defaultValue
arguments
    options.FontSize {mustBeInteger,mustBePositive} = 16
end


%% Define our persistent variables
persistent magicLaTeXCleanup;
persistent persistentOldDefaults;


%% Magic LaTeX everywhere
% Execute only if we want an output or we are calling the persistent
% variant and it's the first time ever we call the function.
if ( ( isempty(persistentOldDefaults) && nargout==0 ) || nargout~=0 )
    % Get current defaults
    grootfields = {'defaultAxesTickLabelInterpreter', ...
        'defaultAxesTickLabelInterpreter', ...
        'defaultAxesFontSizeMode', ...
        'defaultAxesFontSize', ...
        'DefaultAxesLineWidth', ...
        'defaultColorbarTickLabelInterpreter', ...
        'defaultColorbarFontSizeMode', ...
        'defaultColorbarFontSize', ...
        'defaultFigureUnits', ...
        'defaultLineLineWidth', ...
        'defaultLegendInterpreter', ...
        'defaultFigurePaperPositionMode', ...
        'defaultPolaraxesTickLabelInterpreter', ...
        'defaultTextInterpreter', ...
        'defaultTextarrowshapeInterpreter', ...
        'defaultTextboxshapeInterpreter', ...
        'defaultFigurePosition'};
    for grootfield=grootfields
        if verLessThan('matlab','9.3')
            % MATLAB R2017a and earlier
            oldDefaults = setfield(oldDefaults,grootfield{1},get(groot,grootfield{1}));
        else
            % MATLAB R2017b and later
            oldDefaults.(grootfield{1}) = get(groot,grootfield{1});
        end
    end
end

% Store the first-ever settings we get
if isempty(persistentOldDefaults)
    persistentOldDefaults = oldDefaults;
end

% Execute only if we want an output or we are calling the persistent
% variant and it's the first time we call the persistent variant (but not
% necessarily the first time we call the function overall).
if ( ( isempty(magicLaTeXCleanup) && nargout==0 ) || nargout~=0 )
    % LaTeX everywhere
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultAxesFontSizeMode', 'manual');
    set(groot, 'defaultAxesFontSize', options.FontSize);
    set(groot, 'DefaultAxesLineWidth', 1.0);
    set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
    set(groot, 'defaultColorbarFontSizeMode', 'manual');
    set(groot, 'defaultColorbarFontSize', options.FontSize);
    set(groot, 'defaultFigureUnits', 'points');
    set(groot, 'defaultLineLineWidth', 2.0);
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultFigurePaperPositionMode', 'auto');
    set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultTextarrowshapeInterpreter', 'latex');
    set(groot, 'defaultTextboxshapeInterpreter', 'latex');
    defaultFigurePosition = get(groot, 'defaultFigurePosition');
    set(groot, 'defaultFigurePosition', [defaultFigurePosition(1:2),560,420]);
    clear defaultFigurePosition;
end

% Contruct the output object, if any.
switch nargout
    case 0
        % Create a persistent cleanup object; magic LaTeX survives until
        % magicLaTeX itself is cleared.
        if isempty(magicLaTeXCleanup)
            magicLaTeXCleanup = onCleanup(@()cleanupFun(persistentOldDefaults));
        end
    case 1
        varargout = {onCleanup(@()cleanupFun(oldDefaults))};
    otherwise
        error('Too many output variables requested!');
end


end


%% Cleanup function
function cleanupFun(oldDefaults)
    set(groot,oldDefaults);
end