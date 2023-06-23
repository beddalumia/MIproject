function saveFigureAsPDF(fig,filename,trim,transp)
%SAVEFIGUREASPDF Save a figure as PDF, with optimal borders and background
% 
%   SAVEFIGUREASPDF(FIG, FILENAME) saves the figure with handle FIG as a
%     PDF file called FILENAME. The figure is temporarily trimmed in
%     order to ensure optimal borders, and it's provided with a
%     transparent background.
%   SAVEFIGUREASPDF(FIG,FILENAME,TRIM,TRANSP) saves the figure with handle
%     FIG as a PDF file called FILENAME. The optional logical inputs TRIM
%     and TRANSP enable/disable the trim operation and the transparent
%     background color, respectively. If TRANSP is set to false, in
%     particular, the figure will use a white background.
% 
%   Example:
%     Plot random data and save to PDF.
%       plot(randn(1,10));
%       SAVEFIGUREASPDF(gcf,'randn.pdf');
%
%   See also TRIMFIGURE, EXPORTGRAPHICS, SAVEAS, PRINT. 

% (c) 2023 Matteo Seclì <secli.matteo@gmail.com>. All Rights Reserved.“
% Last revision: 15-Apr-2023

if nargin < 4
    transp = true;
    if nargin < 3
        trim = true;
        if nargin < 2
            eid = 'saveFigureAsPDF:tooFewInputs';
            msg = 'Too few inputs.';
            throw(MException(eid,msg));
        end
    end
end

if ~all(ishandle(fig))
    eid = 'saveFigureAsPDF:invalidHandle';
    msg = 'Invalid handle.';
    throw(MException(eid,msg));
end

% Save current settings
prevPaperPositionMode = fig.PaperPositionMode;
prevPaperSize         = fig.PaperSize;
prevColor             = fig.Color;
prevInvertHardcopy    = fig.InvertHardcopy;

% Trim
if trim
    fig.PaperPositionMode = 'auto';
    fig.PaperSize = fig.PaperPosition(3:4);
end

% Choose background color
if transp
    figColor = 'None';
else
    figColor = 'white';
end

% Temporarily set the new color settings
fig.Color          = figColor;
fig.InvertHardcopy = 'Off';

% Save while temporarily deactivating a print warning
w_id = 'MATLAB:print:FigureTooLargeForPage';
w    = warning('query',w_id);
warning('off',w_id);
saveas(fig,filename,'pdf');
warning(w.state,w_id);

% Restore old settings
fig.PaperPositionMode = prevPaperPositionMode;
fig.PaperSize         = prevPaperSize;
fig.Color             = prevColor;
fig.InvertHardcopy    = prevInvertHardcopy;

end

