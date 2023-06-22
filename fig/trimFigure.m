function trimFigure(fig)
%TRIMFIGURE Trim figure borders for faithful PDF export
%   Sets proper figure borders in order to ensure optimal PDF output.
%
%   Usage:
%
%     TRIMFIGURE, by itself, trims the current figure.
%     TRIMFIGURE(FIG) trims the figure with handle FIG.
% 
%   Examples:
%     - Plot random data, trim, and save to PDF.
%         plot(randn(1,10));
%         trimFigure;
%         saveas(gcf,'randn.pdf');
%
%   See also FIGURE, SAVEAS, PRINT.

% (c) 2019-2023 Matteo Seclì <secli.matteo@gmail.com>. All Rights Reserved.“
% Last revision: 13-Apr-2023

if nargin < 1
    fig = gcf;
end

fig.PaperPositionMode = 'auto';
fig.PaperSize = fig.PaperPosition(3:4);

end

