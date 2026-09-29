function saveFigAll(fig, basePath, formats)
% SAVEFIGALL Save a figure to multiple output formats in one call.
%
%   saveFigAll(fig, basePath)
%   saveFigAll(fig, basePath, formats)
%
%   Inputs:
%     fig      - figure handle (e.g. gcf)
%     basePath - full path + filename, WITHOUT extension
%                e.g. fullfile(saveDir_Fig4, [tv 'Example Excited cell ' num2str(n)])
%     formats  - (optional) cellstr of formats to save, any of:
%                'fig', 'png', 'svg'
%                Default: {'fig','png','svg'} (matches original script behavior)
%
%   Example:
%     saveFigAll(gcf, fullfile(saveDir_Fig4, [tv 'Example Excited cell ' ...
%         num2str(n) ' ' Params.ExpName{F}]));

    if nargin < 3
        formats = {'fig','png','svg'};
    end

    for k = 1:numel(formats)
        switch lower(formats{k})
            case 'fig'
                saveas(fig, basePath, 'fig');
            case 'png'
                % screen2png does not accept a figure handle in the
                % original codebase's usage, so keep the same call style
                screen2png(basePath);
            case 'svg'
                print(fig, basePath, '-dsvg');
            otherwise
                warning('saveFigAll:unknownFormat', ...
                    'Unknown format "%s" — skipped.', formats{k});
        end
    end
end