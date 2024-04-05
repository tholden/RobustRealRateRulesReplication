% This script is only used to build the release package. It requires PANDOC (https://pandoc.org/) and a LATEX interpreter (like MIKTEX, https://miktex.org/) to be installed in the system path.

close all;
clear variables;

warning( 'off', 'MATLAB:DELETE:FileNotFound' );
try
    delete Release.zip;
catch
end
try
    delete private/autocorrNoMean.m; % This is generated automatically by Main.m based on a file distributed with MATLAB.
catch
end
warning( 'on', 'MATLAB:DELETE:FileNotFound' );

warning( 'off', 'MATLAB:MKDIR:DirectoryExists' );
try
    rmdir Release s;
catch
end
try
    mkdir Release;
catch
end
try
    mkdir Release/Inputs;
catch
end
try
    mkdir Release/private;
catch
end

% Copy the core files.
copyfile Inputs/*.xlsx Release/Inputs/ f;
copyfile private/*.m Release/private/ f;
copyfile Main.m Release/ f;
copyfile *.md Release/ f;

% Compile the PDF version of the README.
!pandoc -f gfm -t pdf --pdf-engine=xelatex --shift-heading-level-by=-1 -V lang=en-GB -V boxlinks -V hyperrefoptions:pdfborderstyle="{/S/U/W 1}" -V hyperrefoptions:allbordercolors="{0 0 0}" -V papersize=a4 -V geometry:margin=1.25in -V mainfont="TeX Gyre Pagella" -o Release/README.pdf README.md

% Produce CSV versions of all XLSX spreadsheets.
try
    mkdir Release/InputsAsCSV;
catch
end

FileList = dir( 'Inputs' );

for i = 1 : numel( FileList )

    FileName = FileList( i ).name;

    if ~endsWith( FileName, '.xlsx', 'IgnoreCase', true )
        continue
    end

    NewFileOrFolderName = [ 'Release/InputsAsCSV/' FileName( 1 : ( end - 5 ) ) ];

    FilePath = [ 'Inputs/' FileList( i ).name ];

    SheetNames = sheetnames( FilePath );

    if ~isscalar( SheetNames )

        try
            mkdir( NewFileOrFolderName );
        catch
        end

    end

    for j = 1 : numel( SheetNames )

        SheetName = SheetNames{ j };

        Table = readtable( FilePath, 'Sheet', SheetName, 'VariableNamingRule', 'preserve' );

        if isscalar( SheetNames )
            Separator = '';
        else
            Separator = [ '/' SheetName ];
        end

        writetable( Table, [ NewFileOrFolderName Separator '.csv' ], 'WriteVariableNames', ~all( startsWith( Table.Properties.VariableNames, 'Var', 'IgnoreCase', false ) ) );

    end

end

cd Release/;
zip ../Release.zip *;
cd ..;

try
    rmdir Release s;
catch
end

warning( 'on', 'MATLAB:MKDIR:DirectoryExists' );
