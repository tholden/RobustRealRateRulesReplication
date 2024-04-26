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
try
    delete private/dynareNoMean.m; % This is generated automatically by Main.m based on a file distributed with Dynare.
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
    mkdir Release/SmetsWouters2007;
catch
end
try
    mkdir Release/private;
catch
end

% Copy the core files.
copyfile Inputs/*.xlsx Release/Inputs/ f;
copyfile SmetsWouters2007/*.mod Release/SmetsWouters2007/ f;
copyfile SmetsWouters2007/*.mat Release/SmetsWouters2007/ f;
copyfile SmetsWouters2007/*.txt Release/SmetsWouters2007/ f;
copyfile private/*.m Release/private/ f;
copyfile Main.m Release/ f;
copyfile *.md Release/ f;

% Compile the PDF version of the README.
!pandoc -f gfm -t pdf --pdf-engine=xelatex --shift-heading-level-by=-1 -V lang=en-GB -V boxlinks -V hyperrefoptions:pdfborderstyle="{/S/U/W 1}" -V hyperrefoptions:allbordercolors="{0 0 0}" -V papersize=a4 -V geometry:margin=1.25in -V mainfont="TeX Gyre Pagella" -o Release/README.pdf README.md

cd Release/;
zip ../Release.zip *;
cd ..;

try
    rmdir Release s;
catch
end

warning( 'on', 'MATLAB:MKDIR:DirectoryExists' );
