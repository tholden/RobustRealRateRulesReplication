% This script is only used to build the release package. It requires PANDOC (https://pandoc.org/) to be installed in the system path.

close all;
clear variables;

warning( 'off', 'MATLAB:DELETE:FileNotFound' );
try
    delete Release.zip;
catch
end
try
    delete private/autocorrNoMean.m;
catch
end
warning( 'on', 'MATLAB:DELETE:FileNotFound' );

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

copyfile Inputs/*.xlsx Release/Inputs/ f;
copyfile private/*.m Release/private/ f;
copyfile Main.m Release/ f;
copyfile *.md Release/ f;

!pandoc -f gfm -t pdf -o Release/README.pdf README.md

cd Release/;
zip ../Release.zip *;
cd ..;

try
    rmdir Release s;
catch
end
