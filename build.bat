@ECHO OFF

SET CompilerFlags=
SET CompilerFlags=%CompilerFlags% -std=c++0x
SET CompilerFlags=%CompilerFlags% -O3
SET CompilerFlags=%CompilerFlags% -D TIMER
SET CompilerFlags=%CompilerFlags% -g3
REM SET CompilerFlags=%CompilerFlags% -Wall
REM SET CompilerFlags=%CompilerFlags% -Werror
SET CompilerFlags=%CompilerFlags% -I"C:/MinGW/lib/gcc/mingw32/4.8.1/include"
SET CompilerFlags=%CompilerFlags% -I"inc"

SET LinkerFlags=
SET LinkerFlags=%LinkerFlags% -fno-use-linker-plugin
SET LinkerFlags=%LinkerFlags% -B"C:/MinGW/lib/gcc/mingw32/4.8.1/"

SET SourceFiles=
SET SourceFiles=%SourceFiles% "src/main.cpp"

@ECHO ON

g++ %CompilerFlags% %LinkerFlags% -o"bin/main.exe" %SourceFiles% && "bin/main.exe"
