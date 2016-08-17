@ECHO OFF
setlocal

rem SET MINGW_PATH=C:/Program Files/mingw-w64/x86_64-6.1.0-posix-seh-rt_v5-rev0/mingw64/lib/gcc/x86_64-w64-mingw32/6.1.0/
rem 
rem SET CompilerFlags=
rem SET CompilerFlags=%CompilerFlags% -c
rem SET CompilerFlags=%CompilerFlags% -std=c++0x
rem SET CompilerFlags=%CompilerFlags% -O3
rem SET CompilerFlags=%CompilerFlags% -g3
rem SET CompilerFlags=%CompilerFlags% -Wall
rem SET CompilerFlags=%CompilerFlags% -Werror
rem SET CompilerFlags=%CompilerFlags% -I"%MINGW_PATH%include/"
rem SET CompilerFlags=%CompilerFlags% -fno-use-linker-plugin
rem SET CompilerFlags=%CompilerFlags% -DTIMER
rem 
rem g++ %CompilerFlags% -o"bin/main.o" "src/main.cpp" || exit /B

del "bin\main.exe"
del "main.obj"

SET VS_PATH=C:/Program Files (x86)/Microsoft Visual Studio 14.0/VC/

call "%VS_PATH%vcvarsall.bat" x86

SET CompilerFlags=
SET CompilerFlags=%CompilerFlags% /c
SET CompilerFlags=%CompilerFlags% /Od
SET CompilerFlags=%CompilerFlags% /WX
SET CompilerFlags=%CompilerFlags% /Wall
SET CompilerFlags=%CompilerFlags% /wd4100
SET CompilerFlags=%CompilerFlags% /wd4365
SET CompilerFlags=%CompilerFlags% /wd4514
SET CompilerFlags=%CompilerFlags% /D VS
SET CompilerFlags=%CompilerFlags% /D TEST
SET CompilerFlags=%CompilerFlags% /D TIMER

cl %CompilerFlags% /Fo: "bin\main.obj" "src\main.cpp" 2>&1
rem || exit /B

SET CompilerFlags=
SET CompilerFlags=%CompilerFlags% /O2
SET CompilerFlags=%CompilerFlags% /D VS
SET CompilerFlags=%CompilerFlags% /D TEST

cl %CompilerFlags% /Fe: "bin\main.exe" "src\main.cpp" || exit /B

"bin\main.exe" -t || exit /B

echo == SELF TEST OK ==========================

SET CompilerFlags=
SET CompilerFlags=%CompilerFlags% /Od
SET CompilerFlags=%CompilerFlags% /DEBUG
SET CompilerFlags=%CompilerFlags% /D VS
SET CompilerFlags=%CompilerFlags% /D TIMER

cl %CompilerFlags% /Fe: "bin\main.exe" "src\main.cpp" || exit /B
