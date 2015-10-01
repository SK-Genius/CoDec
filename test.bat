ECHO OFF
SET comp=8

SET list=LenaColor
REM SET list=bird circles crosses camera slope text bridge goldhill test test2 Color18V Siemensstern Lena LenaColor

FOR %%n IN (%list%) DO (
	del bin\%%n.dat bin\%%n_.bmp
)

ECHO %comp%

FOR %%n IN (%list%) DO (
	echo --------------------------------------------------------------------------------
	echo %%n
	bin\main.exe -e%comp% bin/%%n.bmp bin/%%n.dat
	ls -s ./bin/%%n.dat
	bin\main.exe -d bin/%%n.dat bin/%%n_.bmp
)

