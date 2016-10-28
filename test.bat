@ECHO OFF

SET comp=0

SET list=bird circles crosses camera slope text bridge goldhill test test2 Color18V Siemensstern LenaColor Lena

FOR %%n IN (%list%) DO (
	del bin\%%n.dat bin\%%n_.bmp
)

ECHO %comp%

FOR %%n IN (%list%) DO (
	echo --------------------------------------------------------------------------------
	echo %%n
	SET X=%%n
	bin\main.exe -e%comp% bin/%%n.bmp bin/%%n.dat
	dir /B bin\%%n.dat
	bin\main.exe -d bin/%%n.dat bin/%%n_.bmp
)
