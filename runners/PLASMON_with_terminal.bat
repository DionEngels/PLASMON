::[Bat To Exe Converter]
::
::YAwzoRdxOk+EWAjk
::fBw5plQjdCyDJGyX8VAjFB9dQwqHAES0A5EO4f7+r6fHlnhdcesxfIfUzLGAHNcW+E7YWLQ1lkpDjMMDAltZcBOndxw9ulJhuWCAC86fvEHoSUfp
::YAwzuBVtJxjWCl3EqQJgSA==
::ZR4luwNxJguZRRnk
::Yhs/ulQjdF+5
::cxAkpRVqdFKZSDk=
::cBs/ulQjdF+5
::ZR41oxFsdFKZSDk=
::eBoioBt6dFKZSDk=
::cRo6pxp7LAbNWATEpCI=
::egkzugNsPRvcWATEpCI=
::dAsiuh18IRvcCxnZtBJQ
::cRYluBh/LU+EWAnk
::YxY4rhs+aU+JeA==
::cxY6rQJ7JhzQF1fEqQJQ
::ZQ05rAF9IBncCkqN+0xwdVs0
::ZQ05rAF9IAHYFVzEqQJQ
::eg0/rx1wNQPfEVWB+kM9LVsJDGQ=
::fBEirQZwNQPfEVWB+kM9LVsJDGQ=
::cRolqwZ3JBvQF1fEqQJQ
::dhA7uBVwLU+EWDk=
::YQ03rBFzNR3SWATElA==
::dhAmsQZ3MwfNWATElA==
::ZQ0/vhVqMQ3MEVWAtB9wSA==
::Zg8zqx1/OA3MEVWAtB9wSA==
::dhA7pRFwIByZRRnk
::Zh4grVQjdCyDJGyX8VAjFB9dQwqHAES0A5EO4f7+r6fHlnhdcesxfIfUzLGAHNcW+E7YWLQ1lkpDjMMDAltZcBOndxw9ulJhuWCAC/ewkjzMZWa28UIkE1pggnHdhSU6bJ1tgsZj
::YB416Ek+ZG8=
::
::
::978f952a14a936cc963da21a135fa983
@echo off
:PythonCheck
python -V >NUL
if %errorLevel% == 0 (
	set py_call=python
	goto Run
) else (
	goto PythonCheck2
)

:PythonCheck2
py -V >NUL
if %errorLevel% == 0 (
	set py_call=py
	goto Run
) else (
	echo Failure: You don't have Python installed. Run the installer script.
	pause
	exit /b 1
) 

:Run
%cd%\venv\Scripts\activate.bat & %py_call% main_gui.py & exit /b 0
if %errorLevel% == 0 (
	exit /b 0
) else (
	echo Failure: Cannot run program. Rerun your installer if you touched the venv folder. Otherwise contact administrator.
	pause
	exit /b 1
) 