@echo off
:PythonCheck
python -V >NUL
if %errorLevel% == 0 (
	set py_call=python
	goto Clone
) else (
	goto PythonCheck2
)

:PythonCheck2
py -V >NUL
if %errorLevel% == 0 (
	set py_call=py
	goto PythonCheckDir
) else (
	echo Failure: You don't have Python installed. Run the installer script.
	exit /b 1
) 


:PythonCheckDir
for /f %%i in ('%py_call% -c "import sys; import os; print(os.path.dirname(sys.executable))"') do set py_dir=%%i
if %errorLevel% == 0 (
	goto Run
) else (
	echo Failure: Cannot find python installation dir. Please contact administrator.
	exit /b 1
) 

:Run
%py_dir%\PLASMON_venv\Scripts\activate.bat & %py_call% main_gui.py & exit
exit
if %errorLevel% == 0 (
	exit /b 0
) else (
	echo Failure: Cannot run program. Contact administrator.
	exit /b 1
) 