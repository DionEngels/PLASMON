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