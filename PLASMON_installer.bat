@echo off
set python_version=3.8.6
set python_url=https://www.python.org/ftp/python/%python_version%/python-%python_version%.exe
goto CheckPermissions

:CheckPermissions
net session >nul 2>&1
if %errorLevel% == 0 (
	echo Success: Administrative permissions confirmed
	goto CheckGit
) else (
	echo Failure: Current permissions inadequate. Run as administrator
	pause
	exit /b 1
)

:CheckGit
git version >NUL
if %errorLevel% == 0 (
	echo Success: You have Git istalled
	goto PythonCheck
) else (
	echo Failure: Git not installed, which I cannot do from this script. Please install Git on your PC, you can use all the standard settings. Go here: https://gitforwindows.org/
	pause
	exit /b 1
)

:PythonCheck
python -V >NUL
if %errorLevel% == 0 (
	set py_call=python
	goto AskDirectory
) else (
	goto PythonCheck2
)

:PythonCheck2
py -V >NUL
if %errorLevel% == 0 (
	set py_call=py
	goto AskDirectory
) else (
	echo Failure: You don't have Python installed. Will download it now.
	goto PythonDownload
) 

:PythonDownload
powershell -Command "(New-Object Net.WebClient).DownloadFile('%python_url%', '%cd%/PythonInstaller.exe')"
if %errorLevel% == 0 (
	echo Success: Python downloaded.
	goto PythonInstall
) else (
	echo Failure: Cannot download Python. Check your internet. If you are doing this before 2023 this should work. If it really doesnt, contact the administrator.
	pause
	exit /b 1
)

:PythonInstall
powershell -Command "%cd%/PythonInstaller.exe -quiet InstallAllUsers=1 Include_test=0"
if %errorLevel% == 0 (
	echo Success: Python installed.
	goto PythonCheck
) else (
	echo Failure:  Cannot install Python. Check if you don't have any remains left and try again.
	pause
	exit /b 1
)

:AskDirectory
echo Success: Python ready to use. 
echo Enter where you want the program to be installed (copy the full route to the directory).
echo Be sure to put "" around the directory. You can also select your current installation directory, the program will update it.
set /p directory="Your desired PLASMON directory (program will check for a PLASMON folder in this directory): "
set directory=%directory%\PLASMON
echo Chosen directory: %directory%
goto CheckDirectory

:CheckDirectory
if exist %directory% (
	echo Directory exists. Will update directory
	goto Update
) else (
	echo Directory does not exist, will make directory and download PLASMON there.
	goto Clone
)

:Clone
git clone https://github.com/DionEngels/PLASMON.git %directory%
if %errorLevel% == 0 (
	echo Success: Directory cloned.
	goto PythonVirtualEnv
) else (
	echo Failure: Cannot clone PLASMON into that directory. Check internet and also check if the directory does not already exist. It has be to be empty.
	pause
	exit /b 1
)

:PythonVirtualEnv
for /f %%i in ('%py_call% -c "import sys; import os; print(os.path.dirname(sys.executable))"') do set py_dir=%%i
%py_call% -m venv %py_dir%\PLASMON_venv
REM powershell -Command "[Environment]::SetEnvironmentVariable('PATH', '${env:path};${%directory%\PLASMON\venv}', 'Machine')"
if %errorLevel% == 0 (
	echo Success: Virtual environment created in Python installation
	goto PythonVirtualEnvActivate
) else (
	echo Failure: Cannot create virtual environment. Contact administrator.
	pause
	exit /b 1
)

:Update
cd %directory%
git fetch
git pull
if %errorLevel% == 0 (
	echo Success: Updated directory
	goto PythonCheckVirtualEnv
) else (
	echo Failure: Cannot update directory. Contact administrator.
	pause
	exit /b 1
)

:PythonCheckVirtualEnv
for /f %%i in ('%py_call% -c "import sys; import os; print(os.path.dirname(sys.executable))"') do set py_dir=%%i

if exist %py_dir%\PLASMON_venv (
	echo Virtual environment found. Will activate and update.
	goto PythonVirtualEnvActivate
) else (
	echo Virtual environment not found. Will create.
	goto PythonVirtualEnv
)

:PythonVirtualEnvActivate
cd %directory%\PLASMON
%py_dir%\PLASMON_venv\Scripts\activate.bat & pip install -r requirements.txt > NUL & color 07 & echo Success: Enabled virtual environment and updated packages. Installation complete. & pause
if not %errorLevel% == 0 (
	echo Cannot enable virtual environment. Contact administrator.
	pause
	exit /b 1
)
