# PLASMON

Created by Dion Engels. Contact: d.j.engels@student.tue.nl
Program with full GUI to analyze time traces and HSM on .nd2 files. No Python knowledge needed, all work is done in the GUI.

Created to work together with SPectrA. SPectrA is a MATLAB app to view and post-process the results of PLASMON.

Link to SPectrA: https://github.com/DionEngels/SPectrA

Installation Instructions
- 
- Go to releases here on Github.
- Download the newest release.
- Install via installer.
- Check out the guide that is also included in the release! This has a lot of useful information.
- Open PLASMON.exe from the installation directory.
- Done!

Troubleshooting
-
#### In case the program does not open
- The program has Microsoft's C++ redistributions as a requirement. Usually you have those, but if the program does not open, please install those.
- Check which version you need for your PC. It depends on your hardware. 
- To check this: search "This PC".
- Right click "Properties".
- Check if the processor is named as x64 or x86 (or something else).
- Go to the Microsoft download site for Visual C++. At point of writing this is: https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads.
- Download the appropriate versions. You need all possible versions. Just start at the top (2019 at point of writing), and work down until your architecture is no longer supported.
- The program should now run.
- In case it does not, result to the general error procedure below.

#### How to report an error during program runtime
- Make a screenshot of the error message.
- Go to My Documents/PLASMON and get the Logging.txt file.
- If the error occured after your pressed "Run", get the Settings.txt file from your results directory.
- Make a WeTransfer or OneDrive transfer of your datafile.
- Send this all to d.j.engels@student.tue.nl with any other useful information about the error.
