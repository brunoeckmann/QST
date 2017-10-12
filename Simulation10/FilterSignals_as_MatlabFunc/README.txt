README - HOW TO USE filter_signals.cpp AS MATLAB-FUNCTION
=========================================================
1) Check if mex is running
--------------------------
Compile with the following command

    mex HelloWorld.cpp

the HelloWorld.cpp file and check if it prints the output proberly.

2) Install boost on machine
---------------------------
Download and install the C++ Boost Libs from
    
    http://www.boost.org/users/download/

3) Test Boost in Terminal
---------------------------
Go to terminal and compile with

    c++ TestBoost.cpp

the Boost test and check if the version is printed properly.

4) Compile filter_signals_mex.cpp
---------------------------------
Go to MATLAB terminal and compile filter_signals_mex.cpp via command window with
    mex file.cpp -I/usr/local/include
for Mac with the Boost lib in the include folder and for Windows via
    mex file.cpp -IC:\pathToBoostFolder
where the boost folder has to be specified (the programm needs the lib multi_array.hpp)

5) Run test_FilterSignals_3322.m in Matlab
------------------------------------------
The output of test_FilterSignals_3322.m must now be the same as the terminal output
of filter_signals.exe with behaviour_noise.txt as input