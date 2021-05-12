Build with CMake
=====================

1. Clone repository
-------------------
::

   git clone https://github.com/andreabonvini/pointprocess.git

2. Build project
------------------

Now you can modify :code:`main.cpp` and build your project by following the instructions given below.

If you're not a contributor to the project you should comment out the last line in :code:`CMakeLists.txt` ::

   # add_subdirectory(ext/Google_tests)

Otherwise you should clone the GoogleTest repository in order to correctly setup the testing environment::

   cd ext/Google_tests
   mkdir lib
   cd lib
   git clone https://github.com/google/googletest
   cd ../../../

Finally: ::

   cd pointprocess
   mkdir build
   cd build
   cmake ..
   make

Your executable is in the :code:`build/` folder.