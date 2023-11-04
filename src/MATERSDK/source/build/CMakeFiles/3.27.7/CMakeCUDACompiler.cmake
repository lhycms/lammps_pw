set(CMAKE_CUDA_COMPILER "/share/app/cuda/cuda-11.8/bin/nvcc")
set(CMAKE_CUDA_HOST_COMPILER "")
set(CMAKE_CUDA_HOST_LINK_LAUNCHER "/opt/rh/devtoolset-7/root/usr/bin/g++")
set(CMAKE_CUDA_COMPILER_ID "NVIDIA")
set(CMAKE_CUDA_COMPILER_VERSION "11.8.89")
set(CMAKE_CUDA_DEVICE_LINKER "/share/app/cuda/cuda-11.8/bin/nvlink")
set(CMAKE_CUDA_FATBINARY "/share/app/cuda/cuda-11.8/bin/fatbinary")
set(CMAKE_CUDA_STANDARD_COMPUTED_DEFAULT "14")
set(CMAKE_CUDA_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_CUDA_COMPILE_FEATURES "cuda_std_03;cuda_std_11;cuda_std_14;cuda_std_17")
set(CMAKE_CUDA03_COMPILE_FEATURES "cuda_std_03")
set(CMAKE_CUDA11_COMPILE_FEATURES "cuda_std_11")
set(CMAKE_CUDA14_COMPILE_FEATURES "cuda_std_14")
set(CMAKE_CUDA17_COMPILE_FEATURES "cuda_std_17")
set(CMAKE_CUDA20_COMPILE_FEATURES "")
set(CMAKE_CUDA23_COMPILE_FEATURES "")

set(CMAKE_CUDA_PLATFORM_ID "Linux")
set(CMAKE_CUDA_SIMULATE_ID "GNU")
set(CMAKE_CUDA_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_CUDA_SIMULATE_VERSION "7.3")



set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CUDA_COMPILER_LOADED 1)
set(CMAKE_CUDA_COMPILER_ID_RUN 1)
set(CMAKE_CUDA_SOURCE_FILE_EXTENSIONS cu)
set(CMAKE_CUDA_LINKER_PREFERENCE 15)
set(CMAKE_CUDA_LINKER_PREFERENCE_PROPAGATES 1)
set(CMAKE_CUDA_LINKER_DEPFILE_SUPPORTED )

set(CMAKE_CUDA_SIZEOF_DATA_PTR "8")
set(CMAKE_CUDA_COMPILER_ABI "ELF")
set(CMAKE_CUDA_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CUDA_LIBRARY_ARCHITECTURE "")

if(CMAKE_CUDA_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CUDA_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CUDA_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CUDA_COMPILER_ABI}")
endif()

if(CMAKE_CUDA_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CUDA_COMPILER_TOOLKIT_ROOT "/share/app/cuda/cuda-11.8")
set(CMAKE_CUDA_COMPILER_TOOLKIT_LIBRARY_ROOT "/share/app/cuda/cuda-11.8")
set(CMAKE_CUDA_COMPILER_TOOLKIT_VERSION "11.8.89")
set(CMAKE_CUDA_COMPILER_LIBRARY_ROOT "/share/app/cuda/cuda-11.8")

set(CMAKE_CUDA_ARCHITECTURES_ALL "35-real;37-real;50-real;52-real;53-real;60-real;61-real;62-real;70-real;72-real;75-real;80-real;86-real;87-real;89-real;90")
set(CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR "35-real;50-real;60-real;70-real;80-real;90")
set(CMAKE_CUDA_ARCHITECTURES_NATIVE "86-real")

set(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES "/share/app/cuda/cuda-11.8/targets/x86_64-linux/include")

set(CMAKE_CUDA_HOST_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES "/share/app/cuda/cuda-11.8/targets/x86_64-linux/lib/stubs;/share/app/cuda/cuda-11.8/targets/x86_64-linux/lib")
set(CMAKE_CUDA_HOST_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

set(CMAKE_CUDA_IMPLICIT_INCLUDE_DIRECTORIES "/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/mkl/include;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/pstl/include;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/pstl/stdlib;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/tbb/include;/opt/rh/devtoolset-7/root/usr/include/c++/7;/opt/rh/devtoolset-7/root/usr/include/c++/7/x86_64-redhat-linux;/opt/rh/devtoolset-7/root/usr/include/c++/7/backward;/opt/rh/devtoolset-7/root/usr/lib/gcc/x86_64-redhat-linux/7/include;/usr/local/include;/opt/rh/devtoolset-7/root/usr/include;/usr/include")
set(CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES "stdc++;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES "/share/app/cuda/cuda-11.8/targets/x86_64-linux/lib/stubs;/share/app/cuda/cuda-11.8/targets/x86_64-linux/lib;/opt/rh/devtoolset-7/root/usr/lib/gcc/x86_64-redhat-linux/7;/opt/rh/devtoolset-7/root/usr/lib64;/lib64;/usr/lib64;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin;/share/app/intel2020u4/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8;/opt/rh/devtoolset-7/root/usr/lib")
set(CMAKE_CUDA_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

set(CMAKE_CUDA_RUNTIME_LIBRARY_DEFAULT "STATIC")

set(CMAKE_LINKER "/opt/rh/devtoolset-7/root/usr/bin/ld")
set(CMAKE_AR "/opt/rh/devtoolset-7/root/usr/bin/ar")
set(CMAKE_MT "")
