add_test([=[AlignedArrayTest.Size]=]  /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/bin/core/test_AlignedArray [==[--gtest_filter=AlignedArrayTest.Size]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[AlignedArrayTest.Size]=]  PROPERTIES WORKING_DIRECTORY /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/core SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test([=[AlignedArrayTest.Resize]=]  /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/bin/core/test_AlignedArray [==[--gtest_filter=AlignedArrayTest.Resize]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[AlignedArrayTest.Resize]=]  PROPERTIES WORKING_DIRECTORY /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/core SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test([=[AlignedArrayTest.GetElementByIndex]=]  /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/bin/core/test_AlignedArray [==[--gtest_filter=AlignedArrayTest.GetElementByIndex]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[AlignedArrayTest.GetElementByIndex]=]  PROPERTIES WORKING_DIRECTORY /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/core SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  test_AlignedArray_TESTS AlignedArrayTest.Size AlignedArrayTest.Resize AlignedArrayTest.GetElementByIndex)
