add_test([=[AngularPartTest.get_level]=]  /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/bin/matersdk/feature/mtp/test_angular_part [==[--gtest_filter=AngularPartTest.get_level]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[AngularPartTest.get_level]=]  PROPERTIES WORKING_DIRECTORY /data/home/liuhanyu/hyliu/code1/lammps_pw/src/MATERSDK/source/build/matersdk/feature/mtp SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  test_angular_part_TESTS AngularPartTest.get_level)
