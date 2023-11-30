/** ********************************************************************************************************************
 *
 *  @file main.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief Main function for unit tests
 *
 *  https://gitlab.com/cunity/cunit
 *
 * ****************************************************************************************************************** */

#include <time.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

/* Test suite prototypes */
#include "tests/tests.h"

/* Python logging prototypes */
#include "../log.h"

/* Including templates.h has so undesired effects, so we'll define the prototypes
   as we need them */
int init_rand (int seed);

/** *******************************************************************************************************************
 *
 *  @brief Entry point for unit tests
 *
 * ****************************************************************************************************************** */

int
main (void)
{
  /* I don't want any output from Python, thank you very much */
  Log_set_verbosity (0);

  /* Initialize the CUnit test registry */
  if (CU_initialize_registry () != CUE_SUCCESS)
  {
    return CU_get_error ();
  }

  /* Initialise some stuff for Python -- e.g. RNG */
  init_rand ((int) time (NULL));

  /* Create test suites */
  create_compton_test_suite ();
  create_define_wind_test_suite ();
  create_matrix_test_suite ();

  /* Run the test suites */
  CU_basic_set_mode (CU_BRM_VERBOSE);
  CU_basic_run_tests ();

  const int num_tests_failed = CU_get_number_of_tests_failed ();

  if (num_tests_failed > 0)
  {
    printf ("\033[1;31m%d test(s) failed\n\033[1;0m", num_tests_failed);
  }
  else
  {
    printf ("\033[1;32mAll tests ran successfully\n\033[1;0m");
  }

  /* Clean up the CUnit registry */
  CU_cleanup_registry ();

  return num_tests_failed;
}
