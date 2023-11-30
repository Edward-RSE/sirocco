/** ********************************************************************************************************************
 *
 *  @file test_define_wind.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date November 2023
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../python.h"
#include "../assert.h"

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *
 * ****************************************************************************************************************** */

void
test_define_wind (void)
{

}

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *
 * ****************************************************************************************************************** */

int
suite_teardown (void)
{

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *
 * ****************************************************************************************************************** */

int
suite_init (void)
{

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for define wind
 *
 * @details
 *
 * This function will create a test suite for functions related to creating the wind and plasma grid
 *
 * ****************************************************************************************************************** */

void
create_define_wind_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Define Wind", suite_init, suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Define Wind` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "Test Define Wind", test_define_wind) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Define Wind` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
