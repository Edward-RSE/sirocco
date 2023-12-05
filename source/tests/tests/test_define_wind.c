/** ********************************************************************************************************************
 *
 *  @file test_define_wind.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date November 2023
 *
 *  @brief Unit tests for `define_wind`
 *
 *  TODO: checking actual against test data should be modular, rather than copy and paste
 *
 * ****************************************************************************************************************** */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../python.h"
#include "../assert.h"

char *PYTHON_ENV;
char TEST_CWD[LINELENGTH];
char ATOMIC_DATA_TARGET[LINELENGTH];
char ATOMIC_DATA_DEST[LINELENGTH];

#define PATH_SEPARATOR '/'
#define TEST_DATA_LENGTH 2056

/** *******************************************************************************************************************
 *
 * @brief Free a pointer and set to NULL
 *
 * @param [in]  void** ptr  An address to the pointer to free and set to null
 *
 * @details
 *
 * To use this function, you need to pass the address of a pointer cast as void**, e.g.
 *
 *      free_and_null((void**) &wmain);
 *
 * ****************************************************************************************************************** */

void
free_and_null (void **ptr)
{
  if (ptr != NULL & *ptr != NULL)
  {
    free (*ptr);
    *ptr = NULL;
  }
}

/** *******************************************************************************************************************
 *
 * @brief Get the last component in a file path.
 *
 * @details
 *
 * This will return the file name and extension within a file path.
 *
 * ****************************************************************************************************************** */

const char *
get_last_component (const char *path)
{
  const char *last_separator = strrchr (path, PATH_SEPARATOR);
  return last_separator ? last_separator + 1 : path;
}

/** *******************************************************************************************************************
 *
 * @brief Find the location of the atomic data, relative to a model.
 *
 * @details
 *
 * This function exists to modify the atomic data in the parameter file to make it an absolute path to where
 * it has been set at suite initialisation.
 *
 * ****************************************************************************************************************** */

int
set_atomic_data_filename (void)
{
  char temp_filepath[LINELENGTH];
  char atomic_data_filepath[LINELENGTH];

  if (strlen (ATOMIC_DATA_DEST) + strlen (get_last_component (geo.atomic_filename)) >= LINELENGTH)
  {
    perror ("New atomic data filepath will be too long for internal buffer");
    return EXIT_FAILURE;
  }

  strncpy (temp_filepath, ATOMIC_DATA_DEST, LINELENGTH - 1);
  temp_filepath[LINELENGTH - 1] = '\0';

  size_t len1 = strlen (temp_filepath);
  if (len1 > 0 && temp_filepath[len1 - 1] != PATH_SEPARATOR)
  {
    strncat (temp_filepath, "/", LINELENGTH - len1 - 1);
  }
  snprintf (atomic_data_filepath, LINELENGTH, "%s%s", temp_filepath, get_last_component (geo.atomic_filename));

  if (strlen (atomic_data_filepath) >= LINELENGTH)
  {
    perror ("Buffer overflow when creating new atomic data filepath");
    return EXIT_FAILURE;
  }

  strcpy (geo.atomic_filename, atomic_data_filepath);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after a model has run as a test case.
 *
 * @param [in]  root_name  The root name of the model in $PYTHON/source/tests/test_data/define_wind
 *
 * @details
 *
 * ****************************************************************************************************************** */

int
cleanup_model (const char *root_name)
{
  int n_row;
  int n_plasma;
  char parameter_filepath[LINELENGTH];

  PlasmaPtr plasma_cell;
  MacroPtr macro_cell;

  (void) root_name;

  snprintf (parameter_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/%s.pf", PYTHON_ENV, files.root);
  if (cpar (parameter_filepath) != 1)   /* cpar returns 1 when something is "normal" */
  {
    return EXIT_FAILURE;
  }

  /* free domains */
  free_and_null ((void **) &zdom);

  /* free dynamic grid properties */
  free_and_null ((void **) &wmain);

  /* NPLASMA + 1 is the dummy plasma cell */
  for (n_plasma = 0; n_plasma < NPLASMA + 1; ++n_plasma)
  {
    plasma_cell = &plasmamain[n_plasma];
    free (plasma_cell->density);
    free (plasma_cell->partition);
    free (plasma_cell->ioniz);
    free (plasma_cell->recomb);
    free (plasma_cell->scatters);
    free (plasma_cell->xscatters);
    free (plasma_cell->heat_ion);
    free (plasma_cell->heat_inner_ion);
    free (plasma_cell->cool_rr_ion);
    free (plasma_cell->lum_rr_ion);
    free (plasma_cell->inner_recomb);
    free (plasma_cell->inner_ioniz);
    free (plasma_cell->cool_dr_ion);
    free (plasma_cell->levden);
    free (plasma_cell->recomb_simple);
    free (plasma_cell->recomb_simple_upweight);
    free (plasma_cell->kbf_use);
  }

  free_and_null ((void **) &plasmamain);
  free_and_null ((void **) &photstoremain);
  free_and_null ((void **) &matomphotstoremain);        /* This one doesn't care about if macro atoms are used or not */

  if (nlevels_macro > 0)
  {
    for (n_plasma = 0; n_plasma < NPLASMA + 1; n_plasma++)
    {
      macro_cell = &macromain[n_plasma];
      free (macro_cell->jbar);
      free (macro_cell->jbar_old);
      free (macro_cell->gamma);
      free (macro_cell->gamma_old);
      free (macro_cell->gamma_e);
      free (macro_cell->gamma_e_old);
      free (macro_cell->alpha_st);
      free (macro_cell->alpha_st_old);
      free (macro_cell->alpha_st_e);
      free (macro_cell->alpha_st_e_old);
      free (macro_cell->recomb_sp);
      free (macro_cell->recomb_sp_e);
      free (macro_cell->matom_emiss);
      free (macro_cell->matom_abs);
      free (macro_cell->cooling_bf);
      free (macro_cell->cooling_bf_col);
      free (macro_cell->cooling_bb);

      if (macro_cell->store_matom_matrix == TRUE)
      {
        for (n_row = 0; n_row < nlevels_macro + 1; ++n_row)
        {
          free (macro_cell->matom_matrix[n_row]);
        }
        free (macro_cell->matom_matrix);
      }
    }

    free_and_null ((void **) &macromain);
  }

  NDIM2 = 0;
  NPLASMA = 0;                  /* We probably don't need to do this, but better safe than sorry */

  /* free atomic data elements */
  free_and_null ((void **) &ele);
  free_and_null ((void **) &ion);
  free_and_null ((void **) &xconfig);
  free_and_null ((void **) &line);
  free_and_null ((void **) &auger_macro);


  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Initialise all the necessary parameters to define the wind
 *
 * @param [in]  root_name  The root name of the model in $PYTHON/source/tests/test_data/define_wind
 *
 * @details
 *
 * All the required parameters to define the wind, plasma and macro grid are initialised. There are a number of things
 * which aren't initialised (e.g. binary_basics, get_spec_type, init_observers, etc.) as these are not required to
 * initialise the grids.
 *
 * ****************************************************************************************************************** */

int
initialise_model_for_define_wind (const char *root_name)
{
  int n_dom;
  char rdchoice_answer[LINELENGTH];
  char rdchoice_choices[LINELENGTH];
  char parameter_filepath[LINELENGTH];

  geo.run_type = RUN_TYPE_NEW;

  /* Set up parameter file, that way we can get all the parameters from that
   * instead of defining them manually */
  strcpy (files.root, root_name);
  snprintf (parameter_filepath, LINELENGTH, "%s/source/tests/test_data/define_wind/%s.pf", PYTHON_ENV, files.root);
  if (opar (parameter_filepath) != 2)   /* opar returns 2 when reading for the parameter file */
  {
    fprintf (stderr, "Unable to read from parameter file %s.pf", files.root);
    return EXIT_FAILURE;
  }

  zdom = calloc (MaxDom, sizeof (domain_dummy));        /* We'll allocate MaxDom to follow python */
  if (zdom == NULL)
  {
    fprintf (stderr, "Unable to allocate space for domain structure\n");
    return EXIT_FAILURE;
  }
  init_geo ();

  /* Now when we call the initialisation functions or use rdXXX, the rdchoice_choices
   * for the parameter will come from the parameter file */
  rdint ("Wind.number_of_components", &geo.ndomain);
  if (geo.ndomain > MaxDom)
  {
    fprintf (stderr, "Using more domains (%d) in model than MaxDom (%d)\n", geo.ndomain, MaxDom);
    return EXIT_FAILURE;
  }
  strncpy (rdchoice_answer, "star", LINELENGTH);
  snprintf (rdchoice_choices, LINELENGTH, "%d,%d,%d,%d,%d", SYSTEM_TYPE_STAR, SYSTEM_TYPE_CV, SYSTEM_TYPE_BH, SYSTEM_TYPE_AGN,
            SYSTEM_TYPE_PREVIOUS);
  geo.system_type = rdchoice ("System_type(star,cv,bh,agn,previous)", rdchoice_choices, rdchoice_answer);

  /* These routines don't seem to depend on the atomic data or anything which
   * depends on the atomic data */
  const double star_lum = get_stellar_params ();
  get_bl_and_agn_params (star_lum);
  get_disk_params ();

  /* We have to be a bit creative with the atomic data, to munge the correct
   * filepath with what's in the parameter file */
  rdstr ("Atomic_data", geo.atomic_filename);
  if (set_atomic_data_filename ())
  {
    return EXIT_FAILURE;
  }

  /* We should now be able to initialise everything else which seems to have
   * some dependence on the ionisation settings or atomic data */
  init_ionization ();
  setup_atomic_data (geo.atomic_filename);
  for (n_dom = 0; n_dom < geo.ndomain; ++n_dom)
  {
    get_domain_params (n_dom);
    get_wind_params (n_dom);
  }
  setup_windcone ();
  DFUDGE = setup_dfudge ();

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test an SV wind model for an AGN system with macro atoms.
 *
 * @details
 *
 * This uses the data from $PYTHON/source/tests/test_data/define_wind/agn_macro.pf and
 * $PYTHON/source/tests/test_data/define_wind/agn_macro.grid.txt. The latter was created using the Python script in the
 * $PYTHON/source/tests/test_data/define_wind directory.
 *
 * TODO: we don't check anything to do with macro atoms (we don't output this from windsave2table)
 *
 * ****************************************************************************************************************** */

void
test_sv_agn_macro_wind (void)
{
  int n;
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = initialise_model_for_define_wind ("agn_macro");
  if (init_error)
  {
    CU_FAIL_FATAL ("Unable to initialise AGN Macro model");
  }

  define_wind ();

  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/agn_macro.grid.txt", PYTHON_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    CU_FAIL_FATAL ("Unable to open test data for AGN Macro");
  }

  /* For the CV model, we'll be looking at these parameters */
  int i, j, inwind;
  double x, z, xcen, zcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /* Here's what the header of the file should look like:
     * # i j x z xcen zcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */
    const short n_read = sscanf (test_data_line,
                                 "%d %d %le %le %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &j, &x, &z, &xcen, &zcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4, &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v, &dvds_max, &gamma);
    if (n_read != 29)
    {
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d wmain */
    wind_ij_to_n (0, i, j, &n);
    wind_cell = &wmain[n];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[0], x, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[2], z, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[0], xcen, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[2], zcen, FRACTIONAL_ERROR);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, FRACTIONAL_ERROR);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, FRACTIONAL_ERROR);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, FRACTIONAL_ERROR);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[8] / (n_h * ele[2].abun), c4, FRACTIONAL_ERROR);
  }

  fclose (fp);
  cleanup_model ("agn_macro");
}

/** *******************************************************************************************************************
 *
 * @brief Test a SV wind model for a CV system.
 *
 * @details
 *
 * This uses the data from $PYTHON/source/tests/test_data/define_wind/cv.pf and
 * $PYTHON/source/tests/test_data/define_wind/cv.grid.txt. The latter was created using the Python script in the
 * $PYTHON/source/tests/test_data/define_wind directory.
 *
 * ****************************************************************************************************************** */

void
test_sv_cv_wind (void)
{
  int n;
  FILE *fp;
  char test_data_line[TEST_DATA_LENGTH];
  char test_data_filename[LINELENGTH];

  WindPtr wind_cell;
  PlasmaPtr plasma_cell;

  const int init_error = initialise_model_for_define_wind ("cv");
  if (init_error)
  {
    CU_FAIL_FATAL ("Unable to initialise CV model");
  }

  /* With the defined, we can try and create the wind */
  define_wind ();

  /* And now we can compare our created grid to the "ground truth" grid */
  snprintf (test_data_filename, LINELENGTH, "%s/source/tests/test_data/define_wind/cv.grid.txt", PYTHON_ENV);
  fp = fopen (test_data_filename, "r");
  if (fp == NULL)
  {
    CU_FAIL_FATAL ("Unable to open test data for CV model");
  }

  /* For the CV model, we'll be looking at these parameters */
  int i, j, inwind;
  double x, z, xcen, zcen;
  double v_x, v_y, v_z;
  double vol, rho, ne, h1, c4;
  double t_e, t_r;
  double dv_x_dx, dv_x_dy, dv_x_dz;
  double dv_y_dx, dv_y_dy, dv_y_dz;
  double dv_z_dx, dv_z_dy, dv_z_dz;
  double div_v, dvds_max;
  double gamma;

  /* Skip the first line */
  if (fgets (test_data_line, TEST_DATA_LENGTH, fp) == NULL)
  {
    CU_FAIL_FATAL ("Unable to read first line of test data");
  }

  while (fgets (test_data_line, TEST_DATA_LENGTH, fp) != NULL)
  {
    /* Here's what the header of the file should look like:
     * # i j x z xcen zcen inwind v_x v_y v_z vol rho ne t_e t_r h1 c4 dv_x_dx dv_y_dx dv_z_dx dv_x_dy dv_y_dy dv_z_dy
     * dv_x_dz dv_y_dz dv_z_dz div_v dvds_max gamma
     */
    const short n_read = sscanf (test_data_line,
                                 "%d %d %le %le %le %le %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                                 &i, &j, &x, &z, &xcen, &zcen, &inwind, &v_x, &v_y, &v_z, &vol, &rho, &ne, &t_e, &t_r, &h1, &c4, &dv_x_dx,
                                 &dv_y_dx, &dv_z_dx, &dv_x_dy, &dv_y_dy, &dv_z_dy, &dv_x_dz, &dv_y_dz, &dv_z_dz, &div_v, &dvds_max, &gamma);
    if (n_read != 29)
    {
      CU_FAIL_FATAL ("Test data is in an invalid format");
    }

    /* Convert wind indices into an n in 1d wmain */
    wind_ij_to_n (0, i, j, &n);
    wind_cell = &wmain[n];
    plasma_cell = &plasmamain[wind_cell->nplasma];

    /* cell positions */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[0], x, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->x[2], z, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[0], xcen, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xcen[2], zcen, FRACTIONAL_ERROR);
    CU_ASSERT_EQUAL_FATAL (wind_cell->inwind, inwind);

    /* The default behaviour of Python's output tools (e.g. windsave2table) is
     * to ignore file which are not fully in the wind. So we shall also ignore
     * them here */
    if (wind_cell->inwind != W_ALL_INWIND)
    {
      continue;
    }

    /* velocities */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[0], v_x, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[1], v_y, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v[2], v_z, FRACTIONAL_ERROR);
    /* velocity gradients */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][0], dv_x_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][1], dv_x_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[0][2], dv_x_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][0], dv_y_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][1], dv_y_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[1][2], dv_y_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][0], dv_z_dx, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][1], dv_z_dy, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->v_grad[2][2], dv_z_dz, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->div_v, div_v, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->dvds_max, dvds_max, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (wind_cell->xgamma, gamma, FRACTIONAL_ERROR);

    /* Some things (plasma properties) are stored in plasma cells */
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->rho, rho, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->ne, ne, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_e, t_e, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->t_r, t_r, FRACTIONAL_ERROR);

    /* Ion abundances are tested in their number density relative to Hydrogen.
     * This is the default output option in windsave2table */
    const double n_h = rho2nh * plasma_cell->rho;
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[0] / (n_h * ele[0].abun), h1, FRACTIONAL_ERROR);
    CU_ASSERT_DOUBLE_FRACTIONAL_EQUAL_FATAL (plasma_cell->density[8] / (n_h * ele[2].abun), c4, FRACTIONAL_ERROR);
  }

  fclose (fp);
  cleanup_model ("cv");
}

/** *******************************************************************************************************************
 *
 * @brief Clean up after the `define_wind` tests have run.
 *
 * @details
 *
 * This will deallocate the dynamic memory allocated for the domain structure and will also remove the symbolic link
 * to the atomic data from the test directory.
 *
 * ****************************************************************************************************************** */

int
suite_teardown (void)
{
  if (unlink (ATOMIC_DATA_DEST) != EXIT_SUCCESS)
  {
    perror ("Unable to unlink test data symbolic link");        /* We won't worry too hard about this */
  }

  free (zdom);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Initialise the testing suite for testing `define_wind`
 *
 * @details
 *
 * To be able to create a plasma grid, we need atomic data. The main purpose of this function is to figure out where
 * the tests are being run and to create a symbolic link to the atomic data and to   create file paths to the atomic
 * data. This function will also initialise some global properties in Python, allocating space for the domain structure
 * and, initialise the geometry structure and initialising MPI if that is required. Note that the unit tests should only
 * ever be run using a single MPI rank.
 *
 * Also sets some global variables (to this file) associated with environment variables.
 *
 * ****************************************************************************************************************** */

int
suite_init (void)
{
  struct stat sb;

#ifdef MPI_ON
  MPI_Comm_rank (MPI_COMM_WORLD, &rank_global);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi_global);
#else
  rank_global = 0;
  np_mpi_global = 1;
#endif

  /* Find the PYTHON env var and the directory where tests are called. We need
   * these to create a symbolic link to the atomic data required for the
   * tests */
  PYTHON_ENV = getenv ("PYTHON");
  if (PYTHON_ENV == NULL)
  {
    fprintf (stderr, "Failed to find PYTHON environment variable\n");
    return EXIT_FAILURE;
  }
  if (getcwd (TEST_CWD, LINELENGTH) == NULL)
  {
    perror ("Failed to find current working directory for tests");
    return EXIT_FAILURE;
  }

  /* Set global variables for atomic data */
  snprintf (ATOMIC_DATA_TARGET, LINELENGTH, "%s/xdata", PYTHON_ENV);
  if (!(stat (ATOMIC_DATA_TARGET, &sb) == EXIT_SUCCESS && S_ISDIR (sb.st_mode)))
  {
    perror ("Unable to find atomic data directory");
    return EXIT_FAILURE;
  }
  snprintf (ATOMIC_DATA_DEST, LINELENGTH, "%s/data", TEST_CWD);
  if (symlink (ATOMIC_DATA_TARGET, ATOMIC_DATA_DEST) != EXIT_SUCCESS)
  {
    /* If the symlink exists, we'll try not worry about it as if something is
     * wrong with the atomic data it'll be caught later */
    if (errno != EEXIST)
    {
      perror ("Unable to created symbolic link for atomic data for test case");
      return EXIT_FAILURE;
    }
  }

  /* Now initialise the things Python will need to work, and things which aren't
   * specific to the model being tested such as the domain allocation */
  rel_mode = REL_MODE_FULL;


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

  if ((CU_add_test (suite, "SV: Cataclysmic Variable", test_sv_cv_wind) == NULL) ||
      (CU_add_test (suite, "SV: AGN Macro", test_sv_agn_macro_wind) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Define Wind` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
