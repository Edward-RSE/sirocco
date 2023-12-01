/** ********************************************************************************************************************
 *
 * @file assert.h
 * @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 * @date August 2023
 *
 * @brief Additional custom assertions for testing
 *
 * ****************************************************************************************************************** */

#ifndef ASSERT_H
#define ASSERT_H

#define TRUE 1
#define FALSE 0
#define CU_BUFFER_SIZE 1024
#define EPSILON 1.0e-6
#define FRACTION 1e-3

#define CU_ASSERT_DOUBLE_EQUAL_FRAC_FATAL(actual, expected, granularity)                                               \
{                                                                                                                      \
  if ((double) actual < EPSILON && (double) expected < EPSILON)                                                        \
  {                                                                                                                    \
    ;                                                                                                                  \
  }                                                                                                                    \
  else                                                                                                                 \
  {                                                                                                                    \
    CU_assertImplementation((((double) actual - (double) expected) / (double) expected) <= granularity, __LINE__,      \
    "CU_ASSERT_DOUBLE_EQUAL_FRAC_FATAL", __FILE__, "", TRUE);                                                          \
  }                                                                                                                    \
}

#define CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL(actual, expected, size, granularity)                                            \
{                                                                                                                      \
    int i;                                                                                                             \
    int not_equal_count = 0;                                                                                           \
    for (i = 0; i < size; ++i)                                                                                         \
    {                                                                                                                  \
        if (fabs((double)(actual[i]) - (double)(expected[i])) >= fabs((double)(granularity)))                          \
        {                                                                                                              \
            not_equal_count += 1;                                                                                      \
        }                                                                                                              \
    }                                                                                                                  \
  CU_assertImplementation(not_equal_count == 0, __LINE__, "CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL", __FILE__, "", TRUE);       \
}

#define CU_ASSERT_DOUBLE_ARRAY_EQUAL(actual, expected, size, granularity)                                                  \
{                                                                                                                      \
    int i;                                                                                                             \
    int not_equal_count = 0;                                                                                           \
    for (i = 0; i < size; ++i)                                                                                         \
    {                                                                                                                  \
        if (fabs((double)(actual[i]) - (double)(expected[i])) >= fabs((double)(granularity)))                          \
        {                                                                                                              \
            not_equal_count += 1;                                                                                      \
        }                                                                                                              \
    }                                                                                                                  \
  CU_assertImplementation(not_equal_count == 0, __LINE__, "CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL", __FILE__, "", FALSE);      \
}

#define CU_FAIL_MSG_FATAL(...)                                                                                         \
{                                                                                                                      \
    char buffer[CU_BUFFER_SIZE];                                                                                       \
    snprintf(buffer, CU_BUFFER_SIZE, __VA_ARGS__);                                                                     \
    CU_FAIL_FATAL(buffer);                                                                                             \
}

#endif
