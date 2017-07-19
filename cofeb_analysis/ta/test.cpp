/**
* @Author: Jenkins Alec <alec>
* @Date:   2017-07-09T13:25:23-07:00
* @Project: LTSPM analysis
* @Last modified by:   alec
* @Last modified time: 2017-07-09T13:31:11-07:00
*/

#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" void fun()
#else
void fun(void)
#endif
{
    printf( "%i", 12 );
}
