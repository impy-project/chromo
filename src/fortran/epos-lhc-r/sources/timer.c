/* timer.c */
#include <sys/time.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>

void timer_(double *etime)
{
    struct timeval tv;
    struct tms tu;
    gettimeofday(&tv, NULL);
    times(&tu);
    etime[0] = (double)tu.tms_utime;
    etime[1] = (double)tu.tms_stime;
    etime[2] = (double)tv.tv_sec;
    etime[3] = (double)tv.tv_usec / 1000.0;
    etime[4] = 0.0;
}
