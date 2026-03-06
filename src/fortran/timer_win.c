/* timer_win.c - Windows-compatible replacement for epos-lhc-r/sources/timer.c */
#include <windows.h>

void timer_(double *etime)
{
    FILETIME creation, exit, kernel, user;
    LARGE_INTEGER freq, now;

    /* Wall clock time */
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&now);
    double wall_sec = (double)now.QuadPart / (double)freq.QuadPart;

    /* Process CPU times */
    double user_time = 0.0, sys_time = 0.0;
    if (GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user)) {
        ULARGE_INTEGER u, k;
        u.LowPart = user.dwLowDateTime;
        u.HighPart = user.dwHighDateTime;
        k.LowPart = kernel.dwLowDateTime;
        k.HighPart = kernel.dwHighDateTime;
        /* Convert from 100-ns units to clock ticks (approximate) */
        user_time = (double)u.QuadPart / 1e7;
        sys_time  = (double)k.QuadPart / 1e7;
    }

    etime[0] = user_time;
    etime[1] = sys_time;
    etime[2] = (double)(long long)wall_sec;
    etime[3] = (wall_sec - (long long)wall_sec) * 1000.0;
    etime[4] = 0.0;
}
