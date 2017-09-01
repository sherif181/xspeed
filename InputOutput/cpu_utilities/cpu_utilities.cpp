/*
 * cpu_utilities.cpp
 *
 *  Created on: 24-Apr-2015
 *      Author: amit
 */

#include "cpu_utilities.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/times.h"
#include "sys/vtimes.h"

static clock_t lastCPU, lastSysCPU, lastUserCPU;
static int numProcessors;

void init_cpu_usage() {
	FILE* file;
	struct tms timeSample;
	char line[128];

	lastCPU = times(&timeSample);
	lastSysCPU = timeSample.tms_stime;
	lastUserCPU = timeSample.tms_utime;

	file = fopen("/proc/cpuinfo", "r");
	numProcessors = 0;
	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "processor", 9) == 0)
			numProcessors++;
	}
	fclose(file);
}

double getCurrent_ProcessCPU_usage() {
	struct tms timeSample;
	clock_t now;
	double percent;

	now = times(&timeSample);
	if (now <= lastCPU || timeSample.tms_stime < lastSysCPU
			|| timeSample.tms_utime < lastUserCPU) {
		//Overflow detection. Just skip this value.
		percent = -1.0;
	} else {
		percent = (timeSample.tms_stime - lastSysCPU)
				+ (timeSample.tms_utime - lastUserCPU);
		percent /= (now - lastCPU);
		percent /= numProcessors;
		percent *= 100;
	}
	lastCPU = now;
	lastSysCPU = timeSample.tms_stime;
	lastUserCPU = timeSample.tms_utime;

	return percent;
}
