// Author: Hang-Hyun Jo
// Description: functions for some basic operations

// for event sequence analysis
long read_timings(char *folder, char *filename, long *timings, long n_max){
    long i, timing, timing0; 
    char output[500]; 
    FILE *file_in; 

    for(i = 0; i < n_max; i ++) timings[i] = 0;

    sprintf(output, "%stimings_%s.txt", folder, filename);
    file_in = fopen(output, "r"); 
    i = 0;
    timing0 = -1;
    while(fscanf(file_in, "%ld", &timing) && !feof(file_in)){
        if(timing == timing0) continue;
        timings[i] = timing;
        i ++;
        if(i == n_max) printf("data size exceeding the limit in read_timings()\n");
        timing0 = timing;
    }
    fclose(file_in); 

    return i;
}

// exchange one form to another
long get_iet_sequence_from_timings(long *timings, long n, long *iet_sequence){
    long i, timing0, timing, iet, iet_max;

    for(i = 0; i < n - 1; i ++) iet_sequence[i] = 0;

    iet_max = -1;
    timing0 = -1;
    for(i = 0; i < n; i ++){
        timing = timings[i];
        if(timing0 >= 0){
            iet = timing - timing0;
            iet_sequence[i - 1] = iet;
            if(iet > iet_max) iet_max = iet;
        }
        timing0 = timing;
    }

    return iet_max;
}



