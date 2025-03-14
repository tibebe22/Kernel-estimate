// Author: Hang-Hyun Jo
// Description: main functions

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

void get_ready(char *folder, char *filename, long *timings, long *iet_sequence, long **tree, long num_event, long isTreeRead){
    long num_iet;

    num_iet = num_event - 1;

    // get interevent time sequence
    get_iet_sequence_from_timings(timings, num_event, iet_sequence);

    // get tree
    if(isTreeRead == 0){ // detect tree after obtaining inter-event time sequence
        get_tree(iet_sequence, num_iet, tree);
        print_tree(folder, filename, tree, num_iet);
        printf("get and print tree done\n");
    }
    else{ // read tree from the file
        read_tree(folder, filename, tree, num_iet);
        printf("read tree done\n");
    }
}

void burst_analysis_option(char *folder, char *filename, long ens, long *timings, long *iet_sequence, long **tree, long num_event, char *options){
    long i, j, k, num_iet, train_Dts[6], num_Dt;
    double autocorrel_params[3], logbin_params[3];
    char filename_ens[500];

    num_iet = num_event - 1;
    if(ens >= 0) sprintf(filename_ens, "%s_ens%ld", filename, ens);
    else sprintf(filename_ens, "%s", filename);

    // kernel 2D
    if(options[1] == '1'){
        printf("get kernel2D\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_kernel2DLog(folder, filename_ens, tree, num_iet, logbin_params);
    }
}
