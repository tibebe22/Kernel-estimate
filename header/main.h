// Author: Hang-Hyun Jo
// Description: main functions

void get_ready(char *folder, char *filename, long *timings, long *iet_sequence, long **tree, long num_event, long isTreeRead){
    long num_iet;

    /* read timings
    timings = vector_long(0, num_event_max);
    num_event = read_timings(timings, num_event_max, filename, folder);
    */

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

