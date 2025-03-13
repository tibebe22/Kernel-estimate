// Author: Hang-Hyun Jo
// Description: functions for some basic operations

// get distribution
void get_distr_long(long *sequence, long num_item, long *distr, long x_min, long x_max){
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = sequence[i];
        distr[x] ++;
    }
}

void get_distr_double(double *sequence, long num_item, double resolution, long *distr, long x_min, long x_max){ // x is in [x_min, x_max]
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = (long)(sequence[i] / resolution + 0.5);
        if(x < x_min || x > x_max) printf("out of range\n");
        distr[x] ++;
    }
}

long find_max(long *sequence, long num_item){
    long i, x, x_max;

    x_max = -1;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] > x_max)
            x_max = sequence[i];
    }

    return x_max;
}

long find_min(long *sequence, long num_item){
    long i, x, x_min;

    x_min = 1e8;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] < x_min)
            x_min = sequence[i];
    }

    return x_min;
}

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

void get_timing_sequence_from_iet_sequence(long *iet_sequence, long num_iet, long *timing_sequence){
    long i, timing;

    for(i = 0; i <= num_iet; i ++) timing_sequence[i] = 0;

    timing = 0;
    for(i = 0; i < num_iet; i ++){
        timing += iet_sequence[i];
        timing_sequence[i + 1] = timing;
    }
}

void get_iet_sequence_from_tree_direct(long **tree, long *iet_sequence, long num_iet){
    long i, loc;

    for(i = 0; i < num_iet; i ++){
        loc = tree[i][0];
        iet_sequence[loc] = tree[i][1];
    }
}

void get_iet_sequence_from_tree(long **tree, long root_rank, long root_loc, long *iet_sequence, long num_iet){
    long rank_left, rank_right, train0, train1, loc;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(rank_left >= 0){
        train1 = tree[rank_left][5]; // train size of right branch of "rank_left"
        loc = root_loc - train1;
        iet_sequence[loc] = tree[rank_left][1];
        get_iet_sequence_from_tree(tree, rank_left, loc, iet_sequence, num_iet);
    }
    if(rank_right >= 0){
        train0 = tree[rank_right][4]; // train size of left branch of "rank_right"
        loc = root_loc + train0;
        iet_sequence[loc] = tree[rank_right][1];
        get_iet_sequence_from_tree(tree, rank_right, loc, iet_sequence, num_iet);
    }
}






