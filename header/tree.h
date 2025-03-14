// Authors: Hang-Hyun Jo and Tibebe Birhanu
// Description: functions for detecting and processing tree structure

/************ DETECTING TREE *****************/

void search_branch(long **tree, long root_rank, long list_size){
    long root_loc, bound_left, bound_right, k, k_loc;

    if(root_rank == list_size - 1)
        return;

    root_loc = tree[root_rank][0];
    bound_left = tree[root_rank][4];
    bound_right = tree[root_rank][5];

    // find the "left" branch
    if(root_loc > 0){
        for(k = root_rank + 1; k < list_size; k ++){
            k_loc = tree[k][0];
            if(k_loc < root_loc && k_loc >= bound_left){
                tree[root_rank][2] = k; // rank of left branch
                tree[k][4] = bound_left; // inherit bound_left from the root
                tree[k][5] = root_loc - 1; // new bound_right
                break;
            }
        }
    }

    // find the "right" branch
    if(root_loc < list_size - 1){
        for(k = root_rank + 1; k < list_size; k ++){
            k_loc = tree[k][0];
            if(k_loc > root_loc && k_loc <= bound_right){
                tree[root_rank][3] = k; // rank of right branch
                tree[k][4] = root_loc + 1; // new bound_left
                tree[k][5] = bound_right; // inherit bound_right from the root
                break;
            }
        }
    }
}

void merge_trains(long **tree, long root_rank){
    long rank_left, rank_right, train;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left == -1){
        if(rank_right == -1){ // no branches
            tree[root_rank][4] = 1;
            tree[root_rank][5] = 1;
        }
        else{ // only right branch
            tree[root_rank][4] = 1;
            tree[root_rank][5] = tree[rank_right][4] + tree[rank_right][5];
        }
    }
    else{
        if(rank_right == -1){  // only left branch
            tree[root_rank][4] = tree[rank_left][4] + tree[rank_left][5];
            tree[root_rank][5] = 1;
        }
        else{ // both branches
            tree[root_rank][4] = tree[rank_left][4] + tree[rank_left][5];
            tree[root_rank][5] = tree[rank_right][4] + tree[rank_right][5];
        }
    }
}

void get_tree(long *iet_sequence, long num_iet, long **tree){
    long i, iet, loc;
    double **iet_sequence_sort;

    iet_sequence_sort = matrix_double(0, num_iet - 1, 0, 1);

    for(i = 0; i < num_iet; i ++){ // rank of iet
        iet_sequence_sort[i][0] = i;
        iet_sequence_sort[i][1] = iet_sequence[i] + genrand_real2() * 0.1;
        // <-- avoid the degeneracy
        tree[i][0] = -1; // loc (=original index)
        tree[i][1] = -1; // iet
        tree[i][2] = -1; // rank of left branch
        tree[i][3] = -1; // rank of right branch
        tree[i][4] = 0; // left bound
        tree[i][5] = num_iet - 1; // right bound
    }

    QuickSort3_descend_double(iet_sequence_sort, 1, 1, 0, num_iet - 1);

    for(i = 0; i < num_iet; i ++){
        tree[i][0] = (long)iet_sequence_sort[i][0]; // loc
        tree[i][1] = (long)iet_sequence_sort[i][1]; // iet
    }

    for(i = 0; i < num_iet; i ++){
        search_branch(tree, i, num_iet);
        if(i % 10000 == 0) printf("searching %ldth node\n", i);
    }
    printf("search done\n");

    for(i = 0; i < num_iet; i ++){ // rank
        tree[i][4] = 0; // train size of left branch of node i
        tree[i][5] = 0; // train size of right branch of node i
    }

    for(i = num_iet - 1; i >= 0; i --) merge_trains(tree, i);
    printf("merge done\n");

    free_matrix_double(iet_sequence_sort, 0, num_iet - 1, 0, 1);
}

void print_tree(char *folder, char *filename, long **tree, long num_iet){
    long i, j;
    char output[500];
    FILE *tree_out;

    sprintf(output, "%stree_%s.txt", folder, filename);
    tree_out = fopen(output, "w");

    for(i = 0; i < num_iet; i ++){
        fprintf(tree_out, "%ld", i);
        for(j = 0; j < 6; j ++){
            fprintf(tree_out, " %ld", tree[i][j]);
        }
        fprintf(tree_out, "\n");
    }
    fclose(tree_out);
}

void read_tree(char *folder, char *filename, long **tree, long num_iet){
    long i, j, x;
    char output[500], temp[500];
    FILE *tree_in;

    sprintf(output, "%stree_%s.txt", folder, filename);
    tree_in = fopen(output, "r");

    for(i = 0; i < num_iet; i ++) for(j = 0; j < 6; j ++) tree[i][j] = 0;

    for(i = 0; i < num_iet; i ++){
        fscanf(tree_in , "%s", temp);
        for(j = 0; j < 6; j ++){
            fscanf(tree_in , "%ld", &x);
            tree[i][j] = x;
        }
    }
    fclose(tree_in);
}

/************ QUANTIFYING TREE: kernel  *********/

//get the distribution of Ms
void get_Ms(long **tree, double **Ms, long num_iet, long num_event_temp){
    long i, j, k, train1, train2;
    
    for(i = 1; i <= num_event_temp; i++){
        for(j = 1; j <= num_event_temp; j++){
            Ms[i][j] = 0;
        }
    }
     
    for(i = 0; i< num_iet; i++){
        train1 = tree[i][4];
        train2 = tree[i][5];
        if(train1 <= num_event_temp && train2 <= num_event_temp)
            Ms[train1][train2]++;
    }
}

// get train distribution for entire time step
void get_train_distribution(long **tree, double **Q_all, long *train_distr, long num_event, long num_event_temp){
    long i, j, k, train1, train2, train, train_min, num_iet;
    double n_1;
    num_iet = num_event - 1;
    for(i = 0; i < num_iet; i++){
        for(j = 1; j <= num_event_temp; j++){
            Q_all[i][j] = 0;
        }
    }
     
    for(i = 1; i < num_event; i ++) train_distr[i] = 0;
    train_distr[num_event] = 1;
    train_min = num_event;
        
    for(i = 0; i < num_iet; i++){
        train1 = tree[i][4];
        train2 = tree[i][5];
        train = train1 + train2;

        if(train1 < train_min) train_min = train1;
        if(train2 < train_min) train_min = train2;

        train_distr[train] --;
        train_distr[train1] ++;
        train_distr[train2] ++;
        n_1 = 1. /(double)(i + 2);
        for(j = train_min; j <= num_event_temp; j++){
            Q_all[i][j] = n_1 * train_distr[j];      
        }
    }   
}

//find the sum of the denumenator part of the kernel equation
double get_sum_denominator(double **Q_all, double **prev_kernel, long train_min, long num_event_temp, long s){
    long i, j, k;
    double summ;
       
    summ = 0.;
    for(j = train_min; j <= num_event_temp; j++){
        for(k = train_min; k <= num_event_temp; k++){
            if(Q_all[s][j] && Q_all[s][k] && prev_kernel[j][k])
                summ += Q_all[s][j] * Q_all[s][k] * prev_kernel[j][k];
        }
    }
    return summ;
}

//the log likelihood 
double get_likelihood(long **tree, double **Q_all, double **P, long num_iet, long num_event_temp){
    long i, j, k, train_min, train1, train2, num_event;
    double sum_temp, summ1, summ2, likeli;
       
    summ1 = 0.;
    summ2 = 0.;
    sum_temp = 0.;
    num_event = num_iet + 1;
    train_min = num_event;
    for(i = 0; i < num_iet; i++){
        train1 = tree[i][4];
        train2 = tree[i][5];
        if(train1 < train_min)train_min = train1;
        if(train2 < train_min)train_min = train2;
        for(j = train_min; j <= num_event_temp; j++){
            for(k = train_min; k <= num_event_temp; k++){
                if(Q_all[i][j] && Q_all[i][k] && P[j][k]) sum_temp += Q_all[i][j] * Q_all[i][k] * P[j][k];  
            } 
        }      
        if(sum_temp) summ2 += log(sum_temp);
        if(train1 <= num_event_temp && train2 <= num_event_temp && P[train1][train2])
            summ1 += log(P[train1][train2]);
    }
    likeli = summ1 - summ2;
    return likeli;
}

//update the kernels
void update_kernel(double **kernel2D, double **prev_kernel, long num_event_temp){
    long i, j, k;
    for(i = 1; i <= num_event_temp; i++){
        for(j = 1; j <= num_event_temp; j++){
            prev_kernel[i][j] = kernel2D[i][j];
        }
    }
}

// estimate K(b,b') for conserving the correlations between b and b'
void get_kernel2D(long **tree, long num_iet, long num_event_temp, double **kernel2D){
    long i, j, k, num_event, train_min, train, train1, train2;
    long *train_distr;
    double n_1, **QQ, **Ms, **Q_all, **prev_kernel, sum_denom;
    double lcurrent1, lcurrent2, lprev1, lprev2, curr_loglikeli, prev_loglikeli, diff, epsilon;
    
    epsilon = 1e-4;
    diff = 1.;
    num_event = num_iet + 1;
    train_distr = vector_long(1, num_event);
    if (num_event_temp > 0 && num_event_temp <= SIZE_MAX / sizeof(double)) {
        QQ = matrix_double(1, num_event_temp, 1, num_event_temp);
        Ms = matrix_double(1, num_event_temp, 1, num_event_temp);
        Q_all = matrix_double(0, num_iet - 1, 1, num_event_temp);
        prev_kernel = matrix_double(1, num_event_temp, 1, num_event_temp); 
    } 
    else {
        printf("Invalid size for allocation, potential overflow detected.\n");
    }
     
    for(i = 1; i <= num_event_temp; i ++){ // i = 0 is for Q(b), otherwise kernel2D[b][b']
        for(j = 1; j <= num_event_temp; j ++){
            prev_kernel[i][j] = 1;
        }
    }
    //get the distribution of Ms
    get_Ms(tree, Ms, num_iet, num_event_temp);
    
    //get train distribution for entire time step
    get_train_distribution(tree, Q_all, train_distr, num_event, num_event_temp);
    
    //printf("kernel 1\n");
    while(diff > epsilon){ 
        train_min = num_event;
        for(i = 1; i <= num_event_temp; i ++){ // i = 0 is for Q(b), otherwise kernel2D[b][b']
            for(j = 1; j <= num_event_temp; j ++){
                kernel2D[i][j] = 0;
                QQ[i][j] = 0;
            }
        } 
         
        for(i = 0; i < num_iet; i ++){
            train1 = tree[i][4];
            train2 = tree[i][5];
            if(train1 < train_min) train_min = train1;
            if(train2 < train_min) train_min = train2;
            sum_denom = get_sum_denominator(Q_all, prev_kernel, train_min, num_event_temp, i); 
            for(j = train_min; j <= num_event_temp; j ++){
                for(k = train_min; k <= num_event_temp; k ++){
                    if(Q_all[i][j] && Q_all[i][k] && sum_denom){
                        QQ[j][k] += (Q_all[i][j] * Q_all[i][k]) / sum_denom;                    
                    }
                }
            }
        }
    
        for(i = train_min; i <= num_event_temp; i ++){
            for(j = train_min; j <= num_event_temp; j ++){
                if(Ms[i][j] && QQ[i][j])
                    kernel2D[i][j] = Ms[i][j] / QQ[i][j];
            }
        }
        
        //normalize the kernel
        for(i = train_min; i <= num_event_temp; i ++){
            for(j = train_min; j <= num_event_temp; j ++){
                if(kernel2D[1][1])
                    kernel2D[i][j] = kernel2D[i][j]/kernel2D[1][1];
            }
        }
    
        prev_loglikeli = get_likelihood(tree, Q_all, prev_kernel, num_iet, num_event_temp);
        printf("the previous loglikelihood %lf\n", prev_loglikeli);
        curr_loglikeli = get_likelihood(tree, Q_all, kernel2D, num_iet, num_event_temp);
        printf("the current loglikelihood %lf\n", curr_loglikeli);
    
        //the condtion is
        diff = fabs(curr_loglikeli - prev_loglikeli) / ( fabs(prev_loglikeli) + 1 );
        // printf("the condtion value %lf\n", diff); 
    
        //update  the kernels
        update_kernel(kernel2D, prev_kernel, num_event_temp); 
    }
    free_vector_long(train_distr, 1, num_event);
    free_matrix_double(QQ, 1, num_event_temp, 1, num_event_temp);
    free_matrix_double(Ms, 1, num_event_temp, 1, num_event_temp);
    free_matrix_double(Q_all, 0, num_iet - 1, 1, num_event_temp); 
    free_matrix_double(prev_kernel, 1, num_event_temp, 1, num_event_temp);
}

// for binning
long get_mid_bin(double *mid_bin, double binStart, double binSize, long num_event){
    long train2, index2;
    double x_factor, x_factor_2, y_end; 
    
    x_factor = exp(binSize);
    x_factor_2 = exp(-binSize * 0.5);

    index2 = 0; train2 = 1; // initialize
    while(index2 < binStart){ // no binning
        mid_bin[index2] = train2;
        index2 ++;
        train2 ++;
    } 
    y_end = train2 * x_factor; // bin range starting from train2
    while(train2 <= num_event){ // log binning
        if(train2 < y_end){
            train2 ++;
        }
        else{
            mid_bin[index2] = y_end * x_factor_2;
            y_end *= x_factor;
            index2 ++;
        }
    }

    return index2; // number of bins
}
        
// repeating module for a given index1 and train1
void kernel2D_yLog(double **kernel2DLog, long **nums, double **kernel2D, long num_event, double binStart, double x_factor, long index1, long train1){
    long train2, index2;
    double y_end; 
    
    // initialize
    index2 = 0; train2 = 1;

    // no binning
    while(index2 < binStart){ 
        kernel2DLog[index1][index2] += kernel2D[train1][train2];
        nums[index1][index2] ++;
        index2 ++;
        train2 ++;
    } 
    
    // bin range starting from train2
    y_end = train2 * x_factor;

    // log binning
    while(train2 <= num_event){
        if(train2 < y_end){
            kernel2DLog[index1][index2] += kernel2D[train1][train2];
            nums[index1][index2] ++;
            train2 ++;
        }
        else{
            y_end *= x_factor;
            index2 ++;
        }
    }
}
        
// get the log-binning of K_{b,b'}
void get_kernel2DLog(char *folder, char *filename, long **tree, long num_iet, double *logbin_params){
    long i, j, k, num_event, train_min, train, train1, train2, index1, index2, num_event_temp;
    long **nums, num_bin, num;
    double **kernel2D, **kernel2DLog, *mid_bin, binStart, binSize, x_factor, x_end, avg;
    char output[500], prefix[500];
    FILE *file_out;

    num_event = num_iet + 1;
    num_event_temp = 2e3;
    //num_event_temp = (long)(num_event / 2.);
    kernel2D = matrix_double(0, num_event_temp, 1, num_event_temp);
    printf("start get kernel 2d\n");
    get_kernel2D(tree, num_iet, num_event_temp, kernel2D);
    printf("kernel 2d done\n");

    binStart = logbin_params[0];
    binSize = logbin_params[1];
    x_factor = exp(binSize);

    // initialize bins' boundaries and mid values
    mid_bin = vector_double(0, num_event_temp);
    for(i = 0; i <= num_event_temp; i ++) mid_bin[i] = 0;
    num_bin = get_mid_bin(mid_bin, binStart, binSize, num_event_temp);

    // declare arrays for log binning using "num_bin"
    kernel2DLog = matrix_double(0, num_bin, 0, num_bin);
    nums = matrix_long(0, num_bin, 0, num_bin);
    for(i = 0; i <= num_bin; i ++){
        for(j = 0; j <= num_bin; j ++){
            kernel2DLog[i][j] = 0;
            nums[i][j] = 0;
        }
    }

    // initialize
    index1 = 0; train1 = 1;

    // no binning
    while(index1 < binStart){ // no binning
        kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
        index1 ++;
        train1 ++;
    }
    
    // bin range starting from train2
    x_end = train1 * x_factor;

    // log binning
    while(train1 <= num_event_temp){
        if(train1 < x_end){ 
            kernel2D_yLog(kernel2DLog, nums, kernel2D, num_event_temp, binStart, x_factor, index1, train1);
            train1 ++;
        }
        else{
            x_end *= x_factor;
            index1 ++;
        }
    }
        
    // print
    strcpy(prefix, "kernel2D");
    sprintf(output, "%s%s_%sLog.txt", folder, prefix, filename);
    file_out = fopen(output, "w");

    for(i = 0; i < num_bin; i ++){
        for(j = 0; j < num_bin; j ++){
            num = nums[i][j];
            if(num){
                avg = kernel2DLog[i][j] / num;
            }
            else avg = 0;

            fprintf(file_out, "%lf %lf %.10lf\n", mid_bin[i], mid_bin[j], avg);
        }
        fprintf(file_out, "\n");
    }
    fclose(file_out);

    free_matrix_double(kernel2D, 0, num_event_temp, 1, num_event_temp);
    free_vector_double(mid_bin, 0, num_event_temp);
    free_matrix_double(kernel2DLog, 0, num_bin, 0, num_bin);
    free_matrix_long(nums, 0, num_bin, 0, num_bin);
}
