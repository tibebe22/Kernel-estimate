// quick sort

//////////////////////////////////////// Working functions!

// reverse the list

void Reverse_list_double(double **a, long a_column, long low, long high){
	long i, left=low, right=high;
    double tmp; 

	while(left<right){
		for(i=0; i<=a_column; i++){
			tmp=a[left][i]; a[left][i]=a[right][i]; a[right][i]=tmp; 
        }
		left++, right--; 
	}
}

void QuickSort3_ascend_double(double **a, long a_column, long sort_col, long first, long last){
	long key, i, j, k, index, pivot; 
    double temp;

	if(first<last){
        pivot = first;
        i = first;
        j = last;

        while(i < j){
            while(a[i][sort_col] <= a[pivot][sort_col] && i < last) i ++;
            while(a[j][sort_col] > a[pivot][sort_col]) j --;
            if(i < j){
                for(index = 0; index <= a_column; index ++){
                    temp = a[i][index]; a[i][index] = a[j][index]; a[j][index] = temp;
                }
            }
        } 
        
        for(index = 0; index <= a_column; index ++){
            temp = a[pivot][index]; a[pivot][index] = a[j][index]; a[j][index] = temp;
        }
        QuickSort3_ascend_double(a, a_column, sort_col, first, j - 1);
        QuickSort3_ascend_double(a, a_column, sort_col, j + 1, last);
    }
}

void QuickSort3_descend_double(double **a, long a_column, long sort_col, long first, long last){
	QuickSort3_ascend_double(a, a_column, sort_col, first, last); 
	Reverse_list_double(a, a_column, first, last); 
}

