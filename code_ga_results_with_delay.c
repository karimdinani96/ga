//GA Final - Delay Degradation 378K
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>

#define FILE_NM1         "bench/c432.bench"	//c432 c499 c880 c1355 c1908 c3540 c6288	c17
#define NO_OF_INPUT      36		   	//36   41   60   41    33    50    32		5
#define NO_OF_OUTPUT     07		   	//07   32   26   32    25    22    32		2
#define POPULATION_SIZE  100   //Create 100 arrays
#define Wd      0.0   //Weightage for Degradation
#define Wl      1.0   //Weightage for Leakage
#define stop 	200   //Stopping the GA


int time1=0;
struct node{
    int     color;
    int     start_time;
    int     finish_time;
    char    type[30];
    char    name[30];
    int     no_of_input;
    char    name_input[30][30];
    unsigned int     input_val[30];
    unsigned int     output_val_new;
    unsigned int     output_val_prev;
    struct  node *ptr[POPULATION_SIZE]; //front direction
    struct  node *prev[4];  //back direction
    int     no_of_output;
    int     order_evaluation;
    int     node_no;
    float   leak;
    float   degr;
    int     dummy;
};
struct node node_array[25000];

char    arr[1000000]; //To read benchmark ckt
struct  tms buf1, buf2; //Timing Calculation
float   nand2_leak[4],nand3_leak[8],nand4_leak[16],nor2_leak[4],nor3_leak[8],nor4_leak[16],inv_leak[2]; //Leakage Values of each gate
float   nand2_deg[4],nand3_deg[8],nand4_deg[16],nor2_deg[4],nor3_deg[8],nor4_deg[16],inv_deg[2]; //Deg values of each gate
int     count=0; //Used for making ckt
char    output_name[NO_OF_OUTPUT][30];//Used for making ckt
int     output_node[NO_OF_OUTPUT]; //Used for making ckt
int     out=0; //Used for making ckt
int     ar[25000]; //Used for making ckt -> topologicalsort
double   leakage_power = 0.0; //Leakage power calculation
double   degradation_value = 0.0; //Degradation calculation
double   max_deg, max_leak;
int stop_count = 0;

void initilize_node();
void leak_read();
void deg_read();
void make_graph(int start);
void topologicalsort_f();
void ga();
void dfs_visit_f(int i);
void make_interconnect(char str1[30],int c);

void evaluation_leak(char c[NO_OF_INPUT+1]);
void evaluate_output_leak(int len);
void calc_leakage_nand(int len,int no);
void calc_leakage_nor(int len,int no);
void calc_leakage_inv(int len,int no);

void evaluation_deg(char c[NO_OF_INPUT+1]);
void evaluate_output_deg(int len);
void calc_deg_nand(int len,int no);
void calc_deg_nor(int len,int no);
void calc_deg_inv(int len,int no);

void sort(char sort_array[POPULATION_SIZE][NO_OF_INPUT+1],double calculated_value[POPULATION_SIZE]);
void sort_opt(char sort_array[POPULATION_SIZE][NO_OF_INPUT+1],double calculated_value[POPULATION_SIZE],double deg_opt[POPULATION_SIZE],double leak_opt[POPULATION_SIZE]);
void fitness_evaluation(char opt_input[POPULATION_SIZE][NO_OF_INPUT+1],double leak_calc_copy[POPULATION_SIZE],double deg_calc_copy[POPULATION_SIZE],double opt_calc[POPULATION_SIZE]);
void next_gen(char opt_input[POPULATION_SIZE][NO_OF_INPUT+1],char opt_input_nextgen[POPULATION_SIZE][NO_OF_INPUT+1],double opt_calc[POPULATION_SIZE],double opt_calc_nextgen[POPULATION_SIZE]);


//--------------------------------------------------------------------------------
//Assign leakage
float nand2_leak[4]={178,833.9,746.8,1947};
float nand3_leak[8]={98.96,178.4,178.6,834.2,170.4,747.1,735.4,2921};
float nand4_leak[16]={69.87,99.29,99.35,178.6,99.82,178.9,179.2,834.5,94.20,170.7,171,747.5,170.3,735.7,729.4,3895};
float nor2_leak[4]={1667,974,794.5,130.9};
float nor3_leak[8]={2501,974.2,794.9,131.5,780.1,131.9,132.9,69.25};
float nor4_leak[16]={3334,974.5,795.2,131.6,780.4,132.8,133.2,69.65,980,974.5,795.1,131.8,780.4,132.3,133.2,69.71};
float inv_leak[2]={833.6,973.7};

//Assign Degradation
float nand2_deg[4]={0.1934,0.1002,0.1002,0.0093};
float nand3_deg[8]={0.2619,0.1761,0.1761,0.0927,0.1761,0.0927,0.0927,0.0125};
float nand4_deg[16]={0.3266,0.2478,0.2478,0.1692,0.2478,0.1692,0.1692,0.0906,0.2478,0.1692,0.1692,0.0906,0.1692,0.0906,0.0906,0.0158};
float nor2_deg[4]={0.2344,0.1172,0.0122,0.0043};
float nor3_deg[8]={0.3789,0.2526,0.1334,0.1262,0.0206,0.1156,0.0112,0.0046};
float nor4_deg[16]={0.4177,0.3131,0.2147,0.2089,0.1214,0.2004,0.1136,0.1082,0.3031,0.2089,0.1104,0.1044,0.0170,0.0956,0.0092,0.0038};
float inv_deg[2]={0.1333,0.0064};

//--------------------------------------------------------------------------------


//
int     br[25000];
char    brr[3000];
unsigned int transition=0;
char    best[30][NO_OF_INPUT+1];
//


int main()
{
    int     start = 0, i = 0;
    struct  node *ptr = NULL;
    FILE    *fp;
    times(&buf1); //get process time
    
    initilize_node();
    
    fp = fopen(FILE_NM1,"r");
    if(!feof(fp)){
        do{
            fscanf(fp,"%c",&arr[i]);
            i++;
            }while(!feof(fp));
        arr[i]='\0';
        fclose(fp);
    }
    else
        perror("FF");
    i=0;

    start=0;
    make_graph(start);
	
    ga();
	
    times(&buf2);
    printf("\nTime = %f\n", (float)((buf2.tms_utime - buf1.tms_utime) + (buf2.tms_stime - buf1.tms_stime)) / 1000000);
    return 1;
}

void initilize_node()
{
    int i;
    for(i=0;i<25000;i++){
        node_array[i].color=0;
        node_array[i].dummy=0;
        node_array[i].no_of_input=0;
        node_array[i].no_of_output=0;
        node_array[i].output_val_new=0;
    }
}


//Make the ckt----------------------------------------------------------------------------------------------------------------------------------------------------
void make_graph(int start)
{
    int     i,j,k,l,max_fan_in=0,max_fan_out=0;
    float   max_asap=0.0;
    char    str1[30],str2[30];
    struct  node    *ptr=NULL;

    for(i=start;arr[i]=='i';){
        while(arr[i]!='\n'){
            j=0;
            while(arr[i]!='(')
                str1[j++]=arr[i++];
            i++;
            str1[j]='\0';
            j=0;
            while(arr[i]!=')')
                str2[j++]=arr[i++];
            i++;
            str2[j]='\0';

            strcpy(node_array[count].type,str1);
            strcpy(node_array[count].name,str2);
            node_array[count].node_no=count;
            count++;
            while(arr[i]!='\n')i++;
        }
        i++;
    }

    while(arr[i++]=='\n');
    i--;

    while(arr[i]=='o'){
        j=0;
        while(arr[i++]!='(');
            while(arr[i]!=')')
                output_name[out][j++]=arr[i++];
        output_name[out][j]='\0';
        out++;
        while(arr[i]!='\n')
            i++;
        i++;
    }

    while(arr[i++]=='\n');
    i--;

    while(arr[i]!='\0'){
        while(1){
            j=0;
            while(arr[i]!='=')
                str1[j++]=arr[i++];
            i=i+2;
            j--;
            str1[j]='\0';
            j=0;
            while(arr[i]!='(')
                str2[j++]=arr[i++];
            i++;
            str2[j]='\0';
            strcpy(node_array[count].type,str2);
            strcpy(node_array[count].name,str1);
            node_array[count].node_no=count;

            node_array[count].no_of_input=1;
            node_array[count].output_val_new=0;
            j=k=0;
            while(arr[i]!=')'){
                if(arr[i]==','){
                    i=i+2;
                    str1[j]='\0';
                    j=0;
                    strcpy(node_array[count].name_input[k++],str1);
                    node_array[count].no_of_input++;
                }
                else
                    str1[j++]=arr[i++];
            }
            str1[j]='\0';
            strcpy(node_array[count].name_input[k],str1);
            count++;
            while(arr[i++]!='\n');
                break;
        }
    }
    for(i=0;i<count;i++){
        for(l=0;l<NO_OF_OUTPUT;l++)
            if(!(strcmp(output_name[l],node_array[i].name))){
                output_node[l]=i;
                break;
            }
    }

    for(i=NO_OF_INPUT;i<count;i++)
        for(j=0;j<node_array[i].no_of_input;j++)
            make_interconnect(node_array[i].name_input[j],i);
    
	topologicalsort_f();
}
	
void make_interconnect(char  str1[30],int    c)
{
    int i;
    for(i=0;i<count;i++)
        if(!(strcmp(str1,node_array[i].name))){
            node_array[i].ptr[node_array[i].no_of_output]=&node_array[c];
            node_array[c].prev[node_array[c].dummy++]=&node_array[i];
            node_array[i].no_of_output++;
            break;
    }
}

void topologicalsort_f()
{
    int i,j,temp1,temp2;
    for(i=0;i<count;i++){
        ar[i]=i;
        if(node_array[i].color==0)
            dfs_visit_f(i);
    }
    for(i=count-1;i>=0;i--){
        for(j=1;j<=i;j++){
            if(node_array[j-1].finish_time<node_array[j].finish_time){
                temp1=node_array[j-1].finish_time;
				temp2=ar[j-1];
                node_array[j-1].finish_time=node_array[j].finish_time;
                ar[j-1]=ar[j];

                node_array[j].finish_time=temp1;
                ar[j]=temp2;
            }
        }
    }
}

void dfs_visit_f(int i)
{
    int j;

    node_array[i].color=1;
    time1=time1+1;
    node_array[i].start_time=time1;
    for(j=0;j<node_array[i].no_of_output;j++)
        if(node_array[i].ptr[j]->color==0)
            dfs_visit_f(node_array[i].ptr[j]->node_no);
    time1=time1+1;
    node_array[i].finish_time=time1;
    node_array[i].color=2;
}


//Ckt is done----------------------------------------------------------------------------------------------------------------------------------------------------


//Genetic Algorithm------------------------------------------------------------------------------------------------------------------------------------------------
void ga()
{
    char  leak_input[POPULATION_SIZE][NO_OF_INPUT+1]; //Inputs for Leakage
    double  leak_calc[POPULATION_SIZE],leak_calc_copy[POPULATION_SIZE]; //Leakages value array
    double  zeroleak,oneleak; //All 0's and all 1's leakage
	
    char  deg_input[POPULATION_SIZE][NO_OF_INPUT+1]; //Inputs for Degradation
    double  deg_calc[POPULATION_SIZE],deg_calc_copy[POPULATION_SIZE]; //Degradation values array
    double  zerodeg,onedeg; //All 0's and all 1's Degradation

    char  opt_input[POPULATION_SIZE][NO_OF_INPUT+1], opt_input_nextgen[POPULATION_SIZE][NO_OF_INPUT+1]; //Optimized Inputs
    double  opt_calc[POPULATION_SIZE], opt_calc_nextgen[POPULATION_SIZE]; //Optimized values array
    double  zeroopt,oneopt; //Optimized values of all 0's and all 1's
//_copy are in unsorted format

    int   i,j,k;
    int   kk=0;
    float   random;

    //All 0's and all 1's leakage and degradation is calculated------------------
    //First Calculate all 0's and all 1's leakage
    srand((unsigned)time(NULL));
	//Assign one array to all 0's input and one array to all 1's input
	for(i=0;i<NO_OF_INPUT;i++){
	    leak_input[0][i]='0';
        leak_input[1][i]='1';
	}

    leak_input[0][i]='\0';
    leak_input[1][i]='\0';
	//Evaluate all 0's and all 1's leakage
    evaluation_leak(leak_input[0]);
    zeroleak = leakage_power; //All 0's leakage value
    leakage_power = 0.0;
    evaluation_leak(leak_input[1]);
    oneleak = leakage_power; //All 1's leakage value
	
	//First Calculate all 0's and all 1's Degradation
	//Assign one array to all 0's input and one array to all 1's input
	for(i=0;i<NO_OF_INPUT;i++){
	    deg_input[0][i]='0';
        deg_input[1][i]='1';
	}

    deg_input[0][i]='\0';
    deg_input[1][i]='\0';
    //Evaluate all 0's and all 1's degradation
    evaluation_deg(deg_input[0]);
    zerodeg = degradation_value; //All 0's degradation value
    degradation_value = 0.0;
    evaluation_deg(deg_input[1]);
    onedeg = degradation_value; //All 1's degradation value
    printf("All 1's degradation - %lf, All 0's Degradation - %lf\n",onedeg,zerodeg);
    printf("All 1's leakage - %lf, All 0's leakage - %lf\n",oneleak,zeroleak);
//----------------------------------------------------------------------------

	//Initial Population generation for leakage
    	for(i=2;i<POPULATION_SIZE;i++){
		for(j=0;j<NO_OF_INPUT;j++){
			random=drand48();
			if(random>0.5){
				leak_input[i][j]='1';
				deg_input[i][j]='1';
			}
			else{
				leak_input[i][j]='0';
				deg_input[i][j]='0';
			}
		}
		leak_input[i][j]='\0';
		for(k=0;k<i;k++)
			if(!(strcmp(leak_input[k],leak_input[i])))
				i-=1;
	}

	//Sort initial population acc to leakage
	printf("Leakage values\n");
	for(kk=0;kk<count;kk++)
        node_array[kk].leak=0.0;
	for(i=0;i<POPULATION_SIZE;i++){
        	leakage_power=0.0;
        	evaluation_leak(leak_input[i]);
        	leak_calc[i]=leakage_power;
        	leak_calc_copy[i]=leakage_power;
		printf("%s	%lf  \n",leak_input[i],leak_calc[i]);
	}
	for(int i=0;i<POPULATION_SIZE;i++)
		strcpy(opt_input[i],leak_input[i]);
	sort(leak_input,leak_calc);
	max_leak = leak_calc[POPULATION_SIZE-1];

	//Sort initial population acc to degradation
	printf("Degradation values\n");
	for(kk=0;kk<count;kk++)
        node_array[kk].degr=0.0;
	for(i=0;i<POPULATION_SIZE;i++){
        	degradation_value=0.0;
        	evaluation_deg(deg_input[i]);
        	deg_calc[i]=degradation_value;
        	deg_calc_copy[i]=degradation_value;
		printf("%s	%lf  \n",deg_input[i],deg_calc[i]);
    	}
	sort(deg_input,deg_calc);
	max_deg = deg_calc[POPULATION_SIZE-1];

	printf("Minimum leakage = %s  %lf\n",leak_input[0],leak_calc[0]);
	printf("Maximum leakage = %s  %lf\n",leak_input[POPULATION_SIZE-1],leak_calc[POPULATION_SIZE-1]);
	printf("Minimum degradation = %s	%lf  \n",deg_input[0],deg_calc[0]);
	printf("Maximum degradation = %s	%lf  \n",deg_input[POPULATION_SIZE-1],deg_calc[POPULATION_SIZE-1]);
	//All 0 and all 1 optimized values
    	zeroopt = (Wd*(zerodeg/max_deg)+Wl*(zeroleak/max_leak));
    	oneopt = (Wd*(onedeg/max_deg)+Wl*(oneleak/max_leak));
	printf("All 0 cooptimized = %lf  \n",zeroopt);
	printf("All 1 cooptimized = %lf  \n",oneopt);
	
	printf("Cooptimized values of Initial Population\n\n");
	//Calculate the cooptimized values and sort the arrays
	fitness_evaluation(opt_input,leak_calc_copy,deg_calc_copy,opt_calc);
	for(int i=0;i<POPULATION_SIZE;i++){
		printf("Inputs : %s	%lf  \n",opt_input[i],opt_calc[i]);
	}
	//Initial Population is Done
	
	//Next Generations are made and evaluated
	printf("GA Starts\n");
	for(i=0;i>=0;i++){
		next_gen(opt_input,opt_input_nextgen,opt_calc,opt_calc_nextgen);
		printf("Generation: %d  input: %s  %lf\n",i,opt_input_nextgen[0],opt_calc_nextgen[0]);
		if(stop_count == stop) break;
	}

}
//GA is over-------------------------------------------------------------------------------------------------------------------------------------------------------


//Leakage Evaluation of ckt-----------------------------------------------------------------------------------------------------------------------------------------------
void evaluation_leak(char c[NO_OF_INPUT+1])
{
    int i=0,j=0,k,l;

    for(j=0;j<NO_OF_INPUT;j++){
        if(c[j]=='1')
            node_array[j].input_val[0]=1;
        else if(c[j]=='0')
            node_array[j].input_val[0]=0;
    }

    evaluate_output_leak(1);
}

void evaluate_output_leak(int len)
{
    int i,j;
    unsigned int k;
    struct node *p;
    for(i=0;i<count;i++)
        node_array[i].order_evaluation=0;
    for(i=0;i<count;i++){
        p=&node_array[ar[i]];//Check mark

		if(!(strcmp(p->type,"nor")) || !(strcmp(p->type,"NOR"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=p->input_val[0];
			for(j=1;j<p->no_of_input;j++)
				p->output_val_new=p->output_val_new|p->input_val[j];
			p->output_val_new=~(p->output_val_new);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_leakage_nor(len,p->node_no);
		}
		else if(!(strcmp(p->type,"nand")) || !(strcmp(p->type,"NAND"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=p->input_val[0];
			for(j=1;j<p->no_of_input;j++)
				p->output_val_new=p->output_val_new&p->input_val[j];
			p->output_val_new=~(p->output_val_new);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_leakage_nand(len,p->node_no);
		}

		else if(!(strcmp(p->type,"not")) || !(strcmp(p->type,"NOT"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=~(p->input_val[0]);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_leakage_inv(len,p->node_no);
		}
		else{
			p->output_val_new=p->input_val[0];
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
		}
    }
}

//Calculate Each Gate leakage-------------------------------------------------------------------
void calc_leakage_nand(int len,int no)
{
    int i,j,k;
    unsigned int digit[32];

    for(j=len-1;j>=0;j--){
        digit[j]=0;
        for(i=0;i<node_array[no].no_of_input;i++){
            k=(node_array[no].input_val[i]>>j)&1;
            digit[j]+=((unsigned int)pow(2,i))*k;
        }
    }
    for(i=0;i<len;i++){
        if(node_array[no].no_of_input==2){
            leakage_power+=(double)nand2_leak[digit[i]];
            //printf("\n leakage nand2 %lf",nand2_leak[digit[i]]);
        }
        else if(node_array[no].no_of_input==3){
            leakage_power+=(double)nand3_leak[digit[i]];
            //printf("\n leakage nand3 %lf",nand3_leak[digit[i]]);
        }
        else{
            leakage_power+=(double)nand4_leak[digit[i]];
            //printf("\n leakage nand4 %lf",nand4_leak[digit[i]]);
        }
    }
}

void calc_leakage_nor(int len, int no)
{
    int i,j,k;
    unsigned int digit[32];

    for(j=len-1;j>=0;j--){
        digit[j]=0;
        for(i=0;i<node_array[no].no_of_input;i++){
            k=(node_array[no].input_val[i]>>j)&1;
            digit[j]+=((unsigned int)pow(2,i))*k;
        }
    }
    for(i=0;i<len;i++){
        if(node_array[no].no_of_input==2){
            leakage_power+=(double)nor2_leak[digit[i]];
            //printf("\n leakage nor2 %lf",nor2_leak[digit[i]]);
        }
        else if(node_array[no].no_of_input==3){
            leakage_power+=(double)nor3_leak[digit[i]];
            //printf("\n leakage nor3 %lf",nor3_leak[digit[i]]);
        }
        else{
            leakage_power+=(double)nor4_leak[digit[i]];
            //printf("\n leakage nor4 %lf",nor4_leak[digit[i]]);
        }
    }
}

void calc_leakage_inv(int len, int no)
{
    int i,j,k;
    j=~(node_array[no].output_val_new);
    for(i=len-1;i>=0;i--){
        k=(j>>i)&1;
        leakage_power+=(double)inv_leak[k];
        //printf("\n leakage inv %lf",inv_leak[k]);
    }
}

//Each Gate leakage Calculation is done-----------------------------------------------------------

//Leakage Evaluation of ckt is done--------------------------------------------------------------------------------------------------------------------------



//Degradation Evaluation of ckt-----------------------------------------------------------------------------------------------------------------------------------------------
void evaluation_deg(char c[NO_OF_INPUT+1])
{
    int i=0,j=0,k,l;
	
    for(j=0;j<NO_OF_INPUT;j++){
        if(c[j]=='1')
            node_array[j].input_val[0]=1;
        else if(c[j]=='0')
            node_array[j].input_val[0]=0;
    }

    evaluate_output_deg(1);
}

void evaluate_output_deg(int len)
{
    int i,j;
    unsigned int k;
    struct node *p;
    for(i=0;i<count;i++)
        node_array[i].order_evaluation=0;
    for(i=0;i<count;i++){
        p=&node_array[ar[i]];//Check Mark

		if(!(strcmp(p->type,"nor")) || !(strcmp(p->type,"NOR"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=p->input_val[0];
			for(j=1;j<p->no_of_input;j++)
				p->output_val_new=p->output_val_new|p->input_val[j];
			p->output_val_new=~(p->output_val_new);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_deg_nor(len,p->node_no);
		}
		else if(!(strcmp(p->type,"nand")) || !(strcmp(p->type,"NAND"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=p->input_val[0];
			for(j=1;j<p->no_of_input;j++)
				p->output_val_new=p->output_val_new&p->input_val[j];
			p->output_val_new=~(p->output_val_new);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_deg_nand(len,p->node_no);
		}

		else if(!(strcmp(p->type,"not")) || !(strcmp(p->type,"NOT"))){
			p->output_val_prev=p->output_val_new;
			p->output_val_new=~(p->input_val[0]);
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
			calc_deg_inv(len,p->node_no);
		}
		else{
			p->output_val_new=p->input_val[0];
			for(j=0;j<p->no_of_output;j++)
				p->ptr[j]->input_val[p->ptr[j]->order_evaluation++]=p->output_val_new;
		}
    }
}

//Calculate Each Gate Degradation-------------------------------------------------------------------
void calc_deg_nand(int len,int no)
{
    int i,j,k;
    unsigned int digit[32];

    for(j=len-1;j>=0;j--){
        digit[j]=0;
        for(i=0;i<node_array[no].no_of_input;i++){
            k=(node_array[no].input_val[i]>>j)&1;
            digit[j]+=((unsigned int)pow(2,i))*k;
        }
    }
    for(i=0;i<len;i++){
        if(node_array[no].no_of_input==2){
            degradation_value+=(double)nand2_deg[digit[i]];
            //printf("\n Degradation nand2 %lf",nand2_deg[digit[i]]);
        }
        else if(node_array[no].no_of_input==3){
            degradation_value+=(double)nand3_deg[digit[i]];
            //printf("\n Degradation nand3 %lf",nand3_deg[digit[i]]);
        }
        else{
            degradation_value+=(double)nand4_deg[digit[i]];
            //printf("\n Degradation nand4 %lf",nand4_deg[digit[i]]);
        }
    }
}

void calc_deg_nor(int len, int no)
{
    int i,j,k;
    unsigned int digit[32];

    for(j=len-1;j>=0;j--){
        digit[j]=0;
        for(i=0;i<node_array[no].no_of_input;i++){
            k=(node_array[no].input_val[i]>>j)&1;
            digit[j]+=((unsigned int)pow(2,i))*k;
        }
    }
    for(i=0;i<len;i++){
        if(node_array[no].no_of_input==2){
            degradation_value+=(double)nor2_deg[digit[i]];
            //printf("\n Degradation nor2 %lf",nor2_deg[digit[i]]);
        }
        else if(node_array[no].no_of_input==3){
            degradation_value+=(double)nor3_deg[digit[i]];
            //printf("\n Degradation nor3 %lf",nor3_deg[digit[i]]);
        }
        else{
            degradation_value+=(double)nor4_deg[digit[i]];
            //printf("\n Degradation nor4 %lf",nor4_deg[digit[i]]);
        }
    }
}

void calc_deg_inv(int len, int no)
{
    int i,j,k;
    j=~(node_array[no].output_val_new);
    for(i=len-1;i>=0;i--){
        k=(j>>i)&1;
        degradation_value+=(double)inv_deg[k];
        //printf("\n Degradation inv %lf",inv_deg[k]);
    }
}

//Each Gate degradation Calculation is done-----------------------------------------------------------

//Degradation Evaluation of ckt is done--------------------------------------------------------------------------------------------------------------------------

//Sort the leakage and degradation values
void sort(char sort_array[POPULATION_SIZE][NO_OF_INPUT+1], double calculated_value[POPULATION_SIZE])
{
    int i,j;
    double temp;
    char c[NO_OF_INPUT+1];

    c[0]='\0';
    for(i=0;i<POPULATION_SIZE-1;i++){
        for(j=0;j<POPULATION_SIZE-1-i;j++){
            if(calculated_value[j] >= calculated_value[j+1]){
                temp = calculated_value[j];
                strcpy(c,sort_array[j]);

                calculated_value[j] = calculated_value[j+1];
                strcpy(sort_array[j],sort_array[j+1]);

                calculated_value[j+1] = temp;
                strcpy(sort_array[j+1],c);
            }
	    else
		    continue;
        }
    }
    //for(i=0;i<POPULATION_SIZE;i++){
	    //printf("Input = %s	%lf\n",sort_array[i],calculated_value[i]);
    //}
}

//Sort the co-optimized values with leakage and degradation values
void sort_opt(char sort_array[POPULATION_SIZE][NO_OF_INPUT+1],double calculated_value[POPULATION_SIZE],double deg_opt[POPULATION_SIZE],double leak_opt[POPULATION_SIZE])
{
    int i,j;
    double temp_opt_val;
    double temp_deg_val;
    double temp_leak_val;
    char temp_opt[NO_OF_INPUT+1];

    temp_opt[0]='\0';
    for(i=0;i<POPULATION_SIZE-1;i++){
        for(j=0;j<POPULATION_SIZE-1-i;j++){
            if(calculated_value[j] >= calculated_value[j+1]){
                //Sort the array acc to opt_values
		temp_opt_val = calculated_value[j];
                strcpy(temp_opt,sort_array[j]);

                calculated_value[j] = calculated_value[j+1];
                strcpy(sort_array[j],sort_array[j+1]);

                calculated_value[j+1] = temp_opt_val;
                strcpy(sort_array[j+1],temp_opt);
		
		//Copy the leakage and deg of corresponding opt_values
		temp_deg_val = deg_opt[j];
		temp_leak_val = leak_opt[j];
		deg_opt[j] = deg_opt[j+1];
		leak_opt[j] = leak_opt[j+1];
		deg_opt[j+1] = temp_deg_val;
		leak_opt[j+1] = temp_leak_val;
            }
	    else
		    continue;
        }
    }
    //for(i=0;i<POPULATION_SIZE;i++){
	    //printf("Input = %s	%lf\n",sort_array[i],calculated_value[i]);
    //}
    printf("Leakage: %lf  Degradation: %lf\n",leak_opt[0],deg_opt[0]);
}

//Calculate the cooptimized values and sort the array
void fitness_evaluation(char opt_input[POPULATION_SIZE][NO_OF_INPUT+1],double leak_calc_copy[POPULATION_SIZE],double deg_calc_copy[POPULATION_SIZE],double opt_calc[POPULATION_SIZE]){
	for(int i=0;i<POPULATION_SIZE;i++){
		opt_calc[i] = (Wl*(leak_calc_copy[i]/max_leak)+Wd*(deg_calc_copy[i]/max_deg));
		//printf("%s  %lf\n",opt_input[i],opt_calc[i]);
	}
	//sort(opt_input,opt_calc);
	sort_opt(opt_input,opt_calc,deg_calc_copy,leak_calc_copy);
}


//Next Generation
void next_gen(char opt_input[POPULATION_SIZE][NO_OF_INPUT+1],char opt_input_nextgen[POPULATION_SIZE][NO_OF_INPUT+1],double opt_calc[POPULATION_SIZE],double opt_calc_nextgen[POPULATION_SIZE]){
	
	double leakage[POPULATION_SIZE],degradation[POPULATION_SIZE];
	//We'll do all operation on this then we will copy
	char temp_cross_in1[NO_OF_INPUT+1];
	char temp_cross_in2[NO_OF_INPUT+1];
	char temp_char[NO_OF_INPUT+1];

	//Direct Copy best 10%
	for(int i=0;i<(POPULATION_SIZE* 0.1);i++){
		strcpy(opt_input_nextgen[i],opt_input[i]);
		//printf("%s\n",opt_input_nextgen[i]);
	}
	
	//Multi piont Crossover - 2 point - 80%
	for(int i=(POPULATION_SIZE* 0.1);i<(POPULATION_SIZE* 0.9);i=i+2){
		//Now select the 2 inputs on which we want to do crossover
		//Generate a random number. If no>0.5 then take 2 inputs from best 10 else take from rest. Crosscheck for repeatation of input combinations
		float random = drand48();
		if(random>0.5){
			int random_int = drand48()*100;
		        random_int = random_int %10;
			strcpy(temp_cross_in1,opt_input[random_int]);
			random_int = drand48()*100;
		        random_int = random_int %10;
			strcpy(temp_cross_in2,opt_input[random_int]);
		}
		else{
			int random_int = drand48()*100;
			strcpy(temp_cross_in1,opt_input[random_int]);
			random_int = drand48()*100;
			strcpy(temp_cross_in2,opt_input[random_int]);
		}
		//Generate 3 crossover points in between min count and max count and swap 2nd and 4th part
		int MAX_COUNT = NO_OF_INPUT *(0.8);
		int MIN_COUNT = NO_OF_INPUT *(0.2);
		int point1, point2, point3;
		point1 = rand() %MAX_COUNT + MIN_COUNT;
		point2 = rand() %MAX_COUNT + MIN_COUNT;
		//Sort these 3 points in ascending order - point1 < point2 < point3
		if(point1 > point2){
			int swap;
			swap = point2;
			point2 = point1;
			point1 = swap;
		}
		//printf("%d  %d  %d\n",point1,point2,point3);
		//Copy 1st part, swap 2nd, copy 3rd
		for(int j=point1;j<point2;j++){
			temp_char[j] = temp_cross_in1[j];
			temp_cross_in1[j] = temp_cross_in2[j];
			temp_cross_in2[j] = temp_char[j];
		}
		//Crossover is done. Now we will check whether this input already exists in the existing array or not. If exists, then regenerate.
		int flag = 0;
		for(int j=0;j<POPULATION_SIZE;j++){
			if((strcmp(opt_input[j],temp_cross_in1) == 0) || (strcmp(opt_input[j],temp_cross_in2) == 0) || (strcmp(temp_cross_in1,temp_cross_in2) == 0)){
				i-=2;
				flag = 1;
				break;
			}
			else if((j == POPULATION_SIZE-1) && (flag == 0)){
				strcpy(opt_input_nextgen[i],temp_cross_in1);
				strcpy(opt_input_nextgen[i+1],temp_cross_in2);
			}
		}
	}

	//Multi point Mutation - 10%
	for(int i=(POPULATION_SIZE* 0.9);i<(POPULATION_SIZE);i=i+1){
		//Now select the 1 input on which we want to do mutation
		//Generate a random number. If no>0.5 then take the input from best 10 else take from rest. Crosscheck for repeatation of input combinations
		float random = drand48();
		if(random>0.5){
			int random_int = drand48()*100;
		        random_int = random_int %10;
			strcpy(temp_char,opt_input[random_int]);
		}
		else{
			int random_int = drand48()*100;
			strcpy(temp_char,opt_input[random_int]);
		}
		//Now do the mutation at 2 to 20% of the NO_OF_INPUTS
		//Generate an random no that tells at howmany positions we want to do the mutation then generate those positions
		int MAX_COUNT = NO_OF_INPUT *(0.2);
		//printf("MAX_COUNT = %d\n",MAX_COUNT);
		int random_no = rand() %MAX_COUNT + 2;//Atleast 2 positions for mutation
		//printf("random_no = %d\n",random_no);
		int position[random_no];
		for(int i=0;i<random_no;i++){
			position[i] = rand() %NO_OF_INPUT;
			//printf("position = %d\n",position[i]);
			if(temp_char[position[i]] == '1')
				temp_char[position[i]] = '0';
			else
				temp_char[position[i]] = '1';
		}
		//Mutation is done. Now we will check whether this input already exists in the existing array or not. If exists, then regenerate.
		int flag = 0;
		for(int j=0;j<POPULATION_SIZE;j++){
			if(strcmp(opt_input[j],temp_char) == 0){
				i-=1;
				flag = 1;
				break;
			}
			else if((j == POPULATION_SIZE-1) && (flag == 0)){
				strcpy(opt_input_nextgen[i],temp_char);
			}
		}
	}
	
	//The next gen is raedy. Now Calculate its leakage, degradation and cooptimized value and sort it
	for(int kk=0;kk<count;kk++)
        node_array[kk].leak=0.0;
	for(int i=0;i<POPULATION_SIZE;i++){
        	leakage_power=0.0;
        	evaluation_leak(opt_input_nextgen[i]);
        	leakage[i]=leakage_power;
		//printf("%s	%lf  \n",leak_input[i],leak_calc[i]);
	}
	
	for(int kk=0;kk<count;kk++)
        node_array[kk].degr=0.0;
	for(int i=0;i<POPULATION_SIZE;i++){
        	degradation_value=0.0;
        	evaluation_deg(opt_input_nextgen[i]);
        	degradation[i]=degradation_value;
		//printf("%s	%lf  \n",deg_input[i],deg_calc[i]);
    	}
	fitness_evaluation(opt_input_nextgen,leakage,degradation,opt_calc_nextgen);
	//sort(opt_input_nextgen,opt_calc_nextgen);
	
	//Stopping Criteria
	if(opt_calc_nextgen[0] == opt_calc[0]){
		stop_count++;
	}
	if(opt_calc_nextgen[0] != opt_calc[0]){
		stop_count=0;
	}
	
	for(int i=0;i<POPULATION_SIZE;i++){
		strcpy(opt_input[i],opt_input_nextgen[i]);
		opt_calc[i] = opt_calc_nextgen[i];
	}
}
