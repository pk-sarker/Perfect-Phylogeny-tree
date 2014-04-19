/* **********************************************************
 * ** Perfect Phylogeny Tree from Character State Matrix   **
 * **********************************************************
 * 
 * File:   perfectphy.c
 * Author: Pijus Kumar Sarker
 * 
 * 
 * - Algorithms in Bioinformatics
 * 
 *
 * 
 * 
 * +------------------------------------------------------------------------------------------------------+
 * | PROBLEM 1                                                                                            |
 * | ---------                                                                                            |
 * | Implement the algorithm for testing if a character-state matrix, with 2 states for each character,   |
 * | admits a prefect phylogeny. When it does, display the phylogeny tree. Even when it doesn't, you can  |
 * | display the phylogeny tree for the first i<=i<=m compatible characters.                              |
 * +------------------------------------------------------------------------------------------------------+
 * 
 * *****************
 * *  Test Input   *
 * *****************
 * Number of Taxon(n) : 5
 * Number of Characters in each Taxon (m): 5
 * 
 * Taxon 1 = 1 1 0 0 0
 * Taxon 2 = 0 0 1 0 0
 * Taxon 3 = 1 1 0 0 1
 * Taxon 4 = 0 0 1 1 0
 * Taxon 5 = 0 1 0 0 0
 * 
 * ***************
 * *   Output    *
 * ***************
 * 
 * 
               {1,2,3,4,5}
     -------------00100--------------
    1|                              |0
  01000{1,3,5}                    00100{2,4}


                  {1,3,5}
     -------------01000--------------
    1|                              |0
  11000{1,3}                       01000{5}


                  {2,4}
     -------------00100--------------
    1|                              |0
  00110{4}                        00100{2}


                  {1,3}
     -------------11000--------------
    1|                              |0
  11001{3}                       11000{1}
 * 
 * 
 * Above O/P is the representation of the tree bellow 
 *                 {1,2,3,4,5}
 *               1/           \0
 *           {1,3,5}          {2,4}
 *         1/      \0        1/   \0
 *       {1,3}     {5}      {4}   {2}
 *      1/    \0
 *     {3}     {1}
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


int N,M,isPhylo,vi=0;
int MAT[1000][1000],MAT1[1000][1000], O[1000][1000],v[1000][1000], TR[100][3][100];
int vs[1000][1000],vc=0,MO[1000][1000],s1[1000],s0[1000],p=0;

/*
 * **********************
 * ** Global Variables **
 * **********************
 * N : Number of Taxon
 * M : Number of Characters in each Taxon
 * isPhylo : Flag which indicates whether it is a perfect phylogeny or not
 * MAT[][] : Character-State Matrix
 * MAT1[][] : Character-State Matrix where the columns are ordered in descending base on number of 1's per column 
 * O[][] : O[i] set
 * vs[][] : Vertex set those are already divided by characters or O[i].
 * vc : Number of vertexes in vs[][]
 * MO[][] : O[i] set created for newly formed character-state matrix M'
 * s1[] : Set of taxon that has 1's in a particular column
 * s0[] : Set of taxon that has 0's in a particular column
 *  
 */
int main(int argc, char** argv) {

    int n,m, i,j,temp;
    int tm[3][100];
        
    printf("Number of Taxon (N): ");
    scanf("%d",&n);
    printf("\nNumber of Characters in each Taxon (M): ");
    scanf("%d",&m);
    
    N=n;
    M=m;
    char taxon[n][m];
    
    for(i=1; i<=n; i++){
        temp = i-1;
        printf("\nTaxon %d : ",i);
        scanf("%s",taxon[temp]);
    }
    
    printf("\n\n");
    
    for(i=0; i<N; i++){
        for(j=0; j<M; j++){
            MAT[i][j] = taxon[i][j]-'0';
        }
    }
    
    cal_o(); 
    /*** Display Character-State Matrix ***/
    displayMatrix();
    
    /*** Transformed Character-State Matrix based on number of 1's in columns ***/
    newmat();
    
    /*** Reconstruct the Phylogeny tree ***/
    conustTree();
    return (0);
    
}
/**
 * Function : cal_o()
 * Description: Creates the O[i] set 
 *              where 0<i<M, M = Number of characters in each taxon
 *              And check each pairwise compatibility of all O[i] 
 */
void cal_o(){
    int i,j;
    isPhylo=1;
    for(i=0; i<M; i++){
        int t=1;
        for(j=0; j<N; j++){
            if(MAT[j][i]==1){
                O[i][t] = j+1;
                t++;
            }
        }
        O[i][0] = --t;
    }
    
    for(i=1; i<=M; i++){
        for(j=i+1; j<=M; j++){
            if(disjoint(i,j)==1 || iscontained(i,j)==1){  // check whether O[i] and O[j] is compatible or not
            
            }else{
                 isPhylo=0;
            }
        }
    }
    
    // Display the character-state matrix represents a perfect phylogeny or not
    if(isPhylo==1){
        printf("\n\n ****** Its Perfect Phylogeny ****** \n");
    }else{
        printf("\n\n ****** Not Perfect Phylogeny ****** \n");
    }
}

/**
 * Function : displayMatrix()
 * Description: Display character-state matrix, MAT[i][j]
 *              where 0<i<N and 0<j<M, and 
 *              N = Number of Taxon, M = Number of characters in each taxon
 */
void displayMatrix(){
    int i, j;
    printf("\n\n***** MATRIX (M) *****\n");
    for(i=0; i<N; i++){
        int t=1;
        for(j=0; j<M; j++){
            printf("%d  ",MAT[i][j]);
        }
        printf("\n");
    }
    
}

/**
 * Function : newmat()
 * Description: This function creates a new Character-State Matrix 
 *              where the columns are ordered in descending base on 
 *              number of 1's per column.
 *              
 *              N = Number of Taxon, M = Number of characters in each taxon
 * Variables :
 * MAT1[][] : Transformed matrix
 * cnt[i][] : Number of 1's in i column
 * soor[i][] : Stores the index of column i
 *  
 */
void newmat(){
    int i,j,k,l,n,temp,cnt[1000][100] ,soor[1000][100];
    
    // Initialize MAT1[][] with 0
    for(i=0;i<M; i++){
        MAT1[i][0] = 0;
    }
    
    // Count number of 1's in each column
    for(i=0;i<M; i++){
        k=0;
        for(j=0;j<N;j++){
            if(MAT[j][i]==1){
                k++;
            }
        }
        soor[i][0] = i;
        cnt[i][0] = k;
    }
    
    // Sort the columns in descending order based on maximum number of 1's 
    for(i=0 ; i<M ; i++)
    {
        for(j=0 ; j<M-i-1 ; j++)
        {
            if(cnt[j][0]<cnt[j+1][0])
            {
                temp=cnt[j][0];
                cnt[j][0]=cnt[j+1][0];
                cnt[j+1][0]=temp;
                
                temp=soor[j][0];
                soor[j][0]=soor[j+1][0];
                soor[j+1][0]=temp;
            }
        }
    }
    
    // Create the new character-state matrix, MAT1
    for(i=0; i<M; i++){
        k=soor[i][0];
        for(j=0;j<N;j++){
            MAT1[j][i] = MAT[j][k];
        }
    }
    
    // Display New matrix, M'
    printf("\n\n***** MATRIX (M') *****\n");
    for(i=0; i<N; i++){
        int t=1;
        for(j=0; j<M; j++){
            printf("%d  ",MAT1[i][j]);
        }
        printf("\n");
    }
    
}

/**
 * Function : conustTree()
 * Description: This function creates a new Character-State Matrix 
 *              where the columns are ordered in descending base on 
 *              number of 1's per column.
 *              
 *              N = Number of Taxon, M = Number of characters in each taxon
 * Variables :
 * MAT1[][] : Transformed matrix
 * cnt[i][] : Number of 1's in i column
 * soor[i][] : Stores the index of column i
 *  
 */

void conustTree(){
    int i,j,ii,jj,k=0,f=0,c;
    printf("\n\n ******* CONSTRUCT TREE *******\n");
        s1[0]=0;
        s0[0]=0;
    for(i=1;i<=1000;i++){
        s1[i]=0;
        s0[i]=0;
    }
    for(i=0;i<M; i++){
        k=1;
        for(j=1;j<=N;j++){
            if(MAT1[j-1][i]==1){
                MO[i][k]=j;
                k++;
            }
            
        }
        MO[i][0] = k;
    }
    
    for(i=0; i<M; i++){
        if(vc==0){
            // include all to leaf to vs
            c=0;
            for(j=1; j<=N;j++){
                vs[vc][j]=j;
                c++;
            }
            vs[vc][0]=c;
            vc++;
        }
    
        f = 0;
        int fl=0,vcc=vc;
        for(j=0; j<vc; j++){
            f=isdividable(j,i);
            if(f==1){
                fl++;
                for(k=0; k<=vs[j][0];k++){
                    TR[p][0][k] = vs[j][k];
                }        
                
                for(k=1; k<=s1[0]; k++){
                    vs[j][k] = s1[k];
                    TR[p][1][k] = s1[k];
                }
                TR[p][1][0] = s1[0];
                vs[j][0] = s1[0];
                vs[j][k] = 1;
                
                for(k=1; k<=s0[0]; k++){
                    vs[vc][k] = s0[k];
                    TR[p][2][k] = s0[k];
                }
                TR[p][2][0] = s0[0];
                vs[vc][0] = s0[0];
                vs[vc][k] = 0;
                vcc++;
                p++;
            }    
        }
        vc = vcc;
    }
    
    for(ii=0;ii<p;ii++){
        if(vs[ii][0]>1 && isPhylo==1){
            for(jj=0;jj<=vs[ii][0];jj++){
                TR[p][0][jj] = vs[ii][jj]; 
                
                TR[p][1][0] = 1;
                TR[p][1][1] = vs[ii][1];
                TR[p][1][2] = vs[ii][vs[ii][0]+1];
                TR[p][2][0] = 1;
                TR[p][2][1] = vs[ii][2];
                TR[p][2][2] = vs[ii][vs[ii][0]+1];
                
            }
            p++;
        }
    }
    
    dispTree();
    
}


/************************************************
 * Function: dispTree
 * Description: Print phylogeny tree
 * Input Parameter: NULL
 * 
 * Variables: 
 * TR[][][] : stores the tree
 * For example:
 *          {1,2,3,4,5}
 *           1/     \0
 *       {1,3,5}   {2,4}
 * Then TR stores this like below,
 * TR[i][0][] = {5,1,2,3,4,5,x}  // first value is the length of the set
 * TR[i][1][] = {3,1,3,5,x}
 * TR[i][2][] = {2,2,4,x}
 * 
 * Here the last value(x) indicates the index of TR, where this vertex is divided in two sets
 * if it is leaf vertex then x will be 0, that means it has not more left or right child.
 * 
 * vtx[][][] : represnts the vertex codes,
 * 
 * For Example: if TR[i][1][] = {1,4,x} then vtx[i][1][] = [0 0 1 1 0] (forth taxon code)
 ***********************************************/

void dispTree(){
    int i,j,k,mi,fl,vtx[1000][3][500];
    for(i=0; i<p;i++){
        
        if(i==(p-1)){
            if(TR[i][1][0]==1){
                TR[i][1][2] =0;
            }
            if(TR[i][2][0]==1){
                TR[i][2][2]=0;
            }
        }else{
            for(j=i+1;j<p;j++){
                if(TR[i][1][0]==1){
                TR[i][1][2] =0;
                }else{
                    fl=1;
                    for(k=0;k<TR[i][1][0];k++){
                        if(TR[j][0][k]!=TR[i][1][k]){
                            fl=0;
                            break;
                        }
                    }
                    if(fl==1){
                        TR[i][1][k+1] = j;
                    }
                }

                if(TR[i][2][0]==1){
                    TR[i][2][2]=0;
                }else{
                    fl=1;
                    for(k=0;k<TR[i][2][0];k++){
                        if(TR[j][0][k]!=TR[i][2][k]){
                            fl=0;
                            break;
                        }
                    }
                    if(fl==1){
                        TR[i][2][k+1] = j;
                    }
                }
            }
        }
        
    }
    printf("\n\n");
    
    for(i=p-1; i>-1;--i){
        //Left child - 1
        if(TR[i][1][0]==1){
            k=TR[i][1][1]-1;
            for(j=0; j<M; j++){
                vtx[i][1][j] = MAT[k][j];
            }
        }else{
            k = TR[i][1][0]+1;
            k = TR[i][1][k];
            for(j=0; j<M; j++){
                vtx[i][1][j] = vtx[k][0][j];
            }
        }
        //Right child - 0
        if(TR[i][2][0]==1){
            k=TR[i][2][1]-1;
            for(j=0; j<M; j++){
                vtx[i][2][j] = MAT[k][j];
                vtx[i][0][j] = MAT[k][j]; // Parent
            }
            
        }else{
            k = TR[i][2][0]+1;
            k = TR[i][2][k];
            for(j=0; j<M; j++){
                vtx[i][2][j] = vtx[k][0][j];
                vtx[i][0][j] = vtx[k][0][j];
            }
        }
    }   
    
    for(i=0; i<p;i++){
        printSp(18);
        printf("{");
        for(j=1; j<=TR[i][0][0];j++){
            if(j>1){
                printf(",");
            }
            printf("%d",TR[i][0][j]);
        }
        printf("}\n");
        printSp(5);
        printDash(13);
        for(j=0; j<M; j++){
            printf("%d",vtx[i][0][j]);
        }
        printDash(14);
        printf("\n    1|                              |0\n  ");
        
        for(j=0; j<M; j++){
            printf("%d",vtx[i][1][j]);
        }
        printf("{");
        for(j=1; j<=TR[i][1][0];j++){
            if(j>1){
                printf(",");
            }
            printf("%d",TR[i][1][j]);
        }
        printf("}");
        printSp(23);
        
        //printf("");
        for(j=0; j<M; j++){
            printf("%d",vtx[i][2][j]);
        }
        printf("{");
        for(j=1; j<=TR[i][2][0];j++){
            if(j>1){
                printf(",");
            }
            printf("%d",TR[i][2][j]);
        }
        printf("}\n");
        printf("\n\n");
    }
            
}

/************************************************
 * Function: printSp
 * Description: Print spaces
 * Input Parameter: Number of spaces to print 
 ***********************************************/
void printSp(int j){
    int i;
    for(i=0; i<j; i++){
        printf(" ");
    }
}

/************************************************
 * Function: printDash
 * Description: Print dash symbol
 * Input Parameter: Number of dash to print 
 ***********************************************/
void printDash(int j){
    int i;
    for(i=0; i<j; i++){
        printf("-");
    }
}

/************************************************
 * Function: isdividable
 * Description: It checkes if any taxon set can dividable by any O[i]
 *              If true then it assigns the taxons with 1's in the colum to s1[] 
 *                      and the taxons with 0's in the colum to s0[]
 *               
 * Input Parameter: taxon set index(previously separated vertices sets or childs ),vi 
 *                  and O[i] set
 * Variables:
 * s0c : index of s0[] set
 * s1c : index of s1[] set
 * vs[][] : Vertex set, those are already separated/divided by characters(O[i])
 * MO[][] : New O[i] set after transforming M character-state matrix to M'
 * oi : index of newly formed O[i] set from M'
 * 
 * flag : Represents whether any vertex set vs[][] can be separated in to two groups(0 and 1) or not 
 * 
 ***********************************************/

int isdividable(int vi,int oi){
    int i,j,k=0, flag=0,s0c=1,s1c=1;
    
    if(MO[oi][0]<=vs[vi][0]){  // check the length of vertex set and O[i] set because vertex set like {2,4} can't be separated by O[i] = {2,4,5}
         k=0;
        for(i=1; i<MO[oi][0];i++){
           
            for(j=1;j<=vs[vi][0];j++){
                if(MO[oi][i]==vs[vi][j]){
                    s1[s1c++] = vs[vi][j];
                    k++;
                }
            }
        }
        s1[0] = --s1c;
        
        for(j=1;j<=vs[vi][0];j++){
            
            int ch=1;
            for(i=1;i<MO[oi][0];i++){
                if(vs[vi][j]==MO[oi][i]){
                    ch = 0;
                }
            }
            if(ch==1){
                s0[s0c++] = vs[vi][j];
            }
        }
        s0[0] = --s0c;
        if(k==MO[oi][0]-1){
            flag = 1;
        }else{
            flag = 0;
        }
    }else{
        flag=0;
    }
   
    return flag;
}

/************************************************
 * Function: disjoint
 * Description: Checks whether O[i] and O[j] is disjoint or not
 * Input Parameter: indexes of two O set.
 ***********************************************/

int disjoint(int oi, int oj){
    
    int i,j,disj = 1;
    for(i=1; i<=O[oi-1][0]; i++){
        for(j=1; j<=O[oj-1][0]; j++){
            if(O[oi-1][i]==O[oj-1][j]){
                disj=0;
            }
        }
    }
    return disj;
}

/************************************************
 * Function: iscontained
 * Description: Checks whether O[i] contained in O[j] or O[j] contained in O[i]
 * Input Parameter: indexes of two O set.
 ***********************************************/

int iscontained(int oi, int oj){
    int i,j, mfg=1;
    if(O[oi-1][0]==O[oj-1][0]){  // if both set are equal length, like O[1]={1,3} and O[2]={2,4}
        int flg =0;
        for(i=1;i<=O[oi-1][0]; i++){
            
            for(j=1; j<=O[oj-1][0]; j++){
                if(O[oi-1][i]==O[oj-1][j]){
                    flg++;
                }
            }
            
            
        }
        if(flg!=O[oi-1][0]){
            mfg = 0;
        }
    }else if(O[oi-1][0]>O[oj-1][0]){  // if O[oi] set is bigger that O[oj] set.
        
        int flg =0;
        for(i=1;i<=O[oj-1][0]; i++){
            for(j=1; j<=O[oi-1][0]; j++){
                if(O[oj-1][i]==O[oi-1][j]){
                    flg++;
                }
            }
        }
        
        if(flg!=O[oj-1][0]){
            mfg = 0;
        }        
    }else{     // if O[oj] set is bigger that O[oi] set.
        int flg =0;
        for(i=1;i<=O[oi-1][0]; i++){
            for(j=1; j<=O[oj-1][0]; j++){
                if(O[oi-1][i]==O[oj-1][j]){
                    flg++;
                }
            }
        }
        if(flg!=O[oi-1][0]){
            mfg = 0;
        }
    }
    
    return mfg;
}