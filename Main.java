/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pkg1005108_offline_2;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;



class SCP
{
    int MAX =107;
    double INF =1000000007;
    double EPS= (1e-12);
    int N,M;
    int W[]=new int[100];
    int K[]=new int[100];
    int E[][]=new int[100][100];
    double DP[][]=new double[10000][100];
    boolean output[][]=new boolean[10000][100];
    int mask[]=new int[100];
    double b[]=new double[100];
    
    
    double Approximation()
    {
        double A[][]=new double[MAX+7][MAX+7];
        double Ret[]=new double[1];
        double max=0.0;
        double cost=0.0;
        int f[]=new int[N];
        int i,j;
        for(i=0;i<M;i++)
            for(j=0;j<K[i];j++)
                A[E[i][j]][i]=-1;
        for(i=0;i<N;i++) A[i][M]=-1;
        for(i=0;i<M;i++) A[N][i]=-W[i];
        A[N][M+1]=0;
        Simplex(A,b,Ret);
        for(i=0;i<N;i++)
        {
            f[i]=0;
        }
        for(i=0;i<M;i++)
        {
            for(j=0;j<K[i];j++)
            {
                f[E[i][j]]++;
            }
        }
        for(i=0;i<N;i++)
        {
            if(f[i] > max) max =f[i];
        }
        for(i=0;i<M;i++)
        {
            if(b[i]>=1/max)
            {
                b[i]=1;
                cost = cost+W[i];
            }
            else b[i]=0;
        }
        return cost;
    }
    
    double Dynamic(int coveredMask, int nowConsiderIndex)
    {
       
        if (nowConsiderIndex == M & coveredMask!=(1<<N) - 1) DP[coveredMask][nowConsiderIndex] = INF;
        else if (nowConsiderIndex == M & coveredMask==(1<<N) - 1) DP[coveredMask][nowConsiderIndex] = 0;
        else
        {
            double valChoose= Dynamic( coveredMask | mask[nowConsiderIndex], nowConsiderIndex+1 ) + W[nowConsiderIndex];
            double valNotChoose= Dynamic( coveredMask, nowConsiderIndex+1 );
            if(valChoose > valNotChoose) 
            {
                DP[coveredMask][nowConsiderIndex]= valNotChoose;
                output[coveredMask][nowConsiderIndex]= false;
            }
            else 
            {
                DP[coveredMask][nowConsiderIndex]= valChoose;
                output[coveredMask][nowConsiderIndex]= true;
            }
        }
        return DP[coveredMask][nowConsiderIndex];
    }

    void Pivot(double A[][],int B[],int P[],int r,int c )
    {
        int i,j;
        int temp= P[c];
        P[c]= (int) B[r];
        B[r]=temp;
        A[r][c] = 1/A[r][c];
        for( j=0;j<=M;j++ ) if( j!=c ) A[r][j] *= A[r][c];
        for( i=0;i<=N;i++ )
        {
            if( i!=r )
            {
                for( j=0;j<=M;j++ ) if( j!=c ) A[i][j] -= A[i][c]*A[r][j];
                A[i][c] = -A[i][c]*A[r][c];
            }
        }
    }
 
    double Feasible(double A[][],int B[],int P[] )
    {
        int r = 0,c = 0,i;
        double p,v;
        while( true )
        {
            for( p=INF,i=0;i<N;i++ ) 
            {
                if( A[i][M] < p )
                {    
                    p = A[i][M];
                    r=i;
                }
            }
            if( p > -EPS ) return 1;
            for( p=0,i=0;i<M;i++ ) 
            {
                if( A[r][i]<p ) 
                {
                    p = A[r][i];
                    c=i;
                }
            }
            if( p > -EPS ) return 0;
            p = A[r][M]/A[r][c];
            for( i=r+1;i<N;i++ )
            {
                if( A[i][c] > EPS )
                {
                    v = A[i][M]/A[i][c];
                    if( v<p ) 
                    {
                        r=i;
                        p=v;
                    }
                }
            }
            Pivot(A,B,P,r,c );
        }
    }

    double Simplex( double A[][],double b[],double Ret[] )
    {
        int B[]=new int[MAX+7];
        int P[]=new int[MAX+7];
        int r = 0,c = 0,i;
        double p,v;
        for( i=0;i<M;i++ ) P[i] = i;
        for( i=0;i<N;i++ ) B[i] = M+i;
        if( Feasible(A,B,P )==0 ) return 0;
        while( true )
        {
            for( p=0,i=0;i<M;i++ ) if( A[N][i] > p ) p = A[N][c=i];
            if( p<EPS ){
                for( i=0;i<M;i++ ) 
                {
                    if( P[i]<M )
                        b[P[i]] = 0;
                }
                for( i=0;i<N;i++ ) 
                {
                    if( B[i]<M ) 
                        b[B[i]] = A[i][M];
                }
                Ret[0] = -A[N][M];
                return 1;
            }
            for( p=INF,i=0;i<N;i++ ){
                if( A[i][c] > EPS ){
                   v = A[i][M]/A[i][c];
                    if( v<p ) 
                    {
                        p = v;
                        r = i;
                    }
                }
            }
            if( p==INF ) return -1;
            Pivot(A,B,P,r,c );
        }
    }

}


public class Main {

    public static Scanner in;

    public static void main(String[] args) throws IOException {
        SCP S;
        in = new Scanner(new File("1005108.txt"));
        PrintWriter pw1 = new PrintWriter(new FileWriter("1005108_1.csv"));
        PrintWriter pw2 = new PrintWriter(new FileWriter("1005108_2.csv"));
        PrintWriter pw3 = new PrintWriter(new FileWriter("1005108_3.csv"));
        pw1.print("Input Size");
        pw1.print(",");
        pw1.print("Execution Time");
        pw1.print(",");
        pw1.println();
        pw2.print("Input Size");
        pw2.print(",");
        pw2.print("Execution Time");
        pw2.print(",");
        pw2.println();
        pw3.print("Approximation");
        pw3.print(",");
        pw3.print("Dynamic Programming");
        pw3.print(",");
        pw3.println();
        int T=in.nextInt();
        for(int i=0;i<T;i++)
        {
            S=new SCP();
            S.N=in.nextInt();
            S.M=in.nextInt();
            System.out.println("Input Size : "+S.N);
            System.out.println("Subset Number : "+S.M);
            pw1.print(S.N);
            pw1.print(",");
            pw2.print(S.N);
            pw2.print(",");
            System.out.println();
            for(int j=0;j<S.M;j++)
            {
                S.W[j]=in.nextInt();
                S.K[j]=in.nextInt();
                for(int k=0;k<S.K[j];k++)
                {
                    S.E[j][k]=in.nextInt();
                }
            }
            long startTime = System.nanoTime();
            System.out.println("----------Approximation Algorithm--------");
            double cost1=S.Approximation();
            System.out.print(" Set  : ");
            for(int h=0;h<S.M;h++) if(S.b[h]==1) System.out.print(h+"  ");
            System.out.println();
            System.out.println(" Cost : "+cost1);
            long endTime = System.nanoTime();
            long duration = ((endTime - startTime)/1000);
            System.out.println(" execution time : " + duration );
            pw1.print(duration);
            pw1.print(",");
            pw1.println();
            System.out.println();
            startTime = System.nanoTime();
            System.out.println("------------Dynamic Algorithm------------");
            for(int x=0;x<S.M;x++)
            {
                S.mask[x]=0;
                for(int y=0;y<S.K[x];y++)
                {
                    S.mask[x]=  S.mask[x] | (1<<S.E[x][y]);
                }
            }
            int coveredMask = (0 << S.M);
            double cost2=S.Dynamic(coveredMask,0);
            System.out.print(" Set  : ");
            for(int f=0;f<((1<<S.N) - 1);f++)
                for(int g=0;g<S.M;g++)
                    if(S.output[f][g]==true) 
                    {
                        System.out.print(g+"  ");
                        f=f | S.mask[g];
                    }
            System.out.println();
            System.out.println(" Cost : "+ cost2);
            endTime = System.nanoTime();
            duration = ((endTime - startTime)/1000);
            System.out.println(" execution time : " + duration );
            pw2.print(duration);
            pw2.print(",");
            pw2.println();
            System.out.println("-----------------------------------------");
            System.out.println();
            pw3.print(cost1);
            pw3.print(",");
            pw3.print(cost2);
            pw3.print(",");
            pw3.println();
            System.out.println();
        }
          pw1.flush();
          pw1.close();
          pw2.flush();
          pw2.close();
          pw3.flush();
          pw3.close();
    }
    
}
