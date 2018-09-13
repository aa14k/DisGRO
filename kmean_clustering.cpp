#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>

using namespace std;

void	k_mean_clustering(vector<vector<double> >& matrix,int K,vector<int>& centerMin,vector<int>& idclusterMin);

static	int	myrandom(int max);
static	inline double myabs(double x) { return x>0? x:-x;}


void	k_mean_clustering(vector<vector<double> >& matrix,int K,vector<int>& centerMin,vector<int>& idclusterMin)
{
  /*
    Cluster the conformations into "K" clusters, number from 0 to K-1. 

    INPUT:
      vector<vector<double> >& matrix	--> the NC x NC distance matrix, the smaller the matrix[i][j], the similar they are
      int K --> K value, number of clusters
    OUTPUT:
      vector<int>& centerMin  --> stores the index of the centers of each cluster, the size of centerMin is K
        for examplke, if centerMin[5]=1234, it meas that the center of the 5-th center is at the 1234th-conformation
      vector<int>& idclusterMin --> stores the cluster index which the conformation belongs, the size of 
        idclusterMin is NC, where NC is the number of conformations.
        for example, if idclusterMin[1234]=5, it means that the 1234-th conformation belongs to the 5-th cluster

    Algorithm detail:
        Since the results of K-mean cluster depends on the intial conditions, we run the same k-mean clustering
    method "ntrial" times, and find the best according to some criterion.
  */

  //the same K-mean clustering method is run for ntrial times, and output the best one
  int	ntrial=100, NC=matrix.size();
  
  centerMin.clear(); idclusterMin.clear();
  for(int i=0;i<K;i++) centerMin.push_back(i);
  for(int i=0;i<NC;i++) idclusterMin.push_back(i); //simply to allocate memory

  int	idcluster[NC];
  int	center[K]; //store the K centers
  int	NNN[K]; //store the number of conformations of each cluster
  int	NNN2[K]; //tmp memory

  double	ERRmin=1.0e30;

  for(int it=0;it<ntrial;it++) {
    //cout<<"--------------trail "<<it<<" ------------------------------------------"<<endl;
    //find K centers randomly
    for(int i=0;i<K;i++) {
      center[i]=myrandom(NC);
      int j;
      do { //one conformation cannot be the center of two different clusters
        for(j=0;j<i;j++) if(center[i]==center[j]) break;
        if(j<i) center[i]=myrandom(NC);
      } while(j<i);
    }
  //cout <<"K		"<<K <<endl;    
    double	err0, err=1.e30;
    do {

      //scan the conformations and find which cluster it belongs to
      for(int i=0;i<NC;i++) {
        int imin=0;
        double min=matrix[ center[0] ][ i ],tmp;
        for(int k=1;k<K;k++) {
          tmp=matrix[ center[k] ][ i ];
          if(min>tmp) { min=tmp; imin=k; }
        }
        //ok, conformation-i is closet to the cluster center-imin
        idcluster[i]=imin;
      }
      //calc the new K-center
      err0=err;
      err=0;
      for(int k=0;k<K;k++) {
        //where is the new center of cluster-k? 
        //image the M conformations in cluster-k forms a MxM matrix,
        //we calculate the sum of each column and find the minimum one,
        //which corresponds to the center
        double min=1.0e30;
        int	imin=0;
        for(int i=0;i<NC;i++) {
          if(idcluster[i]!=k) continue;
          double sum=0.0;
          for(int j=0;j<NC;j++) {
            if(idcluster[j]!=k)	continue;
            sum+=matrix[j][i];
          }
          if(min>sum) {min=sum;imin=i;}
        }
        center[k]=imin; //this is the new center of cluster-k
        err+=min;
      }
    }while(myabs(err-err0)>1.0e-3);

    if(ERRmin>err) {
      ERRmin=err;
      for(int i=0;i<K;i++) centerMin[i]=center[i];
      for(int i=0;i<NC;i++) idclusterMin[i]=idcluster[i];
    }

  }
  for(int i=0;i<K;i++) {
    int ic=0;
    for(int j=0;j<NC;j++) if(idclusterMin[j]==i) ic++;
    NNN[i]=ic;
  }
  for(int i=0;i<K;i++) {
    int nmax=-1, jmax;
    for(int j=0;j<K;j++) if(nmax<NNN[j]) {nmax=NNN[j];jmax=j;}
    center[i]=centerMin[jmax];
    for(int j=0;j<NC;j++) if(idclusterMin[j]==jmax) idcluster[j]=i;
    NNN2[i]=NNN[jmax];
    NNN[jmax]=-1;
  }
  
  for(int i=0;i<K;i++) {centerMin[i]=center[i]; NNN[i]=NNN2[i];}
  for(int i=0;i<NC;i++) idclusterMin[i]=idcluster[i];

  for(int i=0;i<K;i++) {
    cout<<"cluster No."<<i<<", Center at "<<centerMin[i]<<", N# "<<NNN[i]; cout<<endl;
  }
  
}


int	myrandom(int max)
{
  static int taginit=0;
  if(taginit==0) {
    srand48(time(NULL));
    taginit=1;
  }
  
  return int(drand48()*max);
}

