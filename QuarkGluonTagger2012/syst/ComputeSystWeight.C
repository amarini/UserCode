

#ifndef COMP_SYST_SIMPLE_H
	#include "ComputeSystSimple.C"
#endif

#include <cstdio>

#ifndef COMP_SYST_WEIGHT_H
#define COMP_SYST_WEIGHT_H

class ComputeSystWeight: public ComputeSystSimple {

	public:
	ComputeSystWeight(std::string fileName="SystDatabase.txt"):ComputeSystSimple(fileName){SetDebug(0);};
	~ComputeSystWeight(){};
	void SetCut(float c){cut=c;};
	void SetDelta(float d){delta=d;};
	float GetWeight(float pt,float rho,float l);
	float GetSystPlus(float pt,float rho,const char eta,float l);
	float GetSystMinus(float pt,float rho,const char eta,float l);
	int SetDebug(int d){debug_=d;};
	private:
	float cut;
	float delta;
	int debug_;
};


float ComputeSystWeight::GetWeight(float pt,float rho,float l){
	if(l>cut+delta)return 1;
	else if(l<cut-delta)return 0;
	else return (1./(2*delta)*(l-cut +delta) );
	}
float ComputeSystWeight::GetSystPlus(float pt,float rho,const char eta,float l){
	float a=GetSystSimple(pt,rho,eta,l);	
	float l1=-1;

	float x=-999;
	bool findX=false;
	
	if(debug_>0) printf("a=%f\n",a);	

	for(std::map< range,histo>::iterator it=eff.begin();it!=eff.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
		if(eta==it->first.eta)
			{
			int N=it->second.N;
			float lmin=it->second.lmin;
			float lmax=it->second.lmax;
			int n = floor( (l-lmin)*N/(lmax-lmin) );
			if(debug_>0)printf("Readed values %d %.1f %.1f n=%d \n",N,lmin,lmax,n);
			if (n<int(it->second.values.size()))
				{
				x=it->second.values.at(n);
				findX=true;
				int r=-1;
				int p=-1;
				float b=(x+a>1)?1.0-0.00001:x+a;
				if(debug_>1)printf("size=%d >= N=%d\n",int(it->second.values.size()),N);
				for(int i=0;i<N;i++)if(  b <= it->second.values.at(i))r=i;
				for(int i=N-1;i>=0;i--)if( b <= it->second.values.at(i))p=i;
				int k=-1;
				if(r==N-1)k=p;
				if(p==0)k=r;
				if((p==-1)&&(r==-1))k=( it->second.values.at(N-1)< it->second.values.at(0))?N-1:0;
				if((p==0)&&(r==N-1))k=( it->second.values.at(N-1)< it->second.values.at(0))?N-1:0;//
			
				if(debug_>1)printf("r=%d p=%d k=%d\n",r,p,k);	
				if(k>=0)	
					l1=lmin + (lmax-lmin)/N*(k+.5);
				else
					fprintf(stderr,"error k<0 p=%d r=%d N=%d\n",p,r,N);
				//fprintf(stderr,"Plus p=%d r=%d N=%d n=%d k=%d a=%.4f eta=%c\n",p,r,N,n,k,a,eta); //DEBUG
			
				}
			else {fprintf(stderr,"Error in the l value - check your config file [lmin,lmax)\n");}
			}
	}
	return l1;
}
float ComputeSystWeight::GetSystMinus(float pt,float rho,const char eta,float l){
	float a=GetSystSimple(pt,rho,eta,l);	
	float l1=-1;

	float x=-999;
	bool findX=false;
	
	for(std::map< range,histo>::iterator it=eff.begin();it!=eff.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
		if(eta==it->first.eta)
			{
			int N=it->second.N;
			float lmin=it->second.lmin;
			float lmax=it->second.lmax;
			int n = floor( (l-lmin)*N/(lmax-lmin) );
				//fprintf(stderr,"Syst Exact=%.2f ==%.2f\n",l,lmin + (lmax-lmin)/N*(n+.5) );//DEBUG
			if(debug_>0)printf("Readed values %d %.1f %.1f n=%d \n",N,lmin,lmax,n);
			if (n<int(it->second.values.size()))
				{
				x=it->second.values.at(n);
				findX=true;
				int r=-1;
				int p=-1;
				float b=(x-a<0)?0:x-a;
				for(int i=0;i<N;i++)if(  b <= it->second.values.at(i))r=i;
				for(int i=N-1;i>=0;i--)if( b <= it->second.values.at(i))p=i;
				int k=-1;
				if(r==N-1)k=p;
				if(p==0)k=r;
				//for construction the 0 bin can not be there :( - alternatively extend binning including by force 1.0 -> 0.0
				if((p==-1)&&(r==-1))k=( it->second.values.at(N-1)< it->second.values.at(0))?N-1:0;
				if((p==0)&&(r==N-1))k=( it->second.values.at(N-1)< it->second.values.at(0))?N-1:0;//
				
				if(k>=0)	
					l1=lmin + (lmax-lmin)/N*(k+.5);
				else
					fprintf(stderr,"error k<0 p=%d r=%d N=%d\n",p,r,N);
				//fprintf(stderr,"Minus p=%d r=%d N=%d n=%d k=%d a=%.4f eta=%c\n",p,r,N,n,k,a,eta); //DEBUG
			
				}
			else {fprintf(stderr,"Error in the l value - check your config file [lmin,lmax)\n");}
			}
	}
	return l1;
}
#endif
// COMP_SYST_WEIGHT_H
