
#include "ComputeSystSimple.C"

class ComputeSystWeight: public ComputeSystSimple {

	public:
	ComputeSystWeight(std::string fileName="SysDatabase.txt"):ComputeSystSimple(fileName){};
	~ComputeSystWeight(){};
	void SetCut(float c){cut=c;};
	void SetDelta(float d){delta=d;};
	float GetWeight(float pt,float rho,float l);
	float GetWeightPlus(float pt,float rho,float l);
	float GetWeightMinus(float pt,float rho,float l);
	private:
	float cut;
	float delta;
};


float ComputeSystWeight::GetWeight(float pt,float rho,float l){
	if(l>cut+delta)return 1;
	else if(l<cut-delta)return 0;
	else return (1./(2*delta)*(l-cut +delta) );
	}
float ComputeSystWeight::GetWeightPlus(float pt,float rho,float l){
	float a=GetSystSimple(pt,rho,l);	
	float l1;

	float x=-999;
	bool findX=false;
	
	for(std::map< range,histo>::iterator it=eff.begin();it!=eff.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
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
				float b=(x+a>1)?1:x+a;
				for(int i=0;i<N;i++)if(  b <= it->second.values.at(i))r=i;
				for(int i=N;i>=0;i++)if( b <= it->second.values.at(i))p=i;
				int k=-1;
				if(r==N-1)k=p;
				if(p==0)k=N;
				
				if(k>=0)	
					l1=lmin + (lmax-lmin)/N*k;
				else
					fprintf(stderr,"error k<0\n");
			
				}
			else {fprintf(stderr,"Error in the l value - check your config file [lmin,lmax)\n");}
			}
	}
	return l1;
}
float ComputeSystWeight::GetWeightMinus(float pt,float rho,float l){
	float a=GetSystSimple(pt,rho,l);	
	float l1;

	float x=-999;
	bool findX=false;
	
	for(std::map< range,histo>::iterator it=eff.begin();it!=eff.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
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
				float b=(x-a<0)?0:x-a;
				for(int i=0;i<N;i++)if(  b >= it->second.values.at(i))r=i;
				for(int i=N;i>=0;i++)if( b >= it->second.values.at(i))p=i;
				int k=-1;
				if(r==N-1)k=p;
				if(p==0)k=N;
				
				if(k>=0)	
					l1=lmin + (lmax-lmin)/N*k;
				else
					fprintf(stderr,"error k<0\n");
			
				}
			else {fprintf(stderr,"Error in the l value - check your config file [lmin,lmax)\n");}
			}
	}
	return l1;
}
