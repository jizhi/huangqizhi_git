/*
ReadMeanPAON4vismtx.cc is the c++ code that I wrote to instead of rdvisip4 to read and merge the .ppf files, get the data that we want.

Using command:
#	./ReadMeanPAON4vismtx bao5path bao6path t1,t2,dt v1,v2 f1,f2,df outname
	./ReadMeanPAON4vismtx bao5path bao6path t1,t2,dt vis f1,f2,df outname

bao5path, bao6path:
	Path of the data in bao5 and bao6 directories.
	For example, bao5path=/Raid-bao5/PAON4/CygA2mar15, bao6path=/Raid-bao6/PAON4/CygA2mar15.

t1,t2,dt
	t1,t2: select the vismtx.ppf files from t1 to t2 (each files 1 second).
	dt: average the data every dt files. Generally we set dt=10.

vis:
	string. vis=all: save all visibilities. vis=ac: save effective visibilities, 8 auto, 12 cross (include 3-4).

f1,f2,df:
	f1,f2: select the data from frequency f1 to f2 (total range 0-4095).
	df: average the data every df frequency bins. Generally we set df=16.

outname:
	Set output file's name.
	(1) outname=CygA.ppf, will save a ppf file to current directory.
	(2) outname=CygA.fits, will save a fits file.
	(3) outname=CygA (without format), will save both ppf and fits.
	(4) outname=../data/CygA.ppf, will save the file to other directory "../data"
*/


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <list>
#include <vector>
#include <complex>
#include <dirent.h>
#include <typeinfo>
#include <fstream>

// Sophya
#include "array.h"
#include "tmatrix.h"
#include "ntuple.h"
#include "ppersist.h"
#include "fitsioserver.h"
#include "fiosinit.h"
#include "swfitsdtable.h"

using namespace std;
using namespace SOPHYA;


//////////////////////////////////////////////////
string StringBlock(string s, int i=0, int j=0) {
	string block;
	if (i<0) i=s.size()+i;
	if (j<0) j=s.size()+j;
	if (j==0) j=s.size()-1;
	block.append(s, i, j-i+1);
	return block;
}

template<class T1, class T2>
int strcmp(T1 a, T2 b) {
	string sa, sb;
	sa = a;
	sb = b;
	return strcmp(sa.c_str(), sb.c_str());
}

string itoa(int n) {
	char buffer[30];
	string fmt = "%d";
	sprintf(buffer, fmt.c_str(), n);
	string valuestr = buffer;
	return valuestr; 
}

template<class T>
TArray<int> Shape(TArray<T> array) {
	int rank = array.Rank();
	TArray<int> shape(rank);
	for (int i=0; i<rank; i++)
		shape(i) = array.Size(i);
	return shape;
}

void Pause() {
	cin.get();
}

int LS(string path) {
	int num=0;
	DIR *dip;
	struct dirent *dit;
	dip = opendir( path.c_str() );
	if (dip == NULL) {
		cout << "Error: No directory \""+path+"\"." << endl;
		exit(0);
	}
	while (dip != NULL) {
		dit = readdir(dip);
		if (dit == NULL) break;
		num++;
	}
	closedir(dip);
	return num-9;
}

TArray< complex<r_4> > ReorderFreq(TArray< complex<r_4> > vismtx) {
	TArray< complex<r_4> > vismtxre(4096,36);
	// Output 1
	for (int i=0; i<2047; i++) {
		vismtxre(Range(i+2049),Range::all(),Range::all()) = vismtx(Range(i),Range::all(),Range::all()); 
	}
	// Frequency in the middle
	vismtxre(Range(2048),Range::all(),Range::all()) = vismtx(Range(2047),Range::all(),Range::all());
	// Output 2
	for (int i=0; i<2047; i++)  {
		vismtxre(Range(2047-i),Range::all(),Range::all()) = vismtx(Range(2048+i),Range::all(),Range::all());
	}
	// Z(N/2).real -> Z(0).image
	vismtxre(Range(0),Range::all(),Range::all()) = vismtx(Range(4095),Range::all(),Range::all());
	return vismtxre;
}

template<class T>
void Saveppf(TArray<T> array, string outname) {
	POutPersist pout(outname);
	TArray<int> shape = Shape(array);
	if (shape[1] == 36) pout << array;
//	else {
//		pout << PPFNameTag("Auto") << (real(array))(Range::all(),Range(0,7),Range::all());
//		pout << PPFNameTag("Cross") << array(Range::all(),Range(8,shape[1]-1),Range::all());
//	}
	cout << outname+"  -->  Saved." << endl;
}

void Savefits(TArray< complex<r_4> > array, string outname, vector< vector<string> > key, vector< vector<MuTyV> > keyvalue, vector< vector<string> > comment) {
	/*
	If don't save the FITS, one reason is the length of keywords out of the range.
	Note that key.size() must =3
	*/
	FitsIOServerInit();
	FitsInOutFile fos("!"+outname, FitsInOutFile::Fits_Create);
	FitsArrayHandler<r_4> fah; 
	TArray<int> shape = Shape(array);
	sa_size_t offset[shape.Size()];
	long long size[shape.Size()];
	vector< TArray<r_4> > vecarray;
	r_4 zo = 0;
	int nv = shape[1];
	if (nv == 36) {
		vecarray.push_back(real(array));
		vecarray.push_back(imag(array));
	}
	else {
		vecarray.push_back( (real(array))(Range::all(),Range(0,7),Range::all())+zo );
		vecarray.push_back( (real(array))(Range::all(),Range(8,nv-1),Range::all())+zo );
		vecarray.push_back( (imag(array))(Range::all(),Range(8,nv-1),Range::all())+zo );
		array.ZeroSize();
	}
	for (int i=0; i<int(vecarray.size()); i++) {
		int n = i;
		if (vecarray.size() == 2) n = n + 1;
		shape = Shape(vecarray[i]);
		for (int j=0; j<shape.Size(); j++) {
			offset[j] = 0;
			size[j] = shape[j];
		}
		fos.CreateImageHDU(FitsTypes::ImageType(vecarray[i](0)), shape.Size(), size);
		for (int j=0; j<int(key[n].size()); j++)
			fos.WriteKey((key[n])[j], (keyvalue[n])[j], (comment[n])[j]);
		fah = vecarray[i];
		fah.WriteAtOffset(fos, offset);
	}
	cout << outname+"  -->  Saved." << endl;
}
//////////////////////////////////////////////////



int main(int narg, char *arg[]) {
cout << "--------------------------------------------------" << endl;
cout << "--------------------------------------------------" << endl << endl;

// Reorder the frequency or not?
bool refreq = true;

// Read the input parameters
string path5=arg[1], path6=arg[2];
if (strcmp(path5[path5.size()-1],'/')!=0) path5 = path5+"/";
if (strcmp(path6[path6.size()-1],'/')!=0) path6 = path6+"/";
int n=0;
for (n=path5.size()-2; n>0; n--)
	if (strcmp(path5[n],'/')==0) break;
string dayname = StringBlock(path5, n+1, -2);

// Count number of files
int Num = LS(path5);

EnumeratedSequence autocrossvis;
// If we want to change the "ac", just need to modify autocrossvis !!!
//autocrossvis = 0,8,15,21,26,30,33,35, 1,2,3,9,10,27,28,29,31,32;
autocrossvis = 0,8,15,21,26,30,33,35, 1,2,3,9,10,16,27,28,29,31,32,34;

int t1,t2,dt, v1=0,v2, f1,f2,df;
TArray<int> acv(autocrossvis.Size());
acv = autocrossvis;
sscanf(arg[3], "%d,%d,%d", &t1,&t2,&dt);
//sscanf(arg[4], "%d,%d", &v1,&v2);
string vis=arg[4];
if (vis == "all") v2 = 35;
else if (vis == "ac") v2 = acv.Size()-1;
else {
	cout << "Error: ./ReadMeanPAON4vismtx bao5path bao6path t1,t2,dt vis f1,f2,df outname, vis is not 'all' nor 'ac'." << endl;
	exit(0);
}
sscanf(arg[5], "%d,%d,%d", &f1,&f2,&df);
// Because of the averaging dt and df, t2 and f2 may be changed
if (t2 > Num/4-1) t2 = Num/4-1;
int Nt=(t2-t1+1)/dt, Nf=(f2-f1+1)/df, Nv=(v2-v1+1);
t2 = t1+dt*Nt-1; f2 = f1+df*Nf-1;

string outdirname=arg[6];
// devide into outdir and outname
n=-1;
for (int i=outdirname.size()-1; i>=0; i--) {
	if (strcmp(outdirname[i],'/') == 0) {
		n = i; break;
} }
string outdir, outname = StringBlock(outdirname, n+1, -1);
if (n == -1) outdir = "./";
else outdir = StringBlock(outdirname, 0, n); // with "/" end
// save file format, 0->ppf, 1->fits, 2->both
int savefmt=2;
if (strcmp(StringBlock(outname,-5,-1),".fits")==0) savefmt=1;
else if (strcmp(StringBlock(outname,-4,-1),".ppf")==0) savefmt=0;
// mkdir outdir
DIR *dip;
dip = opendir(outdir.c_str());
n = 0;
for (int i=0; i<int(outdir.size()); i++) 
	if (strcmp(outdir[i],'.')!=0 && strcmp(outdir[i],'/')!=0) n = 1;
if (n==1 && dip==NULL) system(("mkdir "+outdir).c_str());


// file name vismtx_1/2/3/4_number.ppf
// Result array
TArray< complex<r_4> > array(Nf, Nv, Nt);
DVList info;
string dateobs;
double deltime=0;
uint_8 firstFC=0, firstTT=0, lastFC=0, lastTT=0, npaqsum=0;
for (int i=0; i<Nt; i++) {
	n = t1 + i*dt;
	TArray< complex<r_4> > vismtxmean(4096,36), vismtxv(f2-f1+1,Nv);//, vismtxtvf(Nf,Nv);

	for (int k=0; k<dt; k++) {
		TArray< complex<r_4> > vismtx(4096,36);
		for (int j=0; j<4; j++) {
			TArray< complex<r_4> > vismtx5, vismtx6;
			string ppfname = "vismtx_"+itoa(j)+"_"+itoa(n+k)+".ppf";
			string ppfname5=path5+ppfname, ppfname6=path6+ppfname;
			PInPersist pin5(ppfname5);
			pin5 >> vismtx5; // shape=(2048,9)
			// Get .Info()
			if (j == 0) {
				info = vismtx5.Info();
				deltime = deltime + info.GetD("DELTIME");
				if (i==0 && k==0) { 
					dateobs = info.GetS("DATEOBS");
					firstFC = info.GetI("FirstFC");
					firstTT = info.GetI("FirstTT");
					npaqsum = info.GetI("NPAQSUM");
				}
				if (i==Nt-1 && k==dt-1) { 
					lastFC = info.GetI("LastFC");
					lastTT = info.GetI("LastTT");
				}
			}
			PInPersist pin6(ppfname6);
			pin6 >> vismtx6;
			vismtx(Range(0,2047),Range(9*j,9*(j+1)-1),Range::all()) = vismtx5;
			vismtx(Range(2048,4095),Range(9*j,9*(j+1)-1),Range::all()) = vismtx6;
		} // end for j
		vismtxmean = vismtxmean + vismtx;
		cout << "\r" << dayname << "  Completed  " << (i*dt+k+1) << " / " << (t2-t1+1);
		if (i*dt+k+1 == t2-t1+1) cout << endl;
		cout << flush;
		
	} // end for k

	// Average the time
	vismtxmean = vismtxmean / complex<r_4>(dt,0);

	// Reorder frequency
	if (refreq) vismtxmean = ReorderFreq(vismtxmean);

	// Select visibility and frequency
	if (Nv == 36) vismtxv = vismtxmean(Range(f1,f2),Range::all(),Range::all());
	else {
		for (int j=0; j<acv.Size(); j++) {
			vismtxv(Range::all(),j,Range::all()) = vismtxmean(Range(f1,f2),acv[j],Range::all());
	} }
	vismtxmean.ZeroSize();
	
	TArray< complex<r_4> > vismtxvf(Nf, Nv);
	// Average frequency 
	if (df > 1) { 
		for (int j=0; j<Nf; j++) {
			for (int k=0; k<df; k++) {
				vismtxvf(Range(j),Range::all(),Range::all()) = vismtxvf(Range(j),Range::all(),Range::all()) + vismtxv(Range(k+j*df),Range::all(),Range::all());
		} } 
		vismtxvf = vismtxvf / complex<r_4>(df,0);
	} // end if
	else vismtxvf = vismtxv;
	array(Range::all(),Range::all(),i) = vismtxvf;

} // end for i


// Write .Info()
array.Info().SetS("DATEOBS", dateobs);
array.Info().SetD("DELTIME", deltime);
array.Info().SetI("FirstFC", firstFC);
array.Info().SetI("FirstTT", firstTT);
array.Info().SetI("LastFC",  lastFC);
array.Info().SetI("LastTT",  lastTT);
array.Info().SetD("TotFreq",  250.);
array.Info().SetI("TotFBin",  4096);
array.Info().SetI("NPAQSUM", npaqsum);
array.Info().SetI("TimeBin1", t1);
array.Info().SetI("TimeBin2", t2);
array.Info().SetI("TimeAve", dt);
array.Info().SetI("FreqBin1", f1);
array.Info().SetI("FreqBin2", f2);
array.Info().SetI("FreqAve", df);
array.Info().SetS("AutoVis", "1H, 2H, 3H, 4H, 1V, 2V, 3V, 4V");
array.Info().SetS("CrossVis", "12H, 13H, 14H, 23H, 24H, 12V, 13V, 14V, 23V, 24V");

// FITS key
vector< vector<string> > key, comment;
vector<string> key1, com1;
vector< vector<MuTyV> > value;
vector<MuTyV> value1;
MuTyV mu;

key1.push_back("DATEOBS");
key1.push_back("DELTIME");
key1.push_back("FirstFC");
key1.push_back("FirstTT");
key1.push_back("LastFC");
key1.push_back("LastTT");
key1.push_back("TotFreq");
key1.push_back("TotFBin");
key1.push_back("NPAQSUM");
key1.push_back("TimeBin1");
key1.push_back("TimeBin2");
key1.push_back("TimeAve");
key1.push_back("FreqBin1");
key1.push_back("FreqBin2");
key1.push_back("FreqAve");
key.push_back(key1); // key.size() == 3
key.push_back(key1);
key.push_back(key1);
(key[0]).push_back("Auto");
(key[1]).push_back("CrosReal");
(key[2]).push_back("CrosImag");


com1.push_back("");
com1.push_back("Sum all TimeBins");
com1.push_back("");
com1.push_back("");
com1.push_back("");
com1.push_back("");
com1.push_back("Total band width of the frequency (MHz)");
com1.push_back("Total number of frequency bins");
com1.push_back("NPAQSUM of original vismtx.ppf");
com1.push_back("First TimeBin");
com1.push_back("Last TimeBin");
com1.push_back("Averaged over time per");
com1.push_back("First FreqBin");
com1.push_back("Last FreqBin");
com1.push_back("Averaged over frequency per");
comment.push_back(com1);
comment.push_back(com1);
comment.push_back(com1);
(comment[0]).push_back("");
(comment[1]).push_back("");
(comment[2]).push_back("");

mu = dateobs;
value1.push_back(mu);
mu = deltime;
value1.push_back(mu);
mu = firstFC;
value1.push_back(mu);
mu = firstTT;
value1.push_back(mu);
mu = lastFC;
value1.push_back(mu);
mu = lastTT;
value1.push_back(mu);
mu = 250;
value1.push_back(mu);
mu = 4096;
value1.push_back(mu);
mu = npaqsum;
value1.push_back(mu);
mu = t1;
value1.push_back(mu);
mu = t2;
value1.push_back(mu);
mu = dt;
value1.push_back(mu);
mu = f1;
value1.push_back(mu);
mu = f2;
value1.push_back(mu);
mu = df;
value1.push_back(mu);
value.push_back(value1);
value.push_back(value1);
value.push_back(value1);
mu = "Auto correlations: 1H, 2H, 3H, 4H, 1V, 2V, 3V, 4V";
(value[0]).push_back(mu);
mu = "Real parts of Cross: 12H, 13H, 14H, 23H, 24H, 34H, 12V, ...";
(value[1]).push_back(mu);
mu = "Imag parts of Cross: 12H, 13H, 14H, 23H, 24H, 34H, 12V, ...";
(value[2]).push_back(mu);

cout<<"Real arguments:  "<<t1<<","<<t2<<","<<dt<<"  "<<vis<<"  "<<f1<<","<<f2<<","<<df<<endl;
if (savefmt == 0) Saveppf(array, outdirname);
else if (savefmt == 1) {
	Savefits(array, outdirname, key, value, comment);
}
else if (savefmt == 2) {
	Saveppf(array, outdirname+".ppf");
	Savefits(array, outdirname+".fits", key, value, comment);
}
cout << "--------------------------------------------------" << endl;
cout << "--------------------------------------------------" << endl;

} // end main()
