// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// fitHRG - fits a hierarchical random graph (hrg) model to data
// Copyright (C) 2005-2008 Aaron Clauset
// Copyright (C) 2014 Leto Peel
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 19 April 2006
// Modified     : 19 May 2007
//			 : 19 May 2008 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// This program computes the consensus dendrogram of the sampled set of HRGs. This program requires
// as input a seed HRG model, which can be produced by the fitHRG program. The MCMC samples
// HRG models with probability proportional to their likelihood of generating the observed data, 
// so, by default, the "temperature" of the sampling done here is T=2.0, which we choose because 
// this value gives relatively good results on many of the test cases we have tried. T=1.0 often
// produces too few splits in the consensus tree because the split statistics do not converge 
// sufficiently well over the random walk. This is likely because the sampled set is too large for
// the sample procedure to work well. An optional parameter can adjust T, in the case that the
// split statistics do not converge (e.g., two independent runs of the sampling produce quite 
// different looking consensus trees).
//
// ****************************************************************************************************
// 
//  See http://www.santafe.edu/~aaronc/randomgraphs/ for more information, updates to the code, etc.
//
// ****************************************************************************************************
// *** PROGRAM USAGE NOTES ****************************************************************************
//
// This program reads three files from disk: a .hrg file that seeds the MCMC, a .pairs file that 
// contains the edge list of your network (see fitHRG for formatting requirements), and a -names.lut file 
// that is produced by fitHRG. If any of these files are formatted incorrectly, the program will crash.
//
// ****************************************************************************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include "stdlib.h"
#include "time.h"
#include "math.h"

#include "MersenneTwister.h"
#include "dendro_eq.h"
#include "graph.h"
#include "rbtree.h"

using namespace std;

// ******** Function Prototypes ***************************************************************************

bool		markovChainMonteCarlo();
string    num2str(const unsigned int);
bool		parseCommandLine(int argc, char * argv[]);
bool		readGraphFile();
bool		readLUTFile();
bool		readSplitsFiles();
void		recordHRG(const unsigned int, const double);
void		recordStats();

// ******** Structures and Constants **********************************************************************

struct ioparameters {
	int			n;				// number vertices in input graph
	int			m;				// number of edges in input graph
	int 			ng;				// number of input graphs
	double		T;				// temperature parameter for MCMC sampling

	string		d_dir;			// working directory
	string		f_dendro;			// (in-file) seed dendrograph structure (from MLGraph)
	vector<string>		f_hists;			// (in-file) list of hierarchical statistics files
	string		f_lut;			// (in-file) names look-up table (from MLGraph)
	string		f_pairs;			// (in-file) original adjacency list of graph (from researcher)
	vector<string>	f_ins;		// list of input .pairs files if using multiple graphs
	string		f_splits;			// (out-file) hierarchical statistics
	string		f_consensus;		// (out-file) name of output consensus tree file
	string		s_scratch;		// scratch space for building filenames
	string		s_tag;			// user defined filename tag
	int 			num_hists;		// number of consensus trees
	int			timer;			// timer for reading input
	bool			flag_timer;		// (flag) for timer
	bool			flag_c;			// (flag) for consensus of consensuses
	bool			flag_f;			// (flag) for consensus of a HRG
	bool			flag_name;			// (flag) for 
	string		start_time;		// time simulation was started
	bool			flag_m;			// flag for multiple input graphs 
};

// ******** Global Variables ******************************************************************************

ioparameters	ioparm;				// program parameters
rbtree		namesLUT;				// look-up table; translates .hrg   indices to .pairs indices
rbtree		namesLUTr;			// look-up table; translates .pairs indices to .hrg   indices
dendro*		d;					// inferred dendrograph data structure
unsigned int	t;					// time step number (max = 2^32; rolls over)
double		bestL;				// best likelihood (this iteration)
int			num_samples;			// number of samples to take at equilibrium
int			sample_num;			// sample number
MTRand		mtr;					// Mersenne Twister random number generator instance
char			pauseme;

// ******** Main Loop *************************************************************************************

int main(int argc, char * argv[]) {
	ioparm.n		= 0;					// DEFAULT VALUES for runtime parameters
	ioparm.T		= 2.0;				// 
	ioparm.timer   = 20;				// 
	ioparm.s_tag   = "s1";				// 
	ioparm.flag_f = ioparm.flag_c = ioparm.flag_name = ioparm.flag_m = false;
	time_t t1      = time(&t1);			// 
	num_samples    = 10000;				// 
	t              = 1;					// 

	if (parseCommandLine(argc, argv)) {
		d = new dendro;
		readLUTFile();								// read -names.lut file
		
		if (ioparm.flag_f) {
			readGraphFile();							// read .pairs file
			d->importDendrogramStructure(ioparm.f_dendro);	// read .hrg file
			d->refreshLikelihood();
			bestL = d->getLikelihood();					// 
			ioparm.start_time = asctime(localtime(&t1));
			cout << "\nstep   \tLogL       \tbest LogL\t% complete\tconsensus size\n";
			while (true) {
				if (!(markovChainMonteCarlo())) { return 0; }
				if (t == 4294967294 or t < 0) { t = 1; }	// rollover	
				if (sample_num >= num_samples) { break; }
			}
			recordStats();
		} else if (ioparm.flag_c) {
			readSplitsFiles();
		}
		cout << "Writing consensus tree to file...\n"; 
		d->recordConsensusTree(ioparm.f_consensus);		// record consensus tree
		return 1;
	} else { return 0; }
}

// ******** Function Definitions **************************************************************************

bool markovChainMonteCarlo() {
	double	dL, Likeli;
	bool		flag_taken;
	double	ptest = 1.0/(50.0*(double)(ioparm.n)); // 1.0/50.0;
	int		thresh = 200*ioparm.n;
	bool		flag_go = false;
	
	// Since we're sampling uniformly at random over the equilibrium walk, we just need to do
	// a bunch of MCMC moves and let the sampling happen on its own. We choose the number of 
	// moves to make to be large so as to reduce the overhead of calling this function many
	// times.
	for (int i=0; i<65536; i++) {
		// Make a single MCMC move
		if (!(d->monteCarloMove(dL, flag_taken, ioparm.T))) { return false; }
		
		// Even though we're just sampling the equilibrium space of dendrograms, we might
		// run across one with a higher log-likelihood than the one we started from. If this
		// happens, write it out to disk.
		Likeli = d->getLikelihood();				// get this likelihood
		if (Likeli > bestL) { bestL = Likeli; }		// store the current best likelihood
		
		// We sample the dendrogram space once every n MCMC moves (on average). Depending on
		// the flags on the command line, we sample different aspects of the dendrograph
		// structure.
		if (t > thresh and mtr.randExc() < ptest) {
			sample_num++;
			if (!d->sampleSplitLikelihoods(sample_num)) { return false; }
			if (sample_num > num_samples) { i = 65536; }
			flag_go = true;
		}
		
		// Write some stuff to standard-out to describe the current state of things.
		if (t % 16384 == 1) {
			cout << "[" << t << "]\t" << Likeli << "   \t(" << bestL << ")\t";
			cout << 100.0*(double)(sample_num)/(double)(num_samples) << " %\t";
			if (t < thresh or t % 32768 != 1 or !flag_go)	{ cout << "  \n"; }
			else										{ cout << d->getConsensusSize() << " splits\n"; }
		}
		t++;									// Increment time step
	}
	
	d->refreshLikelihood();						// correct floating-point errors O(n)
	return true;
}

// ********************************************************************************************************

string num2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	string str  = "";
	if (input == 0) { str = "0"; } else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}

// ********************************************************************************************************

bool parseCommandLine(int argc, char * argv[]) {
	int argct = 1;
	string temp, temp2, ext;
	string::size_type pos, pos2;
	bool safeExit = false;
	bool flag_temp;

	if (argc==1) {
		cout << "\n  -- Hierarchical Random Graphs : Consensus Dendrograms --\n";
		cout << "  by Aaron Clauset (copyright 2006-2008)\n\n";
		cout << "  consensusHRG is a command line program that is intended to be used in\n";
		cout << "  conjunction with the fitHRG program. This program takes a seed .hrg file\n";
		cout << "  and then computes the consensus dendrogram of the HRGs sampled at\n";
		cout << "  equilibrium.\n\n";
		cout << "  Note: consensusHRG takes as argument a .hrg file, which it assumes was\n";
		cout << "  produced by fitHRG, imports a corresponding .pairs file that contains the\n";
		cout << "  edge list and the .lut file produced by fitHRG. These files must be\n";
		cout << "  located in the same directory as the .hrg file.\n\n";
		cout << "  Parameterizations:\n";
		cout << "  -f <file>       Input a seed equilibrium .hrg file\n";
		cout << "  -k <int>        (optional) Number of samples, default = 10000\n";
		cout << "  -T <string>     (optional) Temperature of MCMC, default = 2.0\n";
		cout << "  -t <string>     (optional) Label for this run\n";
		cout << "\n";
		cout << "  ./consensusHRG -f seed.hrg\n";
		cout << "  ./consensusHRG -f seed.hrg -k 100000 -t test\n";
		cout << "  ./consensusHRG -f seed.hrg -t test -T 2.0\n";
		cout << "\n";
		return false;
		
	} else {
		while (argct < argc) {
			temp = argv[argct];
			
			if (temp == "-k") {
				argct++;
				if (argct < argc) {
					num_samples = atoi(argv[argct]);
					if (num_samples == 0) { cout << " Warning: malformed modifier for -k; using default.\n"; argct--; num_samples = 10000; } 
				} else { cout << " Warning: missing modifier for -k argument; using default.\n"; argct--; }
				
			} else if (temp == "-T") {
				argct++;
				if (argct < argc) {
					ioparm.T = atof(argv[argct]);
					if (ioparm.T < 0.0) { cout << " Warning: malformed modifier for -T; using default.\n"; argct--; ioparm.T = 2.0; } 
				} else { cout << " Warning: missing modifier for -T argument; using default.\n"; argct--; }
				
			} else if (temp == "-t") { ioparm.s_tag = argv[++argct];
			} else if (temp == "-f" and !ioparm.flag_c) {
				ioparm.flag_f = true;
				temp2 = argv[++argct];
				ext = ".hrg";
				pos = temp2.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .hrg format.\n"; return safeExit; }
				ioparm.f_dendro = temp2;
				ext = "/";
				pos = string::npos;
				for (int i=0; i < temp2.size(); i++) { if (temp2[i] == '/') { pos = i; } }
				if (pos != string::npos) {
					ioparm.d_dir = temp2.substr(0, pos+1);
					temp2 = temp2.substr(pos+1,temp2.size()-pos-1);
				}
				flag_temp = false;
				for (int i=temp2.size()-1; i > 0; i--) { 
					if (temp2[i] == '_' and  flag_temp) { flag_temp = false; pos2 = i; break; } 
					if (temp2[i] == '_' and !flag_temp) { flag_temp = true;  pos  = i; }
				}
				ioparm.s_scratch    = temp2.substr(0,pos);
				safeExit = true;
			} else if (temp == "-m") {
				ioparm.flag_m = true;				
				argct++;
				ioparm.ng = atoi(argv[argct]);
				cout << "Number of input .pairs files: " << ioparm.ng << endl;
				if (ioparm.ng<=0) { cout << " Error: Requires at least one graph file.\n"; return safeExit; }
				for (int c=0; c<ioparm.ng; c++) {
					
					argct++;
					temp = argv[argct];
					ext = ".pairs";
					pos = temp.find(ext,0);
					if (pos == string::npos) { cout << " Error: Input file must claim to be .pairs format.\n"; return safeExit; }
					ioparm.f_ins.push_back(temp);
				}
				//~ ext = "/";
				//~ pos = string::npos;
				//~ for (int i=0; i < temp.size(); i++) { if (temp[i] == '/') { pos = i; } }
				//~ if (pos != string::npos) {
					//~ ioparm.d_dir = temp.substr(0, pos+1);
					//~ temp = temp.substr(pos+1,temp.size()-pos-1);
				//~ }
				//~ // now grab the filename sans extension for building outputs files
				//~ for (int i=0; i < temp.size(); i++) { if (temp[i] == '.') { pos = i; } }
				//~ ioparm.s_scratch = temp.substr(0,pos);
//~ //				ioparm.f_stat    = ioparm.d_dir + ioparm.s_scratch + "-L.xy";
				safeExit         = true;
			} else if (temp == "-c" and !ioparm.flag_f) {
				ioparm.flag_c = true;
				argct++;
				ioparm.num_hists = atoi(argv[argct]);
				cout << "Number of split files: " << ioparm.num_hists << endl;
				if (ioparm.num_hists<=0) { cout << " Error: Requires at least one splits file.\n"; return safeExit; }
				for (int c=0; c<ioparm.num_hists; c++) {
					argct++;
					if (argct>=argc) { cout << " Error: Not enough splits files specified.\n"; return safeExit; }
					temp2 = argv[argct];
					ext = ".splits";
					pos = temp2.find(ext,0);
					if (pos == string::npos) { cout << " Error: Input file must claim to be .splits format.\n"; return safeExit; }
					ioparm.f_hists.push_back(temp2);
					ext = "/";
					pos = string::npos;
					for (int i=0; i < temp2.size(); i++) { if (temp2[i] == '/') { pos = i; } }
					if (pos != string::npos) {
						ioparm.d_dir = temp2.substr(0, pos+1);
						temp2 = temp2.substr(pos+1,temp2.size()-pos-1);
					}
					flag_temp = false;
					for (int i=temp2.size()-1; i > 0; i--) { 
						if (temp2[i] == '_' and  flag_temp) { flag_temp = false; pos2 = i; break; } 
						if (temp2[i] == '_' and !flag_temp) { flag_temp = true;  pos  = i; }
					}
					ioparm.s_scratch    = temp2.substr(0,pos);
					safeExit = true;
				}
			} else if (temp == "-n") {
				ioparm.flag_name = true;
				argct++;
				temp = argv[argct];
				ext = ".lut";
				pos = temp.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .lut format.\n"; return safeExit; }
				ioparm.f_lut = temp;
				safeExit         = true;
			} else { cout << " Warning: ignored argument " << argct << " : " << temp << endl; }
			argct++;
		}
	}
	// now grab the filename sans extension for building outputs files
	if (ioparm.flag_c) {
		ioparm.f_consensus	= ioparm.d_dir + ioparm.s_scratch + "_" + ioparm.s_tag + "-consensuses.tree";
	}else {
		ioparm.f_consensus	= ioparm.d_dir + ioparm.s_scratch + "_" + ioparm.s_tag + "-consensus.tree";
	}
	ioparm.f_splits	= ioparm.d_dir + ioparm.s_scratch + "_" + ioparm.s_tag + "-consensus.splits";
	if (!flag_temp) { ioparm.s_scratch = temp2.substr(0,pos2);  }
	if (!ioparm.flag_m) { ioparm.f_pairs 	= ioparm.d_dir + ioparm.s_scratch + ".pairs"; }
	if (!ioparm.flag_name) { ioparm.f_lut   	= ioparm.d_dir + ioparm.s_scratch + "-names.lut"; }
	
	return safeExit;
}

// ********************************************************************************************************

bool readGraphFile() {

	string in_file;
	int n,m,s,f,a,b;    n = m = 0;
	elementrb *item;
	time_t t1; t1 = time(&t1);
	time_t t2; t2 = time(&t2);
	
	cout<< " HERE";
	//*** DO WE NEED THIS?**************************************************
	// First, we scan through the input file to create a list of unique node names
	// (which we store in the namesLUT), and a count of the number of edges.
	cout << ">> input file scan ( " << ioparm.f_pairs << " )" << endl;
	cout << "   edges: [0]"<<endl;
	ifstream fscan1(ioparm.f_pairs.c_str(), ios::in);
	while (fscan1 >> s >> f) {					// read friendship pair (s,f)
		if (s != f) { m++; }
		
		if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
			cout << "   edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;		// 
		}									// 
		t2=time(&t2);							// 
	}
	fscan1.close();
	cout << "   edges: ["<<m<<"]"<<endl;
	//~ d->g = new graph (ioparm.n);					// make new graph with n vertices
	//**************************************************************************
	// Finally, we reparse the file and added edges to the graph
	
	for (int j=0; j<ioparm.ng; j++) {
	
		m = 0;
		ioparm.flag_timer = false;					// reset timer
		
		cout << ">> input file read ( " << ioparm.f_pairs << " )" << endl;
		cout << "   edges: [0]"<<endl;
		
		if (!ioparm.flag_m) {in_file = ioparm.f_pairs;} else {in_file = ioparm.f_ins.at(j);}
		
		cout << ">> read graph file ( " << in_file << " )" << endl;
		ifstream fin(in_file.c_str(), ios::in);
		while (fin >> s >> f) {
			m++;
			if (s != f) {
				item = namesLUT.findItem(s); a = item->value;
				item = namesLUT.findItem(f); b = item->value;
				if (!(d->gs[j]->doesLinkExist(a,b))) { if (!(d->gs[j]->addLink(a,b))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->gs[j]->getName(a) == "") { d->gs[j]->setName(a, num2str(s)); } }
				if (!(d->gs[j]->doesLinkExist(b,a))) { if (!(d->gs[j]->addLink(b,a))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->gs[j]->getName(b) == "") { d->gs[j]->setName(b, num2str(f)); } }
			}
			if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
				cout << "   edges: ["<<m<<"]"<<endl;
				t1 = t2; ioparm.flag_timer = true;		// 
			}									// 
			t2=time(&t2);							// 
		}
		cout << ">> edges: ["<<m<<"]"<<endl;
		fin.close();
	}
	ioparm.m = d->gs[0]->numLinks();				// store actual number of directional edges created
	ioparm.n = d->gs[0]->numNodes();				// store actual number of nodes used
	cout << "vertices: ["<<ioparm.n<<"]"<<endl;
	
	return true;
}

// ********************************************************************************************************

bool readLUTFile() {
	int s,f,a;
	string st;
	elementrb *item;
	cout << ">> read names look-up table ( " << ioparm.f_lut << " )" << endl;
	ifstream fin(ioparm.f_lut.c_str(), ios::in);
	fin >> st >> st;
	while (fin >> s >> f) {
		cout << f << " " << s << endl;
		namesLUT.insertItem(f,s);		//i.e. maps real -> virtual
		namesLUTr.insertItem(s,f);
		if (s > ioparm.n) { ioparm.n = s; }
	}
	fin.close();
	ioparm.n++;
	cout << "vertices: ["<<ioparm.n<<"]"<<endl;
	
	if (!ioparm.flag_m) { ioparm.ng = 1;}
	d->ng = ioparm.ng;
	d->gs = new graph *[ioparm.ng];
	
	for (int j=0; j<ioparm.ng; j++) {
		d->gs[j] = new graph (ioparm.n);						// make new graph with n vertices
		
	//~ d->g = new graph (ioparm.n);					// make new graph with n vertices
		for (int i=0; i < ioparm.n; i++) {
			item = namesLUTr.findItem(i); a = item->value;
			if (d->gs[j]->getName(i) == "") { d->gs[j]->setName(i, num2str(a)); }
		}
	}
	return true;
}

// ********************************************************************************************************

bool readSplitsFiles() {
	string x,y,st,hists_file;
	for (int i=0; i < ioparm.f_hists.size(); i++) {
		hists_file = ioparm.f_hists.at(i);
		cout << ">> read splits file ( " << hists_file << " )" << endl;
		ifstream fin(hists_file.c_str(), ios::in);
		fin >> st >> st;
		while (fin >> x >> y) {
			//~ cout << x << " : " << y << endl;
			d->addConsensusSplit(x,atof(y.c_str()),(double)ioparm.num_hists);
		}
		fin.close();
	}
	return true;
}

// ********************************************************************************************************

void recordHRG(const unsigned int step, const double thisL) {
	
	time_t t1;
	
	// write hrg to file
	ioparm.f_dendro = ioparm.d_dir + ioparm.s_scratch + "_best-dendro.hrg";
	d->recordDendrogramStructure(ioparm.f_dendro);
	
	// write statistics about hrg to file
	ioparm.f_dendro = ioparm.d_dir + ioparm.s_scratch + "_best-dendro.info";
	
	t1 = time(&t1); 
	ofstream fout(ioparm.f_dendro.c_str(), ios::trunc);
	fout << "---HIERARCHICAL-RANDOM-GRAPH---\n";
	fout << "StartTime     : " << ioparm.start_time;
	fout << "InputFile     : " << ioparm.f_pairs << "\n";
	fout << "Directory     : " << ioparm.d_dir   << "\n";
	fout << "---Basic-Graph-Information---\n";
	fout << "Nodes         : " << ioparm.n		<< "\n";
	fout << "Edges         : " << ioparm.m/2	<< "\n";
	fout << "---HRG-Information---\n";
	fout << "OutputTime    : " << asctime(localtime(&t1));
	fout << "NumStepsMCMC  : " << step			<< "\n";
	fout << "NumPrevBests  : " << 0			<< "\n";
	fout << "LogLikelihood : " << thisL		<< "\n";
	fout << "HRG           : " << ioparm.s_scratch + "_best-dendro.hrg" << "\n";
	fout << "InfoFile      : " << ioparm.s_scratch + "_best-dendro.info"  << "\n";
	fout.close();
	
	return;
}

// ********************************************************************************************************

void	recordStats() {
	
	cout << ">> exported splits ( " << ioparm.f_splits << " )" << endl;
	d->recordSplitHistogram(ioparm.f_splits);
	
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************
