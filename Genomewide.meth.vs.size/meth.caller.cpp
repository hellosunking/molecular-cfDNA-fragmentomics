#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <tr1/unordered_map>
#include <stdlib.h>

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Apr 2020
*/

const unsigned int MAX_SAM_LEN		= 4096;	// maximum length for fragment size and sam line
const unsigned int MAX_CpG_COVER	=  512;	// maximum CpG coverage of a fragment
const unsigned int MIN_ALIGN_SCORE	=	 0;	// minimum alignment score for a read to be considered
const unsigned int MIN_QUAL_SCORE	=   33;	// minimum phred score for a CpG site to be considered
const unsigned int MAX_MERGED_SEQ	=  512;

typedef struct {
	unsigned short wC;	// 'C' on watson chain
	unsigned short wT;	// 'T' on watson chain
	unsigned short wZ;	// neither 'C' nor 'T'; SNPs or sequencing errors
	unsigned short cC;	// 'C' on crick chain
	unsigned short cT;	// 'T' on crick chain
	unsigned short cZ;	// neither 'C' nor 'T'; SNPs or sequencing errors
}meth;

//typedef dense_hash_map<string, int> myhash;
// function declarations, the implementation is at the end of this file
void usage( const char * prg );
void loadgenome( const char * file, unordered_map<string, string> & genome );

void deal_PE( unordered_map<string, string> &g, const char *samfile, const char * outpre );

bool fix_cigar(string &cigar, string &realSEQ, string &realQUAL, string &seq, string &qual);
void callmeth( string &realSEQ, string &realQUAL, string &chr, int pos, bool strand,
				unordered_map<string, string> :: iterator &git, map<string, map<int, meth>*> &methcall );
void write_methcall( map<string, map<int, meth>*> & mc, map<string, int> & chrcount,
				unordered_map<string, string> & g, const char *outfile, bool TAPS );

int main( int argc, char *argv[] ) {
	if( argc != 3 ) {
		usage( argv[0] );
		return 2;
	}
	// load genome
	unordered_map<string, string> g;
	loadgenome( argv[1], g );

	char infile[128];
	char outfile[128];
	for(int i=50; i<=250; ++i ) {
		cerr << "\rProcessing " << i;
		sprintf( infile,  "%s/%d.sam", argv[2], i );
		sprintf( outfile, "%s/%d", argv[2], i );
		deal_PE( g, infile, outfile );
	}
	cerr << "\rDone. Size range 50-250 processed.\n";
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
void deal_PE( unordered_map<string, string> &g, const char *samfile, const char * outpre ) {
	unordered_map<string, string> :: iterator git;
	unordered_map<string, string> :: iterator no_such_chr = g.end();

	map<string, map<int, meth>*> methcall;
	map<string, int> chrcount;

	for( git=g.begin(); git!=g.end(); ++git ) {
		methcall.insert( pair<string, map<int, meth>*>(git->first, new map<int, meth>()) );
		chrcount.insert( pair<string, int>(git->first, 0) );
	}
	map<int, meth>* mp;
	map<int, meth> :: iterator mit;
	meth m;

	// open sam file
	ifstream fin( samfile );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << samfile << " to read!\n";
		exit(200);
	}
	//cout << "Loading alignment " << samfile << " in PE mode ...\n";

	unsigned int count = 0;
	string line1, line2, seqName, chr, cigar1, seq1, qual1, cigar2, seq2, qual2;
	register int flag;
	string mateinfo, matepos, dist, sinfo;   //fields that are ignored; all the sequence are converted to WATSON chain
	register unsigned int pos1, pos2, score;
	stringstream ss;
	string realSEQ1, realQUAL1, realSEQ2, realQUAL2;   //these are CIGAR-processed seq and qual
	string mSEQ, mQUAL; //merged sequence and quality if read1 and read2 has overlap
	line1.resize( MAX_SAM_LEN );
	line2.resize( MAX_SAM_LEN );
	mSEQ.resize( MAX_MERGED_SEQ );
	mQUAL.resize( MAX_MERGED_SEQ );
	
	int p1, p2;
	string *r1, *r2, *q1, *q2;
	bool strand;

	// load sam file
	while( true ) {
		getline( fin, line1 );
		if( fin.eof() ) break;
		getline( fin, line2 );
		//14_R1	83	chr9	73301642	42	36M	=	73301399	-279	TCCTCCTTCTCTCCCTC	HHHHHHHHH	XG:Z:CT
		//14_R2	163	chr9	73301399	42	36M	=	73301642	279	TTTATTTTGATCCTGTA	DDCBA@?>=<;986420.

		ss.clear();
		ss.str( line1 );
		ss >> seqName >> flag >> chr >> pos1 >> score >> cigar1 >> mateinfo >> matepos >> dist >> seq1 >> qual1 >> sinfo;
		//cerr << seqName << '\n';

		if( score < MIN_ALIGN_SCORE ) {
			//cerr << "Discard " << seqName << " due to poor alignment score.\n";
			continue;
		}

		// determine whether the alignemnt is on watson chain or crick chain using the XG:Z: tag
		// which is ALWAYS at the end of read1 for Msuite and Bismark
		if( line1.back() == 'T' ) {  // XG:Z:CT => watson
			strand = true;
		} else if( line1.back() == 'A' ) {	// XG:Z:GA => crick
			strand = false;
		} else {
			cerr << "Error: unknow strand! " << seqName << "\n";
			continue;
		}

		git = g.find( chr );
		if( git == no_such_chr ) continue;   // there is NO such chromosome in the genome!!!

		ss.clear();
		ss.str( line2 );
		ss >> seqName >> flag >> chr >> pos2 >> score >> cigar2 >> mateinfo >> matepos >> dist >> seq2 >> qual2;

		// process CIGAR 1, handle the indels
		realSEQ1.clear();
		realQUAL1.clear();
		if( ! fix_cigar( cigar1, realSEQ1, realQUAL1, seq1, qual1 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar1 << ") in " << seqName << "!\n";
			continue;
		}
		// process CIGAR 2, handle the indels
		realSEQ2.clear();
		realQUAL2.clear();
		if( ! fix_cigar( cigar2, realSEQ2, realQUAL2, seq2, qual2 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar2 << ") in " << seqName << "!\n";
			continue;
		}

		if( pos1 <= pos2 ) {
			p1 = pos1; r1 = & realSEQ1; q1 = & realQUAL1;
			p2 = pos2; r2 = & realSEQ2; q2 = & realQUAL2;
		} else {
			p1 = pos2; r1 = & realSEQ2; q1 = & realQUAL2;
			p2 = pos1; r2 = & realSEQ1; q2 = & realQUAL1;
		}

		if( p1 + r1->size() <= p2 ) { //there is NO overlap
			callmeth( realSEQ1, realQUAL1, chr, pos1, strand, git, methcall );
			callmeth( realSEQ2, realQUAL2, chr, pos2, strand, git, methcall );
		} else {	// there is overlap in read 1 and read 2
			//cerr << "Found overlap in " << seqName << '\n';
			if( p2+r2->size() >= p1+r1->size() ) {	// most case
				mSEQ.clear();
				mQUAL.clear();
				int len = p2 + r2->size() - p1;
				//cerr << "len\t" << len << "\tp1=" << p1 << "\tp2=" << p2 << '\n';
				register int rs = r1->size();
				register int offset = p2 - p1;
				register unsigned int k;
				for( k=0; k != offset; ++k ) {	// read1 only
					mSEQ  += r1->at(k);		// can also use substr
					mQUAL += q1->at(k);
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				for( ; k != rs; ++k ) {	// overlapped region, peak the one with higher quality
					if( q1->at(k) >= q2->at(k-offset) ) {
						mSEQ  += r1->at(k);
						mQUAL += q1->at(k);
					} else {
						mSEQ  += r2->at(k-offset);
						mQUAL += q2->at(k-offset);
					}
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				for( ; k != len; ++k ) {	//read2 only
					mSEQ  += r2->at(k-offset);
					mQUAL += q2->at(k-offset);
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";
				callmeth( mSEQ, mQUAL, chr, p1, strand, git, methcall );
			} else {	// rare case that R1 completely contains R2 => use R1 directly
				callmeth( realSEQ1, realQUAL1, chr, pos1, strand, git, methcall );
			}
		}

		chrcount.find( chr )->second ++;
		++ count;
//		if( ! (count & 0x003fffff) )
//			cout << '\r' << count << " lines loaded.";
	}
	fin.close();
//	cout << '\r' << "Done: " << count << " lines loaded.\n";

	write_methcall( methcall, chrcount, g, outpre, false );
}

/*
 * function implementations
*/
void loadgenome( const char * file, unordered_map<string, string> & genome ) {
	ifstream fin( file );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << file << " !\n";
		exit(200);
	}
	cout << "Loading genome: " << file << '\n';
	string line, chr, tmp="X";	// X is for position-taking
	register unsigned int i, j;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;
		j = line.length();
		if( line[0] == '>' ) {
			chr.clear();
			for( i=1; i!=j; ++i ) {
				if( line[i]==' ' || line[i]=='\t' )break;
				chr += line[i];
			}
			genome.insert( pair<string, string>(chr, tmp) );
//			cout << "=>Adding " << chr << " ...\n";
		} else {
			/*
			for( i=0; i!=j; ++i ) {
				if( line[i] >= 'a' )	// lowercase to uppercase
					line[i] -= 'a' - 'A';
			}*/
			genome[chr] += line;
		}
	}
	fin.close();
}

// update sequence and quality using CIGAR information to support indels
bool fix_cigar( string &cigar, string &realSEQ, string &realQUAL, string &seq, string &qual ) {
	register int i, j, k;
	j = 0;
	register int curr = 0;
	register int cs = cigar.size();
	realSEQ.clear();
	realQUAL.clear();
	for(i=0; i!=cs; ++i) {
		if( cigar[i] <= '9' ) {   // digital
			j *= 10;
			j += cigar[i] - '0';
		} else {	// MUST be M, I, or D
			if( cigar[i] == 'M' ) { // match or mismatch, copy seq and qual
				for(k=0; k!=j; ++k) {
					realSEQ  +=  seq[ curr+k ];
					realQUAL += qual[ curr+k ];
				}
				curr += j;
			} else if ( cigar[i] == 'I' ) { // insertion, discard this part
				curr += j;
			} else if ( cigar[i] == 'D' ) { // deletion, add place holders
				for(k=0; k!=j; ++k) {
					realSEQ  +=  'N';
					realQUAL += '\0';
				}
			} else {	// unsupported CIGAR element
				return false;
			}
			j = 0;
		}
	}
	return true;
}

// call meth from sequence
void callmeth( string &realSEQ, string &realQUAL, string &chr, int pos, bool strand,
				unordered_map<string, string> :: iterator &git, map<string, map<int, meth>*> &methcall ) {
	map<int, meth>* mp;
	map<int, meth> :: iterator mit;
	meth m;
	unsigned int rs = realSEQ.size();
	char c1, c2;
	unsigned int i = 0;
	unsigned int j = pos + i;
	for( ; i!=rs; ++i, ++j) {
		c1 = git->second[j];
		c2 = git->second[j+1];
		if( (c1=='C' || c1=='c') && (c2=='G' || c2=='g') ) {  // meet a CpG site
			if( realQUAL[i] < MIN_QUAL_SCORE ) continue;
//			cerr << "Found CpG: " << pos+i << "(" << i << "), sequenced " << realSEQ[i] << '\n';

			mp  = methcall.find( chr )->second;
			mit = mp->find( j );
			if( mit == mp->end() ) {	// no record, insert one
				//cerr << "New site!\n";
				if( strand ) {	// watson chain, check 'C' site
					if( realSEQ[i] == 'C' ) {
						m.wC=1; m.wT=0; m.wZ=0; m.cC=0; m.cT=0; m.cZ=0;
					} else if( realSEQ[i] == 'T' ) {
						m.wC=0; m.wT=1; m.wZ=0; m.cC=0; m.cT=0; m.cZ=0;
					} else {
						m.wC=0; m.wT=0; m.wZ=1; m.cC=0; m.cT=0; m.cZ=0;
					}
				} else {	// crick chain, check 'G' site
					if( i+1 == rs ) break;	// "C" is at the end of the read, could not infer methylation level
					if( realSEQ[i+1] == 'G' ) {
						m.cC=1; m.cT=0; m.cZ=0; m.wC=0; m.wT=0; m.wZ=0;
					} else if( realSEQ[i+1] == 'A' ) {
						m.cC=0; m.cT=1; m.cZ=0; m.wC=0; m.wT=0; m.wZ=0;
					} else {
						m.cC=0; m.cT=0; m.cZ=1; m.wC=0; m.wT=0; m.wZ=0;
					}
				}
				mp->insert( pair<int, meth>(j, m) );
			} else {	// there is such record
				//cerr << "Existing site!\n";
				if( strand ) {
					if( realSEQ[i] == 'C' ) {
						mit->second.wC ++;
					} else if( realSEQ[i] == 'T' ) {
						mit->second.wT ++;
					} else {
						mit->second.wZ ++;
					}
				} else {
					if( i+1 == rs ) break;
					if( realSEQ[i+1] == 'G' ) {
						mit->second.cC ++;
					} else if( realSEQ[i+1] == 'A' ) {
						mit->second.cT ++;
					} else {
						mit->second.cZ ++;
					}
				}
			}
		}
	}
}

// write meth call into file
void write_methcall( map<string, map<int, meth>*> & mc, map<string, int> & chrcount,
						unordered_map<string, string> & g, const char *outpre, bool TAPS ) {

	string outfile = outpre;
	outfile += ".meth.call";
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "ERROR: write output meth call failed.\n";
		exit(20);
	}

	outfile = outpre;
	outfile += ".meth.log";
	ofstream flog( outfile.c_str() );
	if( flog.fail() ) {
		cerr << "ERROR: write output failed.\n";
		exit(21);
	}

	unordered_map<string, string> :: iterator git;
	map<string, map<int, meth>*> :: iterator mit;	// methcall iterator
	map<int, meth> :: iterator cit;	// methcall iterator for each chromosome
	map<string, int> :: iterator chrit;

	meth m;

	fout << "#chr\tLocus\tTotal\twC\twT\twOther\tContext\tcC\tcT\tcOther\n";
	flog << "#chr\tReads\tTotal.wC\tTotal.wT\tTotal.cC\tTotal.cT\n";
	for( mit=mc.begin(); mit!=mc.end(); ++mit ) {
		git = g.find( mit->first );
		int totalWC=0, totalWT=0, totalCC=0, totalCT=0;	//total C, T
		for( cit=mit->second->begin(); cit!=mit->second->end(); ++cit ) {
			int i = cit->first;
			m = cit->second;
			unsigned int Valid = m.wC+m.wT+m.cC+m.cT;

			fout << mit->first << '\t' << i << '\t' << Valid+m.wZ+m.cZ << '\t'
				 << m.wC << '\t' << m.wT << '\t' << m.wZ << '\t'
				 << git->second[i-1] << git->second[i] << git->second[i+1] << git->second[i+2] << '\t'
				 << m.cC << '\t' << m.cT << '\t' << m.cZ << '\n';

			totalWC += m.wC;
			totalWT += m.wT;
			totalCC += m.cC;
			totalCT += m.cT;
		}

		chrit = chrcount.find( mit->first );
		flog << mit->first << '\t' << chrit->second << '\t'
			 << totalWC << '\t' << totalWT << '\t' << totalCC << '\t' << totalCT << '\n';

		delete mit->second;
	}
	fout.close();
	flog.close();
}

// show usage
void usage( const char * prg ) {
	cerr << "\nUsage: " << prg << " <genome.fa> <in.DIR.name after split>\n\n";
}

