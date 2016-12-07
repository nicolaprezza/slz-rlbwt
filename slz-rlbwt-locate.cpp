
/*
 *  This file is part of slz-rlbwt.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>,
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, and Mathieu Raffinot
 *
 *   slz-rlbwt is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   slz-rlbwt is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

#include <iostream>
#include <slz-rlbwt.h>
#include <utils.h>

string out_file = "";
bool opt = true;

using namespace lzrlbwt;
using namespace std;

void help(){
	cout << "slz-rlbwt-locate: locate all occurrences of the input patterns. Note that this program" << endl;
	cout << "discards the output (i.e. text positions), and should be used only for benchmark purposes." << endl << endl;
	cout << "Usage: slz-rlbwt-locate [options] <index_basename> <patterns_file>" << endl;
	cout << "   -o	<file>          Stream occurrences to path. If specified, occurrences are not stored in RAM." << endl;
	cout <<	"                       If not specified, occurrences are stored in RAM." << endl;
	cout << "   <index_basename>    basename of all index files" << endl;
	cout << "   <patterns_file>     file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_file = string(argv[ptr]);
		ptr++;

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}


void search(string idx_basename, string patterns, ostream& ostr){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

    slz_rlbwt<> idx;

	idx.load_from_file(idx_basename,true);

	auto t2 = high_resolution_clock::now();

	cout << "searching patterns ... " << endl;
	ifstream ifs(patterns);

	//read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
	string header;
	std::getline(ifs, header);

	ulint n = get_number_of_patterns(header);
	ulint m = get_patterns_length(header);

	uint last_perc = 0;

	ulint occ_tot=0;
	ulint occ_tot_count=0;

	//extract patterns from file and search them in the index
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}

		string p = string();

		for(ulint j=0;j<m;++j){
			char c;
			ifs.get(c);
			p+=c;
		}

		//cout << "locating " << idx.count(p) << " occurrences of "<< p << " ... " << flush;

		//occ_tot_count += idx.count(p);
		occ_tot += idx.locate(p,ostr);	//stream occurrences to ostr, return number of occurrences
		//occ_tot += idx.locate(p).size();	//stream occurrences to ostr, return number of occurrences

	}

	double occ_avg = (double)occ_tot / n;

	cout << endl << occ_avg << " average occurrences per pattern" << endl;

	ifs.close();

	auto t3 = high_resolution_clock::now();

	//printRSSstat();

	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	cout << "number of patterns n = " << n << endl;
	cout << "pattern length m = " << m << endl;
	cout << "total number of occurrences  occ_t = " << occ_tot << endl;
	//cout << "total number of counted occurrences = " << occ_tot_count << endl;
	cout << "m * occ_t  = " << occ_tot*m << endl;
	cout << "n*m + occ_t  = " << n*m+occ_tot << endl << endl;

	cout << "Total time : " << search << " milliseconds" << endl;
	cout << "Search time : " << (double)search/n << " milliseconds/pattern (total: " << n << " patterns)" << endl;
	cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;

}

int main(int argc, char** argv){

	int ptr = 1;

	if(argc<3) help();

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	string idx_basename, patterns;

	idx_basename = string(argv[ptr++]);
	patterns = string(argv[ptr++]);

	cout << "Loading slz-rlbwt index" << endl;

	if(out_file.compare("")==0){

		stringstream sstr;
		search(idx_basename,patterns,sstr);

	}else{

		ofstream ofs(out_file);
		search(idx_basename,patterns,ofs);
		ofs.close();

	}

}
