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

using namespace lzrlbwt;
using namespace std;

string out_basename=string();
string input_file=string();
ulint d = slz_rlbwt<>::DEFAULT_SKIP;

void help(){
	cout << "slz-rlbwt-build" << endl << endl;
	cout << "Usage: slz-rlbwt-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>      use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	cout << "   -d <int>           sparsification of LZ77: after each phrase, skip <int> characters before opening a new phrase. "<<endl;
	cout << 		"DEFAULT: " << slz_rlbwt<>::DEFAULT_SKIP <<endl;
	cout << "   <input_file_name>  input text file." << endl;
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

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-d")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -d option." << endl;
			help();
		}

		d = atoi(argv[ptr]);

		if(d<1){
			cout << "Error: -d accepts only arguments greater than 0." << endl;
			help();
		}

		ptr++;

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

	//parse options

    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	cout << "Building slz-rlbwt of input file '" << input_file << "'" << endl;
	cout << "Prefix '" << out_basename << "' will be used for all index files." << endl;

	string input;

	{

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		input = buffer.str();

	}

	slz_rlbwt<> idx = slz_rlbwt<>(input,true,d);
	idx.save_to_file(out_basename, true);

	//printRSSstat();

	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;

}
